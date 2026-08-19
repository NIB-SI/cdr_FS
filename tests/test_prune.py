"""Correlation pruning.

The distance matrices here are built from feature values chosen so that the correlations
are exact - a copy, a sign flip, a shift - which is what lets a threshold test assert on a
boundary rather than on an approximation.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from cases import write_config

from cdr_fs.config import ConfigError, load_config
from cdr_fs.prune import (
    CLUSTER_COLUMNS,
    LINKAGE_COLUMNS,
    aggregate_units,
    cluster_features,
    correlation_distance,
    prune_features,
)

#: Eight values with no pattern, used as the base series every case is derived from.
BASE = np.array([1.0, 4.0, 2.0, 8.0, 5.0, 7.0, 3.0, 6.0])


def units(**columns) -> pd.DataFrame:
    """A unit x feature matrix, one keyword per feature."""
    return pd.DataFrame(columns)


def config_for(tmp_path, **overrides):
    return load_config(write_config(tmp_path, overrides))


def frame_with(values: dict[str, list[float]], *, unit: list[str]) -> pd.DataFrame:
    """An object-level table: one metadata column naming the unit, plus features."""
    return pd.DataFrame({"Metadata_Well": unit, **values})


# ------------------------------------------------------------------- what clusters what


def test_a_copy_and_a_sign_flip_are_both_redundancy():
    """`1 - |r|`, so r = -1 is as redundant as r = +1. That is the point of the absolute value."""
    matrix = units(
        a=BASE,
        b=BASE * 3 + 7,  # r = +1
        c=-BASE,  # r = -1
        d=np.roll(BASE, 3),  # unrelated
    )
    distance, undefined = correlation_distance(matrix)
    assert undefined == 0
    groups, _ = cluster_features(distance, list(matrix.columns), cut=0.1)
    assert sorted(len(group) for group in groups) == [1, 3]
    assert set(groups[0]) == {0, 1, 2}


@pytest.mark.parametrize(
    "cut, merges",
    [
        (0.1, False),  # exactly at the distance: stays apart
        (np.nextafter(0.1, 1.0), True),  # one ulp above: merges
        (0.05, False),
        (0.2, True),
    ],
)
def test_the_cut_is_exclusive_at_the_boundary(cut, merges):
    """Merging happens strictly below the cut, so a pair sitting exactly on it stays apart.

    scikit-learn's `distance_threshold` is the distance at or above which clusters are not
    merged, so this is the published rule. The matrix is written by hand rather than derived
    from a correlation, because the point is bit-exact behaviour at one value and a
    constructed `r = 0.9` is only approximately 0.9.
    """
    distance = np.array([[0.0, 0.1], [0.1, 0.0]])
    groups, _ = cluster_features(distance, ["a", "b"], cut=cut)
    assert (len(groups) == 1) is merges


def test_a_threshold_on_r_becomes_a_cut_on_distance(tmp_path):
    """`[prune] threshold` is on |r| and the cut is `1 - threshold`; nothing else converts."""
    config = config_for(tmp_path, **{"prune.threshold": "0.75"})
    assert 1.0 - config.prune.threshold == 0.25


def test_an_incomputable_correlation_clusters_with_nothing():
    """A constant feature has no correlation with anything, and is kept rather than dropped."""
    matrix = units(a=BASE, b=BASE * 2, flat=np.full(len(BASE), 3.0))
    distance, undefined = correlation_distance(matrix)
    assert undefined == 2
    assert distance[2, 0] == 1.0 and distance[2, 1] == 1.0
    # The diagonal is still zero: a feature is not at distance 1 from itself.
    assert np.diag(distance).tolist() == [0.0, 0.0, 0.0]
    groups, _ = cluster_features(distance, list(matrix.columns), cut=0.1)
    assert sorted(len(group) for group in groups) == [1, 2]


def test_one_feature_is_one_cluster():
    """The n = 1 edge case: scipy has no condensed distance form for a single point."""
    distance, _ = correlation_distance(units(only=BASE))
    groups, tree = cluster_features(distance, ["only"], cut=0.1)
    assert groups == [(0,)]
    assert tree.shape == (0, 4)


def test_cluster_numbering_does_not_depend_on_column_order():
    """Reordering the input must give the same partition, or the run is not reproducible."""
    values = {"a": BASE, "b": BASE * 2, "c": np.roll(BASE, 3), "d": np.roll(BASE, 3) * -4}
    forward = units(**values)
    reversed_ = forward[list(reversed(forward.columns))]

    def partition(matrix):
        distance, _ = correlation_distance(matrix)
        groups, _ = cluster_features(distance, list(matrix.columns), cut=0.1)
        return {frozenset(matrix.columns[i] for i in group) for group in groups}

    assert partition(forward) == partition(reversed_)


# ------------------------------------------------------------------------- aggregation


def test_aggregation_takes_the_median_per_unit():
    frame = frame_with({"f": [1.0, 2.0, 30.0, 4.0, 5.0, 6.0]}, unit=list("AAABBB"))
    aggregated = aggregate_units(frame, ["f"], ["Metadata_Well"], fill_missing="none")
    assert aggregated["f"].tolist() == [2.0, 5.0]


def test_fill_missing_moves_unit_medians_toward_the_global_mean():
    """`column_mean` is the published choice, and it is not a no-op.

    Unit A holds 1 and 2 with one value trimmed away. `none` medians the two that are left;
    `column_mean` first substitutes the overall mean, which drags A's median up.
    """
    frame = frame_with(
        {"f": [1.0, 2.0, np.nan, 10.0, 11.0, 12.0]}, unit=list("AAABBB")
    )
    without = aggregate_units(frame, ["f"], ["Metadata_Well"], fill_missing="none")
    with_fill = aggregate_units(frame, ["f"], ["Metadata_Well"], fill_missing="column_mean")
    assert without["f"].tolist() == [1.5, 11.0]
    assert with_fill["f"].tolist() == [2.0, 11.0]  # mean of the five present values is 7.2


def test_no_aggregation_columns_correlates_the_rows_as_they_are():
    frame = frame_with({"f": [1.0, 2.0, 3.0, 4.0]}, unit=list("AABB"))
    assert aggregate_units(frame, ["f"], []).equals(frame[["f"]])


def test_aggregation_rejects_unknown_columns():
    frame = frame_with({"f": [1.0, 2.0]}, unit=["A", "B"])
    with pytest.raises(KeyError, match="Metadata_Day"):
        aggregate_units(frame, ["f"], ["Metadata_Day"])


# ------------------------------------------------------------------------ the linkage


def test_the_linkage_table_carries_the_tree_and_its_leaf_order():
    """`plots.py` redraws the dendrogram from this table alone, so it must be complete."""
    matrix = units(a=BASE, b=BASE * 2, c=np.roll(BASE, 3), d=-BASE)
    distance, _ = correlation_distance(matrix)
    _, tree = cluster_features(distance, list(matrix.columns), cut=0.1)

    from cdr_fs.prune import linkage_table

    table = linkage_table(tree, list(matrix.columns))
    assert list(table.columns) == list(LINKAGE_COLUMNS)
    leaves = table[table["label"].notna()]
    merges = table[table["label"].isna()]
    assert leaves["label"].tolist() == ["a", "b", "c", "d"]
    assert len(merges) == 3
    rebuilt = merges[["left", "right", "height", "size"]].to_numpy(dtype=float)
    assert np.allclose(rebuilt, tree)


# ------------------------------------------------------------------------- the stage


def test_the_stage_keeps_one_member_per_cluster(tmp_path):
    frame = frame_with(
        {
            "zebra": list(BASE) + list(BASE),
            "apple": list(BASE * 2 + 1) + list(BASE * 2 + 1),
            "other": list(np.roll(BASE, 3)) + list(np.roll(BASE, 3)),
        },
        unit=[f"W{index:02d}" for index in range(8)] * 2,
    )
    config = config_for(tmp_path, **{"prune.aggregate_by": "Metadata_Well"})
    features = ["zebra", "apple", "other"]
    kept, clusters, _, report = prune_features(config, frame, features)

    # Alphabetical inside the cluster, but the returned list is in the input's own order.
    assert kept == ["apple", "other"]
    assert list(clusters.columns) == list(CLUSTER_COLUMNS)
    assert clusters["cluster"].tolist() == [1, 1, 2]
    assert clusters.loc[clusters["representative"], "feature"].tolist() == ["apple", "other"]
    assert clusters["distance_to_representative"].max() == pytest.approx(0.0, abs=1e-12)
    assert (report.features, report.kept, report.removed) == (3, 2, 1)
    assert report.singletons == 1 and report.largest == 2
    assert "pruned 3 to 2 feature(s)" in report.summary()


def test_representative_first_keeps_the_earliest_column(tmp_path):
    frame = frame_with(
        {"zebra": list(BASE) * 2, "apple": list(BASE * 2) * 2},
        unit=[f"W{index:02d}" for index in range(8)] * 2,
    )
    config = config_for(
        tmp_path,
        **{"prune.aggregate_by": "Metadata_Well", "prune.representative": "first"},
    )
    kept, _, _, _ = prune_features(config, frame, ["zebra", "apple"])
    assert kept == ["zebra"]


def test_the_report_names_features_with_no_data(tmp_path):
    """An empty feature survives pruning, so the report has to say it is there."""
    frame = frame_with(
        {"good": list(BASE) * 2, "empty": [np.nan] * 16},
        unit=[f"W{index:02d}" for index in range(8)] * 2,
    )
    config = config_for(tmp_path, **{"prune.aggregate_by": "Metadata_Well"})
    kept, _, _, report = prune_features(config, frame, ["good", "empty"])
    assert kept == ["good", "empty"]
    assert report.features_all_missing == ("empty",)
    assert "entirely missing" in report.summary()


# ---------------------------------------------------------------------- configuration


def test_aggregate_by_defaults_to_the_trim_scope(tmp_path):
    config = config_for(tmp_path)
    assert config.prune.aggregate_by == config.trim.scope
    assert "prune.aggregate_by" in config.defaulted


def test_aggregate_by_can_be_set_to_nothing(tmp_path):
    """Explicitly empty means "correlate the rows as they are", which is not the default."""
    config = config_for(tmp_path, **{"prune.aggregate_by": ""})
    assert config.prune.aggregate_by == ()
    assert "prune.aggregate_by" not in config.defaulted


def test_pruning_without_a_unit_to_aggregate_over_is_rejected(tmp_path):
    with pytest.raises(ConfigError, match=r"\[prune\] aggregate_by"):
        config_for(tmp_path, **{"trim.enabled": "false", "trim.scope": None})


def test_the_published_configuration_reads_as_expected(tmp_path):
    config = config_for(tmp_path)
    assert config.prune.threshold == 0.9
    assert config.prune.linkage == "average"
    assert config.prune.representative == "alphabetical"
    assert config.prune.fill_missing == "column_mean"
