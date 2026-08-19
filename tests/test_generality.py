"""Gate 3a: the same code, a structurally different experiment.

Every other test in this suite is shaped by the RTgill-W1 design - four days, nine
concentrations, `Metadata_*` columns, tab-separated. That leaves one claim untested, and
it is the claim the whole package rests on: that the design comes from the configuration
and from nothing else.

`reshape.py` builds a table that shares its *values* with the published one and none of
its *structure*: no time axis, five exposure levels instead of nine, level labels that do
not sort into the response order, different column names on both sides of the split,
commas instead of tabs, interleaved columns and shuffled rows. This module runs the whole
chain on it and checks the answers are the ones the configuration asks for.

Two bugs came out of writing it, and both are regression-tested below:

* An experiment with no `[schema] group_by` has a single stratum labelled `""`. Written to
  a delimited file that is an empty field, and pandas reads an empty field back as NaN -
  which `groupby` then drops, silently, taking every row with it. `fit` produced an empty
  table, `select` reported "0 of 0", and nothing anywhere said why.
* `prune` on an empty feature list raised `IndexError` out of the clustering, because zero
  features had been treated as one empty cluster rather than as no clusters at all.

What this does **not** show is biological generality. It is one dataset in two shapes, not
two datasets.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
import reshape
from cases import SUBSET, config_text
from reshape import CONTROL, LEVELS

from cdr_fs.cli import main
from cdr_fs.config import ConfigError, load_config
from cdr_fs.schema import read_stage_table

#: 20 features x 5 exposure levels, one stratum.
FEATURES = 20
#: The reshaped metadata: exposure, batch, sample, object_uid and three qc_ columns.
METADATA = 7


@pytest.fixture(scope="module")
def chain(tmp_path_factory):
    """Run every stage once on the reshaped table, and hold on to what it produced."""
    directory = tmp_path_factory.mktemp("reshaped")
    config = reshape.build(directory)
    for stage in ("check", "emd", "fit", "select", "prune", "subset", "plot"):
        assert main([stage, "-c", str(config)]) == 0, stage
    return SimpleNamespace(
        config=config,
        loaded=load_config(config),
        results=directory / "results",
    )


def read(chain, name: str) -> pd.DataFrame:
    return read_stage_table(chain.results / name, ",")


def features_in(path) -> list[str]:
    return path.read_text(encoding="utf-8").split()


# --------------------------------------------------------------- the schema is read, not assumed


def test_the_reshaped_table_shares_no_column_name_with_the_published_one():
    """If it did, a leftover hardcoded name could still be doing the work."""
    published = set(pd.read_csv(reshape.SUBSET, sep="\t", nrows=0).columns)
    reshaped = set(reshape.reshape(pd.read_csv(reshape.SUBSET, sep="\t", nrows=8)).columns)
    assert published & reshaped == set()


def test_check_resolves_the_reshaped_schema(tmp_path, capsys):
    config = reshape.build(tmp_path)
    assert main(["check", "-c", str(config), "--scan"]) == 0
    printed = capsys.readouterr().out
    assert f"{METADATA + FEATURES} = {METADATA} metadata + {FEATURES} feature(s)" in printed
    # Two families, so the breakdown says something rather than naming one prefix.
    assert "alpha_* 10, beta_* 10" in printed
    assert "group_by          (none - a single stratum)" in printed
    assert "configuration is valid" in printed


def test_the_output_files_follow_the_configured_separator(chain):
    """`sep = comma` has to reach the written names, not only the reader."""
    produced = sorted(path.name for path in chain.results.iterdir() if path.suffix != ".png")
    assert produced == [
        "emd.csv",
        "emd_baseline.csv",
        "fit.csv",
        "prune_clusters.csv",
        "prune_linkage.csv",
        "pruned.txt",
        "select_evidence.csv",
        "selected.txt",
        "subset_pruned.csv",
        "subset_pruned_features.csv",
        "subset_pruned_retained.txt",
    ]


# ------------------------------------------------------------------------- the design is followed


def test_the_contrast_set_is_the_five_declared_levels(chain):
    distances = read(chain, "emd.csv")
    assert len(distances) == FEATURES * len(LEVELS)
    assert list(dict.fromkeys(distances["contrast"])) == [
        f"{CONTROL}v{level}" for level in LEVELS
    ]


def test_the_exposure_axis_follows_the_configured_order_not_the_alphabet(chain):
    """The sharp one. Sorted, the labels are high, low, max, mid, trace.

    `[design] levels` says trace, low, mid, high, max, and the fitted axis has to be that.
    A tool that sorted its levels would still produce a full series here, and a plausible
    one - it would simply be answering a different question.
    """
    from cdr_fs.fit import axis_positions, series_from_emd

    assert axis_positions(chain.loaded) == {level: float(i) for i, level in enumerate(LEVELS)}
    assert sorted(LEVELS) != list(LEVELS)  # otherwise this test proves nothing

    distances = read(chain, "emd.csv")
    by_level = {
        (row.feature, row.contrast.removeprefix(f"{CONTROL}v")): row.emd
        for row in distances.itertuples()
    }
    series = {feature: y for feature, _, _, y, _ in series_from_emd(chain.loaded, distances)}
    assert series
    for feature, y in series.items():
        assert list(y) == [by_level[(feature, level)] for level in LEVELS]


def test_every_replicate_pair_contributes_to_the_baseline(chain):
    """Four batches, so six unordered pairs - and one stratum, not four."""
    baseline = read(chain, "emd_baseline.csv")
    assert len(set(baseline["contrast"])) == 6
    assert set(baseline["stratum"]) == {""}


def test_the_fit_uses_the_configured_models_on_five_points(chain):
    fits = read(chain, "fit.csv")
    assert set(fits["model"]) <= set(chain.loaded.fit.models)
    assert "BC5" not in set(fits["model"])
    assert set(fits["n_points"]) == {len(LEVELS)}
    assert set(fits["stratum"]) == {""}


def test_five_levels_will_not_carry_a_five_parameter_model(tmp_path):
    """The model list is a structural fact, and the configuration is made to carry it.

    BC5 has five free parameters and this design fits five points, so the comparison is not
    identifiable. Refusing it up front is the difference between a shorter series and a
    table of meaningless information criteria.
    """
    config = reshape.build(tmp_path, {"fit.models": "BC4,BC5,LL4,WB1.4,Lin,Con"})
    with pytest.raises(ConfigError, match="BC5, which has 5 free parameters"):
        load_config(config)


def test_the_chain_narrows_the_feature_list_at_every_stage(chain):
    selected = features_in(chain.results / "selected.txt")
    pruned = features_in(chain.results / "pruned.txt")
    retained = features_in(chain.results / "subset_pruned_retained.txt")
    assert 0 < len(retained) <= len(pruned) <= len(selected) <= FEATURES
    assert set(retained) <= set(pruned) <= set(selected)


def test_the_subset_carries_the_reshaped_metadata_and_the_retained_features(chain):
    subset = pd.read_csv(chain.results / "subset_pruned.csv", nrows=1)
    retained = features_in(chain.results / "subset_pruned_retained.txt")
    metadata = [column for column in subset.columns if column not in set(retained)]
    assert len(metadata) == METADATA
    # Metadata first, then the features, whatever order the input interleaved them in.
    assert list(subset.columns) == metadata + retained


def test_the_figures_are_named_for_a_single_unnamed_stratum(chain):
    """`fit_<stratum>_part_<n>.png` has no stratum to name here, and must not invent one."""
    drawn = sorted(path.name for path in chain.results.glob("*.png"))
    assert "fit_part_1.png" in drawn
    assert not any(name.startswith("fit_nan") or name.startswith("fit__") for name in drawn)
    assert {"emd.png", "emd_baseline.png", "dendrogram.png"} <= set(drawn)
    assert all((chain.results / name).stat().st_size > 0 for name in drawn)


# ------------------------------------------------------------------------------- the two bugs


def test_a_single_unnamed_stratum_survives_the_round_trip(tmp_path):
    """An empty field is a label here, not a missing value.

    Read with a bare `pandas.read_csv` it comes back as NaN, and `groupby` drops NaN keys -
    so the table empties itself and only for the unstratified design.
    """
    path = tmp_path / "emd.csv"
    path.write_text("feature,stratum,emd\nalpha_01,,1.5\nbeta_01,,2.0\n", encoding="utf-8")

    assert pd.read_csv(path, dtype={"stratum": str})["stratum"].isna().all()

    table = read_stage_table(path, ",")
    assert list(table["stratum"]) == ["", ""]
    assert len(table.groupby(["feature", "stratum"])) == 2


def test_the_unstratified_chain_actually_fits_something(chain):
    """The end-to-end form of the same bug: `fit` used to write a table with no rows."""
    assert len(read(chain, "fit.csv")) > 0
    assert len(read(chain, "select_evidence.csv")) == FEATURES


def test_no_features_means_no_clusters():
    """Zero features used to come back as one empty cluster, which has no member to keep."""
    from cdr_fs.prune import cluster_features

    groups, tree = cluster_features(np.empty((0, 0)), [], cut=0.1, method="average")
    assert groups == []
    assert tree.shape == (0, 4)

    groups, tree = cluster_features(np.zeros((1, 1)), ["only"], cut=0.1, method="average")
    assert groups == [(0,)]
    assert tree.shape == (0, 4)


def test_pruning_an_empty_selection_is_refused_rather_than_crashing(tmp_path, capsys):
    config = reshape.build(tmp_path)
    results = tmp_path / "results"
    results.mkdir()
    (results / "selected.txt").write_text("", encoding="utf-8")

    assert main(["prune", "-c", str(config)]) == 3
    assert "selected.txt is empty" in capsys.readouterr().err


def test_subsetting_an_empty_list_is_refused_rather_than_writing_metadata(tmp_path, capsys):
    config = reshape.build(tmp_path)
    listing = tmp_path / "none.txt"
    listing.write_text("\n\n", encoding="utf-8")

    assert main(["subset", "-c", str(config), "--features", str(listing)]) == 3
    assert "none.txt is empty" in capsys.readouterr().err
    assert not (tmp_path / "results" / "subset_none.csv").exists()


def test_fitting_nothing_is_an_error_not_an_empty_table(tmp_path, capsys):
    """`select` on an empty fit table would report "0 of 0", which reads as a result."""
    config = reshape.build(tmp_path)
    results = tmp_path / "results"
    results.mkdir()
    # A distance table whose only contrast is not one the design declares, so no series is
    # complete and nothing can be fitted.
    (results / "emd.csv").write_text(
        "feature,stratum,contrast,group_a,group_b,replicate,n_a,n_b,emd\n"
        "alpha_01,,vehiclevtrace,vehicle,trace,,10,10,\n",
        encoding="utf-8",
    )
    assert main(["fit", "-c", str(config)]) == 3
    assert "no series was fitted" in capsys.readouterr().err
    assert not (results / "fit.csv").exists()


# ------------------------------------------------------------- the reshaped config is a real one


#: The configuration keys that carry the structure of the experiment. Every one of them
#: must differ from the published design, or this module is testing that design again
#: under a different name. `[trim] enabled` and `[prune] enabled` are deliberately absent:
#: both runs trim and prune, and that is a choice about method rather than about shape.
STRUCTURAL_KEYS = (
    "input.sep",
    "schema.metadata_patterns",
    "schema.condition",
    "schema.group_by",
    "schema.pool_over",
    "design.control",
    "design.levels",
    "design.exclude_from_fit",
    "trim.lower_percentile",
    "trim.upper_percentile",
    "trim.scope",
    "fit.models",
    "subset.drop_missing",
)


@pytest.mark.parametrize("dotted", STRUCTURAL_KEYS)
def test_the_reshaped_configuration_differs_where_it_has_to(dotted):
    """A value inherited from `cases.BASE` would make the whole module vacuous."""
    from cases import BASE

    section, _, key = dotted.partition(".")
    assert reshape.RESHAPED[section].get(key) != BASE[section].get(key)


def test_the_reshaped_table_is_reproducible(tmp_path):
    """Seeded, so a failure can be looked at rather than re-rolled."""
    first = reshape.build(tmp_path / "a").parent / "reshaped.csv"
    second = reshape.build(tmp_path / "b").parent / "reshaped.csv"
    assert first.read_bytes() == second.read_bytes()
    assert config_text(base=reshape.RESHAPED, table="t", output="o").startswith("[input]")
