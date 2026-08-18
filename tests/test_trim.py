"""Extreme-value trimming.

The load-bearing test is `test_matches_the_original_per_group_loop`: the published
pipeline trimmed with a Python loop calling `np.nanpercentile` per group per feature, and
the vectorized implementation must agree with it value for value. Everything else pins the
semantics that loop had by accident of how it was written, and which the rewrite now
states on purpose.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from cdr_fs.trim import trim_extremes

SCOPE = ["day", "well"]
FEATURES = ["f1", "f2", "f3"]


def sample(rows: int = 600, seed: int = 7) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    frame = pd.DataFrame(
        {
            "day": rng.choice(["D1", "D5", "D7"], rows),
            "well": rng.choice(["A01", "A02", "B01", "B02"], rows),
            "f1": rng.normal(size=rows),
            "f2": rng.lognormal(size=rows),
            "f3": rng.normal(scale=50, size=rows),
        }
    )
    frame.loc[rng.choice(rows, 30, replace=False), "f1"] = np.nan
    return frame


def original_loop(frame, features, scope, lower, upper):
    """The trimming loop of `only_trimming_well.py`, as the oracle for this rewrite."""

    def trim_extremes_original(values, lower_percentile, upper_percentile):
        if len(values) == 0:
            return values
        low = np.nanpercentile(values, lower_percentile)
        high = np.nanpercentile(values, upper_percentile)
        return values[(values >= low) & (values <= high)]

    parts = []
    for _, group in frame.groupby(scope, sort=False):
        trimmed = group.copy()
        for feature in features:
            trimmed[feature] = trim_extremes_original(
                group[feature].dropna(), lower, upper
            )
        parts.append(trimmed)
    return pd.concat(parts).loc[frame.index]


@pytest.mark.parametrize(("lower", "upper"), [(2.5, 97.5), (0.0, 100.0), (10.0, 90.0)])
def test_matches_the_original_per_group_loop(lower, upper):
    frame = sample()
    trimmed, _ = trim_extremes(frame, FEATURES, SCOPE, lower, upper)
    expected = original_loop(frame, FEATURES, SCOPE, lower, upper)
    for feature in FEATURES:
        assert np.array_equal(
            trimmed[feature].to_numpy(), expected[feature].to_numpy(), equal_nan=True
        )


def test_removes_values_not_rows():
    frame = sample()
    trimmed, _ = trim_extremes(frame, FEATURES, SCOPE, 2.5, 97.5)
    assert len(trimmed) == len(frame)
    # An object trimmed on one feature keeps its other features.
    dropped = trimmed["f1"].isna() & frame["f1"].notna()
    assert dropped.any()
    assert trimmed.loc[dropped, "f3"].notna().any()


def test_metadata_is_untouched():
    frame = sample()
    trimmed, _ = trim_extremes(frame, FEATURES, SCOPE, 2.5, 97.5)
    for column in SCOPE:
        assert trimmed[column].equals(frame[column])


def test_interval_is_closed_so_a_group_is_never_emptied():
    # The values sitting exactly on the percentile bounds survive, which is why the
    # report's empty cells are always cells that arrived empty.
    frame = pd.DataFrame({"day": ["D1"] * 5, "well": ["A01"] * 5, "f1": [1.0, 2, 3, 4, 5]})
    trimmed, report = trim_extremes(frame, ["f1"], SCOPE, 2.5, 97.5)
    assert trimmed["f1"].notna().sum() >= 1
    assert report.empty_cells == 0


def test_trims_within_each_group_independently():
    # One wild group must not drag another group's bounds with it.
    frame = pd.DataFrame(
        {
            "day": ["D1"] * 100 + ["D5"] * 100,
            "well": ["A01"] * 200,
            "f1": list(np.arange(100.0)) + list(np.arange(100.0) * 1000),
        }
    )
    trimmed, _ = trim_extremes(frame, ["f1"], SCOPE, 5.0, 95.0)
    kept = trimmed.groupby("day", sort=True)["f1"]
    assert kept.min().to_dict() == pytest.approx({"D1": 5.0, "D5": 5000.0})


def test_all_nan_feature_survives_without_error():
    frame = sample()
    frame["f2"] = np.nan
    trimmed, report = trim_extremes(frame, FEATURES, SCOPE, 2.5, 97.5)
    assert trimmed["f2"].isna().all()
    assert "f2" in report.features_with_empty_cells
    assert report.empty_cells == report.groups


def test_report_counts_add_up():
    frame = sample()
    trimmed, report = trim_extremes(frame, FEATURES, SCOPE, 2.5, 97.5)
    assert report.rows == len(frame)
    assert report.features == len(FEATURES)
    assert report.groups == frame.groupby(SCOPE).ngroups
    assert report.values_present == int(frame[FEATURES].notna().to_numpy().sum())
    actually_removed = int(
        (frame[FEATURES].notna() & trimmed[FEATURES].isna()).to_numpy().sum()
    )
    assert report.values_trimmed == actually_removed
    assert report.per_feature.sum() == report.values_trimmed
    assert list(report.per_feature.index) == FEATURES


def test_block_size_does_not_change_the_result():
    frame = sample()
    one, _ = trim_extremes(frame, FEATURES, SCOPE, 2.5, 97.5, block_size=1)
    many, _ = trim_extremes(frame, FEATURES, SCOPE, 2.5, 97.5, block_size=99)
    assert np.array_equal(
        one[FEATURES].to_numpy(), many[FEATURES].to_numpy(), equal_nan=True
    )


def test_inplace_leaves_the_caller_a_choice():
    frame = sample()
    copy, _ = trim_extremes(frame, FEATURES, SCOPE, 2.5, 97.5)
    assert copy is not frame
    assert frame[FEATURES].notna().to_numpy().sum() > copy[FEATURES].notna().to_numpy().sum()

    same = frame.copy()
    result, _ = trim_extremes(same, FEATURES, SCOPE, 2.5, 97.5, inplace=True)
    assert result is same
    assert np.array_equal(
        result[FEATURES].to_numpy(), copy[FEATURES].to_numpy(), equal_nan=True
    )


def test_rows_with_a_missing_scope_value_are_still_trimmed():
    # With pandas' default dropna they would be dropped from the grouped quantiles and
    # so escape trimming entirely, which is silent and wrong.
    frame = sample()
    frame.loc[frame.index[:50], "well"] = None
    _, report = trim_extremes(frame, FEATURES, SCOPE, 2.5, 97.5)
    assert report.rows_with_missing_scope == 50
    assert report.values_trimmed > 0


def test_unknown_column_is_an_error():
    with pytest.raises(KeyError, match="nope"):
        trim_extremes(sample(), ["nope"], SCOPE, 2.5, 97.5)


@pytest.mark.parametrize(("lower", "upper"), [(97.5, 2.5), (-1.0, 50.0), (0.0, 101.0)])
def test_impossible_percentiles_are_an_error(lower, upper):
    with pytest.raises(ValueError):
        trim_extremes(sample(), FEATURES, SCOPE, lower, upper)
