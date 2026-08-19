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


def sample_with_infinities(rows: int = 600, seed: int = 7) -> pd.DataFrame:
    """Like `sample`, but with the infinities real AreaShape_FormFactor columns carry.

    FormFactor is 4*pi*Area / Perimeter^2, so an organelle object whose perimeter rounds to
    zero gives inf. `f2` gets a scattering of them, `f3` one single infinity, which is
    enough to empty every group it touches.
    """
    frame = sample(rows, seed)
    rng = np.random.default_rng(seed + 1)
    frame.loc[rng.choice(rows, 25, replace=False), "f2"] = np.inf
    frame.loc[frame.index[0], "f3"] = np.inf
    frame.loc[frame.index[1], "f3"] = -np.inf
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


@pytest.mark.parametrize("builder", [sample, sample_with_infinities], ids=["finite", "with-inf"])
@pytest.mark.parametrize(("lower", "upper"), [(2.5, 97.5), (0.0, 100.0), (10.0, 90.0)])
def test_matches_the_original_per_group_loop(builder, lower, upper):
    frame = builder()
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


def test_finite_group_is_never_emptied():
    # With finite data the values sitting exactly on the percentile bounds survive, so a
    # group that had data still has data.
    frame = pd.DataFrame({"day": ["D1"] * 5, "well": ["A01"] * 5, "f1": [1.0, 2, 3, 4, 5]})
    trimmed, report = trim_extremes(frame, ["f1"], SCOPE, 2.5, 97.5)
    assert trimmed["f1"].notna().sum() >= 1
    assert report.empty_cells == 0


def test_one_infinity_empties_its_whole_group():
    """The mechanism behind two features missing from the published EMD tables.

    A percentile interpolating into an infinite tail is NaN, nothing satisfies
    `value <= nan`, and the entire (group, feature) cell is discarded. One infinity in a
    group of a thousand is enough. Reproduced deliberately: it is what the published run
    did, and `drop = (v < lo) | (v > hi)` would *not* reproduce it.
    """
    frame = pd.DataFrame(
        {
            "day": ["D1"] * 10 + ["D5"] * 10,
            "well": ["A01"] * 20,
            "f1": [*np.arange(9.0), np.inf, *np.arange(10.0)],
        }
    )
    trimmed, report = trim_extremes(frame, ["f1"], SCOPE, 2.5, 97.5)
    assert trimmed.loc[trimmed["day"] == "D1", "f1"].isna().all()  # emptied
    assert trimmed.loc[trimmed["day"] == "D5", "f1"].notna().any()  # untouched neighbour
    assert report.values_nonfinite == 1
    assert report.features_with_nonfinite == ("f1",)
    assert report.empty_cells == 1


def test_whether_an_infinity_empties_a_cell_depends_on_how_many_there_are():
    """What an infinity does depends on where the percentile's bracketing pair lands.

    numpy interpolates a percentile two different ways. With the fractional position
    `t >= 0.5` it computes `b - (b-a)*(1-t)`, so an infinite `b` gives `inf - inf = NaN`,
    the keep mask is empty, and the whole cell is discarded. With `t < 0.5` it computes
    `a + (b-a)*t`, which gives `+inf`, and then `inf <= inf` holds and the infinity is
    *kept* - it survives trimming and flows on into the distance computation.

    Which branch applies depends on the group size, so "this feature carries infinities"
    does not tell you what happened to it. Hence the report states both the infinities
    found and the cells actually emptied. All of this is inherited behaviour: the
    parametrized oracle test above confirms it matches the original loop exactly.
    """
    frame = sample_with_infinities()
    trimmed, report = trim_extremes(frame, FEATURES, SCOPE, 2.5, 97.5)
    assert report.values_nonfinite == 27
    assert set(report.features_with_nonfinite) == {"f2", "f3"}
    # f2 carries ~2 infinities per group, enough to reach the p97.5 bracket in most.
    assert set(report.features_with_empty_cells) == {"f2"}
    # f3's lone infinities fall outside the bracket, so they are merely dropped.
    assert trimmed["f3"].notna().sum() > 0
    assert not np.isinf(trimmed["f3"].to_numpy()).any()
    # ...but two of f2's land in groups small enough for the t < 0.5 branch, and survive.
    assert np.isinf(trimmed["f2"].to_numpy()).sum() == 2
    assert "infinite value(s)" in report.summary()


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
