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
