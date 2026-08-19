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

from cdr_fs.config import load_config
from cdr_fs.prune import (
    aggregate_units,
    cluster_features,
)

#: Eight values with no pattern, used as the base series every case is derived from.
BASE = np.array([1.0, 4.0, 2.0, 8.0, 5.0, 7.0, 3.0, 6.0])


def config_for(tmp_path, **overrides):
    return load_config(write_config(tmp_path, overrides))


def frame_with(values: dict[str, list[float]], *, unit: list[str]) -> pd.DataFrame:
    """An object-level table: one metadata column naming the unit, plus features."""
    return pd.DataFrame({"Metadata_Well": unit, **values})


# ------------------------------------------------------------------- what clusters what


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


# ------------------------------------------------------------------------- aggregation


def test_aggregation_takes_the_median_per_unit():
    frame = frame_with({"f": [1.0, 2.0, 30.0, 4.0, 5.0, 6.0]}, unit=list("AAABBB"))
    aggregated = aggregate_units(frame, ["f"], ["Metadata_Well"], fill_missing="none")
    assert aggregated["f"].tolist() == [2.0, 5.0]


# ------------------------------------------------------------------------ the linkage


# ------------------------------------------------------------------------- the stage


# ---------------------------------------------------------------------- configuration
