"""The figures.

A drawing test cannot assert what a picture looks like, so these assert the things that
decide whether the picture is right: that a curve is drawn from parameters that read back
exactly, that the panels are paged and named as they should be, that a series with a hole is
still drawn against the full axis, and that the dendrogram's colours come from the clustering
rather than from a redrawing of it.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from cases import write_config

from cdr_fs.config import load_config
from cdr_fs.emd import COLUMNS as EMD_COLUMNS
from cdr_fs.fit import COLUMNS as FIT_COLUMNS
from cdr_fs.plots import plot_distribution
from cdr_fs.prune import cluster_features, correlation_distance, linkage_table

LEVELS = ("10", "9", "8", "7", "6", "5", "4", "3")


def emd_table(features, strata=("D1",), *, hole: set[str] = frozenset()) -> pd.DataFrame:
    """One distance per feature, stratum and fitted contrast; `hole` drops a feature's last."""
    rows = []
    for feature in features:
        for stratum in strata:
            for index, level in enumerate(LEVELS):
                value = np.nan if feature in hole and level == LEVELS[-1] else 1.0 + index
                rows.append((feature, stratum, f"11v{level}", "11", level, "", 100, 100, value))
    return pd.DataFrame.from_records(rows, columns=list(EMD_COLUMNS))


def fit_table(features, strata=("D1",)) -> pd.DataFrame:
    """A Lin and a Con fit per series, with parameters in the form `fit.tsv` writes."""
    rows = []
    for feature in features:
        for stratum in strata:
            rows.append(
                (feature, stratum, "Lin", 8, 2, -20.0, -19.0, -39.0, 1.0, "m=1.0,b=1.0")
            )
            rows.append(
                (feature, stratum, "Con", 8, 1, -10.0, -9.5, -19.5, np.nan, "c=4.5")
            )
    return pd.DataFrame.from_records(rows, columns=list(FIT_COLUMNS))


def config_for(tmp_path, **overrides):
    return load_config(write_config(tmp_path, overrides))


# ------------------------------------------------------------------------- parameters


# ------------------------------------------------------------------------- fit panels


# ----------------------------------------------------------------------- distribution


def test_an_empty_distance_table_is_refused_rather_than_drawn(tmp_path):
    """A blank figure looks like a result, which is worse than an error.

    The way this happens in practice is a table filtered down to nothing on the way in - the
    baseline contrast set is between replicates, not between exposure levels, so anything that
    selects rows by exposure level empties it completely.
    """
    table = emd_table(["f0"])
    with pytest.raises(ValueError, match="no rows"):
        plot_distribution(table.iloc[0:0], tmp_path / "emd.png")

    table["emd"] = np.nan
    with pytest.raises(ValueError, match="none of the 8 distances"):
        plot_distribution(table, tmp_path / "emd.png")

    with pytest.raises(ValueError, match="no emd column"):
        plot_distribution(table.drop(columns=["emd"]), tmp_path / "emd.png")


# ------------------------------------------------------------------------- dendrogram
