"""Laying the distances out along the exposure axis, and fitting them.

The distances are synthesised here rather than computed, so a series can be given an exact
shape - flat, rising, hormetic - and the fit checked against what that shape ought to
produce. `test_golden.py` covers the real numbers.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from cases import write_config

from cdr_fs.config import load_config
from cdr_fs.emd import COLUMNS as EMD_COLUMNS
from cdr_fs.fit import series_from_emd

STRATA = ("D1", "D5", "D7", "D9")
FITTED = ("10", "9", "8", "7", "6", "5", "4", "3")
#: The reference dose vector, so the `dose` axis can be exercised.
DOSE = ",".join(f"{1000 / 1.75**k:.6f}" for k in range(8, -1, -1))


def emd_table(shapes: dict[str, dict[str, np.ndarray]], include_top: bool = True):
    """Build a distance table. `shapes[feature][stratum]` is the 8-point series."""
    rows = []
    for feature, per_stratum in shapes.items():
        for stratum, values in per_stratum.items():
            for level, value in zip(FITTED, values):
                rows.append((feature, stratum, f"11v{level}", "11", level, "", 100, 100, value))
            if include_top:
                # The withheld top dose is measured but must never be fitted.
                rows.append((feature, stratum, "11v2", "11", "2", "", 100, 100, 999.0))
    return pd.DataFrame.from_records(rows, columns=list(EMD_COLUMNS))


def rising(scale: float = 1.0) -> np.ndarray:
    return np.arange(8, dtype=float) * scale


def flat(value: float = 2.0) -> np.ndarray:
    return np.full(8, value)


def config_for(tmp_path, **overrides):
    return load_config(write_config(tmp_path, overrides))


def test_the_withheld_top_dose_is_not_fitted(tmp_path):
    config = config_for(tmp_path)
    (_, _, x, y, _), = list(series_from_emd(config, emd_table({"f": {"D1": rising()}})))
    assert len(y) == 8
    assert 999.0 not in set(y)  # the 11v2 distance is present but excluded from the fit
