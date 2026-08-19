"""Applying a feature list back to the table.

The stage is deliberately dumb - select columns, count what is present - so the tests are
about the two things that are easy to get wrong: which column order comes out, and whether a
feature named in the list but absent from the table is a silent no-op or a reported one.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from cdr_fs.subset import subset_table

METADATA = ["Metadata_Day", "Metadata_Well"]


def table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Metadata_Day": ["D1", "D1", "D5", "D5"],
            "Metadata_Well": ["A01", "A02", "A01", "A02"],
            "zebra": [1.0, 2.0, 3.0, 4.0],
            "apple": [1.0, np.nan, np.nan, np.nan],
            "flat": [7.0, 7.0, 7.0, 7.0],
            "unused": [0.0, 0.0, 0.0, 0.0],
        }
    )


# -------------------------------------------------------------- the missing-data filter


@pytest.mark.parametrize(
    "max_missing, dropped",
    [
        (24, True),  # missing 25% >= 24
        (25, True),  # the boundary is inclusive: at the threshold it goes
        (26, False),
        (100, False),  # only an entirely empty feature reaches 100%
    ],
)
def test_the_threshold_is_inclusive(max_missing, dropped):
    """A feature missing exactly the threshold is dropped, as the published rule had it."""
    frame = pd.DataFrame(
        {"Metadata_Day": list("abcd"), "quarter": [1.0, 2.0, 3.0, np.nan]}
    )
    _, quality, _ = subset_table(
        frame, ["Metadata_Day"], ["quarter"], drop_missing=True, max_missing=max_missing
    )
    assert bool(quality["dropped"].iloc[0]) is dropped


# ------------------------------------------------------------------ named exclusions
