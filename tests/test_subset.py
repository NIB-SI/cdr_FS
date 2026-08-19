"""Applying a feature list back to the table.

The stage is deliberately dumb - select columns, count what is present - so the tests are
about the two things that are easy to get wrong: which column order comes out, and whether a
feature named in the list but absent from the table is a silent no-op or a reported one.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from cdr_fs.subset import QUALITY_COLUMNS, read_feature_list, subset_table

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


def test_columns_come_out_in_the_table_order_not_the_list_order():
    """Two lists over one table then give column-comparable files."""
    subset, _, report = subset_table(table(), METADATA, ["apple", "zebra"])
    assert list(subset.columns) == ["Metadata_Day", "Metadata_Well", "zebra", "apple"]
    assert report.matched == 2
    assert report.rows == 4


def test_a_feature_that_is_not_a_column_is_reported_not_ignored():
    subset, quality, report = subset_table(table(), METADATA, ["zebra", "ghost"])
    assert "ghost" not in subset.columns
    assert report.absent == ("ghost",)
    assert report.requested == 2 and report.matched == 1
    assert "not columns of the table" in report.summary()
    assert quality["feature"].tolist() == ["zebra"]


def test_the_quality_table_is_what_a_missing_data_rule_needs():
    """Per feature: how much data is there, and is it constant. No threshold is applied."""
    _, quality, report = subset_table(table(), METADATA, ["zebra", "apple", "flat"])
    assert list(quality.columns) == list(QUALITY_COLUMNS)
    quality = quality.set_index("feature")
    assert quality.loc["zebra", "n_present"] == 4
    assert quality.loc["apple", "n_present"] == 1
    assert quality.loc["apple", "nonmissing_fraction"] == pytest.approx(0.25)
    assert quality.loc["flat", "n_distinct"] == 1
    assert bool(quality.loc["flat", "constant"]) is True
    assert bool(quality.loc["zebra", "constant"]) is False
    # apple has a single surviving value, which is thin rather than constant - a variance
    # over one observation is undefined, and the published filter read it the same way.
    assert bool(quality.loc["apple", "constant"]) is False
    assert report.constant == ("flat",)
    assert report.values_present == 4 + 1 + 4
    assert report.values_total == 12


def test_a_feature_with_no_data_at_all_is_named():
    frame = table()
    frame["gone"] = np.nan
    _, quality, report = subset_table(frame, METADATA, ["zebra", "gone"])
    assert report.all_missing == ("gone",)
    assert report.constant == ()  # nothing to be constant about
    assert "no data at all" in report.summary()


def test_an_empty_list_yields_metadata_only():
    subset, quality, report = subset_table(table(), METADATA, [])
    assert list(subset.columns) == METADATA
    assert quality.empty
    assert report.matched == 0
    assert report.fraction_present == 0.0


def test_the_list_reader_drops_blanks_and_duplicates(tmp_path):
    path = tmp_path / "features.txt"
    path.write_text("alpha\n\n  beta  \nalpha\n\n", encoding="utf-8")
    assert read_feature_list(path) == ["alpha", "beta"]
