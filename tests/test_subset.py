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


def test_a_mistyped_exclusion_is_warned_about(tmp_path):
    """It fails in the dangerous direction - the feature stays in and nothing says so."""
    from cases import validate, write_config

    config = write_config(
        tmp_path, {"subset.exclude": "counts_RelateLysoCell,counts_RelateLysoCel"}
    )
    _, warnings = validate(config)
    assert any("counts_RelateLysoCel" in warning for warning in warnings)
    assert any("matched exactly" in warning for warning in warnings)


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


# -------------------------------------------------------------- the missing-data filter


def test_nothing_is_dropped_unless_the_filter_is_on():
    subset, quality, report = subset_table(table(), METADATA, ["zebra", "apple"])
    assert "apple" in subset.columns  # 75% missing, and kept
    assert not quality["dropped"].any()
    assert report.dropped == () and report.kept == 2


def test_a_feature_missing_too_much_of_the_table_is_dropped():
    subset, quality, report = subset_table(
        table(), METADATA, ["zebra", "apple"], drop_missing=True, max_missing=30
    )
    assert list(subset.columns) == [*METADATA, "zebra"]
    assert report.dropped == (("apple", 0.75),)
    assert report.matched == 2 and report.kept == 1
    # The dropped feature stays in the quality table: the record of why it went.
    assert quality["feature"].tolist() == ["zebra", "apple"]
    assert quality["dropped"].tolist() == [False, True]
    assert "missing 30% of the table or more" in report.summary()
    assert "apple  75.0% missing" in report.summary()


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


def test_an_entirely_empty_feature_goes_at_any_threshold():
    frame = table()
    frame["gone"] = np.nan
    subset, _, report = subset_table(
        frame, METADATA, ["zebra", "gone"], drop_missing=True, max_missing=100
    )
    assert "gone" not in subset.columns
    assert report.dropped == (("gone", 1.0),)
    # It was filtered out, so it is no longer among the problems left in the subset.
    assert report.all_missing == ()


def test_the_counts_describe_what_survived():
    _, _, report = subset_table(
        table(), METADATA, ["zebra", "apple", "flat"], drop_missing=True, max_missing=30
    )
    assert report.kept == 2
    assert report.values_total == 4 * 2  # zebra and flat, not apple
    assert report.values_present == 8
    assert report.fraction_present == 1.0


# ------------------------------------------------------------------ named exclusions


def test_a_named_feature_is_excluded_whatever_its_quality():
    subset, quality, report = subset_table(
        table(), METADATA, ["zebra", "flat"], exclude=["zebra"]
    )
    assert list(subset.columns) == [*METADATA, "flat"]
    assert report.excluded == ("zebra",)
    assert report.kept == 1
    assert quality.set_index("feature").loc["zebra", "drop_reason"] == "excluded"
    assert "excluded 1 feature(s) by name: zebra" in report.summary()


def test_an_exclusion_takes_precedence_over_the_missing_data_reason():
    """A feature struck out by hand is reported as excluded, not as thin.

    Both would remove it, but only one of them is a decision somebody made.
    """
    _, quality, report = subset_table(
        table(),
        METADATA,
        ["zebra", "apple"],
        drop_missing=True,
        max_missing=30,
        exclude=["apple"],
    )
    assert report.excluded == ("apple",)
    assert report.dropped == ()
    assert quality.set_index("feature").loc["apple", "drop_reason"] == "excluded"


def test_excluding_everything_leaves_metadata():
    subset, _, report = subset_table(
        table(), METADATA, ["zebra", "flat"], exclude=["zebra", "flat"]
    )
    assert list(subset.columns) == METADATA
    assert report.kept == 0


def test_an_exclusion_that_names_nothing_is_harmless_here():
    """`cdr-fs check` warns about it; this stage simply has nothing to remove."""
    _, quality, report = subset_table(
        table(), METADATA, ["zebra"], exclude=["ghost", "unused"]
    )
    assert report.excluded == ()
    assert not quality["dropped"].any()


def test_the_list_reader_drops_blanks_and_duplicates(tmp_path):
    path = tmp_path / "features.txt"
    path.write_text("alpha\n\n  beta  \nalpha\n\n", encoding="utf-8")
    assert read_feature_list(path) == ["alpha", "beta"]
