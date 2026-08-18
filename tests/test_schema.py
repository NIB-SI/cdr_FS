"""The metadata/feature split, checked against the real published column list.

`fixtures/columns_published.txt` holds the 481 column names of the full 3.7 GB table, so
these run in CI without the data. The arithmetic they lock is

    481 columns - 10 metadata = 471 features
                             = 469 rp_norm_* + counts_RelateLysoCell + counts_RelateMitoCell

and, explicitly, that a `^counts_` metadata pattern loses two of them.
"""

from __future__ import annotations

from cases import COLUMNS_PUBLISHED, METADATA_PATTERNS, SUBSET

from cdr_fs.schema import read_header, resolve_schema

PUBLISHED_COLUMNS = COLUMNS_PUBLISHED.read_text(encoding="utf-8").split()

ORGANELLE_COUNTS = ("counts_RelateLysoCell", "counts_RelateMitoCell")
QC_COUNTS = ("counts_Cells", "counts_Cytoplasm", "counts_FilteredNuclei")


def test_published_table_has_481_columns():
    assert len(PUBLISHED_COLUMNS) == 481


def test_published_split_is_10_metadata_and_471_features():
    resolved = resolve_schema(PUBLISHED_COLUMNS, METADATA_PATTERNS)
    assert len(resolved.metadata) == 10
    assert len(resolved.features) == 471
    assert resolved.prefix_breakdown() == [("rp_", 469), ("counts_", 2)]


def test_both_organelle_counts_are_features():
    # The article appends total lysosomal and mitochondrial counts per cell to each
    # cell's profile; counts_RelateLysoCell survives into the published retained list.
    resolved = resolve_schema(PUBLISHED_COLUMNS, METADATA_PATTERNS)
    for column in ORGANELLE_COUNTS:
        assert column in resolved.feature_set
    for column in QC_COUNTS:
        assert column in resolved.metadata_set


def test_bare_counts_pattern_silently_loses_two_retained_features():
    """The counts_ trap, asserted so nobody re-introduces the shortcut.

    Collapsing the three QC columns into `^counts_` is the single most likely way to
    produce a run that looks right and is not: it costs two features, one of which the
    published selection kept.
    """
    lazy = [p for p in METADATA_PATTERNS if not p.startswith("^counts_")] + ["^counts_"]
    resolved = resolve_schema(PUBLISHED_COLUMNS, lazy)
    assert len(resolved.features) == 469
    assert set(ORGANELLE_COUNTS) <= resolved.metadata_set


def test_features_keep_the_table_column_order():
    # Feature order is not cosmetic: the alphabetical cluster representative in `prune`
    # is chosen from it, so re-sorting here would change which features are retained.
    resolved = resolve_schema(PUBLISHED_COLUMNS, METADATA_PATTERNS)
    assert list(resolved.metadata) + list(resolved.features) != PUBLISHED_COLUMNS
    assert [c for c in PUBLISHED_COLUMNS if c in resolved.feature_set] == list(
        resolved.features
    )


def test_patterns_are_searched_not_fullmatched():
    resolved = resolve_schema(["Metadata_Day", "x_Metadata_Day"], ["Metadata_"])
    assert resolved.features == ()
    anchored = resolve_schema(["Metadata_Day", "x_Metadata_Day"], ["^Metadata_"])
    assert anchored.features == ("x_Metadata_Day",)


def test_unused_patterns_are_reported():
    resolved = resolve_schema(["a", "b"], ["^a$", "^zzz$"])
    assert resolved.unused_patterns == ("^zzz$",)


def test_subset_fixture_matches_the_published_schema():
    # The committed subset must stay a faithful slice: same metadata, same split rule.
    resolved = resolve_schema(read_header(SUBSET, "\t"), METADATA_PATTERNS)
    assert len(resolved.metadata) == 10
    assert resolved.prefix_breakdown() == [("rp_", 18), ("counts_", 2)]
