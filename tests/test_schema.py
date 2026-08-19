"""The metadata/feature split, checked against the real published column list.

`fixtures/columns_published.txt` holds the 481 column names of the full 3.7 GB table, so
these run in CI without the data. The arithmetic they lock is

    481 columns - 10 metadata = 471 features
                             = 469 rp_norm_* + counts_RelateLysoCell + counts_RelateMitoCell

and, explicitly, that a `^counts_` metadata pattern loses two of them.
"""

from __future__ import annotations

from cases import COLUMNS_PUBLISHED, METADATA_PATTERNS

from cdr_fs.schema import resolve_schema

PUBLISHED_COLUMNS = COLUMNS_PUBLISHED.read_text(encoding="utf-8").split()

ORGANELLE_COUNTS = ("counts_RelateLysoCell", "counts_RelateMitoCell")
QC_COUNTS = ("counts_Cells", "counts_Cytoplasm", "counts_FilteredNuclei")


def test_published_split_is_10_metadata_and_471_features():
    resolved = resolve_schema(PUBLISHED_COLUMNS, METADATA_PATTERNS)
    assert len(resolved.metadata) == 10
    assert len(resolved.features) == 471
    assert resolved.prefix_breakdown() == [("rp_", 469), ("counts_", 2)]


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
