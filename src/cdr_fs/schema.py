"""Metadata-versus-feature resolution, and reading the table those columns describe.

The split is declared by *inversion*: `[schema] metadata_patterns` names the metadata,
and everything else is a feature. That way adding a measured channel upstream does not
require touching the configuration, which is the common case.

The inversion has one sharp edge, and it is worth stating plainly because it decides
whether a run reproduces the published result. In the RTgill-W1 dataset five columns
begin with `counts_`:

    counts_Cells, counts_Cytoplasm, counts_FilteredNuclei     segmentation QC -> metadata
    counts_RelateLysoCell, counts_RelateMitoCell              organelle counts -> FEATURES

The last two are deliberate features: the article appends total lysosomal and
mitochondrial counts per cell to each cell's profile, and `counts_RelateLysoCell`
survives the whole selection into the published retained list. A convenient-looking
`^counts_` pattern would classify all five as metadata and silently drop a retained
feature. Enumerate the QC columns individually instead.

`ResolvedSchema.prefix_breakdown` exists to make that class of mistake visible: it
groups the resolved features by leading name token, so `counts_ -> 2` shows up next to
`rp_ -> 469` in the output of `cdr-fs check`.

`read_stage_table` is the other reader here: it reads back the tables the stages write,
rather than the input table, and exists because one of their columns does not survive a
naive round trip. See its docstring.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Iterable, Sequence

from cdr_fs.config import ConfigError

if TYPE_CHECKING:  # pragma: no cover - import cost is why this is guarded
    import pandas as pd

__all__ = [
    "ResolvedSchema",
    "read_header",
    "read_stage_table",
    "read_table",
    "resolve_schema",
]


@dataclass(frozen=True)
class ResolvedSchema:
    """The outcome of applying `metadata_patterns` to a table's columns.

    `metadata` and `features` together are exactly `columns`, each in the table's own
    column order - never sorted, because feature order decides tie-breaking downstream
    (the alphabetical cluster representative in `prune`, for one).
    """

    columns: tuple[str, ...]
    metadata: tuple[str, ...]
    features: tuple[str, ...]
    #: Patterns that matched no column at all - almost always a typo.
    unused_patterns: tuple[str, ...]

    @property
    def metadata_set(self) -> frozenset[str]:
        return frozenset(self.metadata)

    @property
    def feature_set(self) -> frozenset[str]:
        return frozenset(self.features)

    def prefix_breakdown(self, separator: str = "_") -> list[tuple[str, int]]:
        """Feature counts grouped by leading name token, largest group first.

        A quick sanity check on the split: the published table resolves to
        `[("rp_", 469), ("counts_", 2)]`, which is the intended 471 and shows at a glance
        that the two organelle counts stayed on the feature side rather than being
        swallowed by a `^counts_` metadata pattern.
        """
        counts: dict[str, int] = {}
        for feature in self.features:
            head = _prefix_of(feature, separator)
            counts[head] = counts.get(head, 0) + 1
        return sorted(counts.items(), key=lambda item: (-item[1], item[0]))

    def features_with_prefix(self, prefix: str, separator: str = "_") -> tuple[str, ...]:
        """The features `prefix_breakdown` counted under `prefix`, in table order."""
        return tuple(
            feature
            for feature in self.features
            if _prefix_of(feature, separator) == prefix
        )


def _prefix_of(name: str, separator: str) -> str:
    """The name up to and including its first separator, or the whole name if there is none.

    Keeping the separator makes the truncation visible in output: `rp_` reads as a family
    of names, where a bare `rp` looks like a column called `rp`.
    """
    head, found, _ = name.partition(separator)
    return head + found


def resolve_schema(
    columns: Sequence[str], patterns: Iterable[str | re.Pattern[str]]
) -> ResolvedSchema:
    """Split `columns` into metadata and features using `patterns`.

    A column is metadata when at least one pattern *searches* successfully in it, so an
    unanchored pattern matches anywhere in the name. Anchor with `^...$` to name one
    column exactly, which is what the published configuration does for every column
    that is not a `Metadata_` prefix.
    """
    compiled = []
    for pattern in patterns:
        if isinstance(pattern, re.Pattern):
            compiled.append(pattern)
        else:
            try:
                compiled.append(re.compile(pattern))
            except re.error as error:
                raise ConfigError(f"{pattern!r} is not a valid regex - {error}") from None

    metadata: list[str] = []
    features: list[str] = []
    used = set()
    for column in columns:
        matched = [index for index, rx in enumerate(compiled) if rx.search(column)]
        if matched:
            used.update(matched)
            metadata.append(column)
        else:
            features.append(column)

    return ResolvedSchema(
        columns=tuple(columns),
        metadata=tuple(metadata),
        features=tuple(features),
        unused_patterns=tuple(
            rx.pattern for index, rx in enumerate(compiled) if index not in used
        ),
    )


def read_header(path: str | Path, sep: str) -> list[str]:
    """Read just the column names. Constant time regardless of how large the table is."""
    path = Path(path)
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.reader(handle, delimiter=sep):
                if not row:
                    continue
                return [name.strip() for name in row]
    except UnicodeDecodeError as error:
        raise ConfigError(f"{path}: header is not valid UTF-8 text ({error})") from None
    except OSError as error:
        raise ConfigError(f"{path}: could not be read - {error}") from None
    raise ConfigError(f"{path}: is empty, so it has no header row")


def read_table(
    path: str | Path,
    sep: str,
    *,
    metadata: Sequence[str],
    features: Sequence[str] | None = None,
    strip: bool = True,
) -> pd.DataFrame:
    """Read the per-object table with metadata as text and features as float64.

    Metadata columns are read as text on purpose. Level labels are compared against the
    strings in `[design]`, so a column of `10, 11` must not become integers on the way
    in - that mismatch between a parsed int and a configured string is a whole class of
    bug this avoids.

    `strip` removes surrounding whitespace from metadata values, reproducing the
    normalisation the published scripts applied. On the published table it changes
    nothing (verified: zero affected rows), but it costs little and a hand-edited table
    is a real possibility.
    """
    import numpy as np
    import pandas as pd

    path = Path(path)
    selected = list(metadata) + (list(features) if features is not None else [])
    dtypes: dict[str, object] = {name: str for name in metadata}
    if features is not None:
        dtypes.update({name: np.float64 for name in features})

    try:
        frame = pd.read_csv(
            path,
            sep=sep,
            usecols=selected if features is not None else None,
            dtype=dtypes,
        )
    except ValueError as error:
        raise ConfigError(f"{path}: could not be read as a table - {error}") from None

    if features is not None:
        # usecols ignores the order it is given; restore the declared order so that
        # downstream feature indexing is predictable.
        frame = frame[selected]
    if strip:
        for column in metadata:
            frame[column] = frame[column].str.strip()
    return frame


def read_stage_table(path: str | Path, sep: str) -> pd.DataFrame:
    """Read a table written by one of the stages - `emd.tsv`, `fit.tsv`, and so on.

    Use this rather than a bare `pandas.read_csv`, because the empty string is a
    meaningful value in these tables and does not survive the round trip on its own. An
    experiment with no `[schema] group_by` has exactly one stratum, whose label is `""`;
    written to a delimited file that is an empty field, and pandas reads an empty field as
    NaN whatever dtype it is given. Every stage downstream groups by `stratum`, and
    `groupby` drops NaN keys by default - so the whole table would vanish, quietly, and
    only for the unstratified design. Filling those back to `""` on the way in is what
    keeps a single-stratum run working.
    """
    import pandas as pd

    frame = pd.read_csv(path, sep=sep, dtype={"feature": str, "stratum": str})
    for column in ("feature", "stratum"):
        if column in frame.columns:
            frame[column] = frame[column].fillna("")
    return frame
