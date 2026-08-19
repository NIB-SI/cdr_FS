"""Apply a feature list back to the object-level table.

The last step of a run: take the features that survived - from `select`, or from `prune` when
it is enabled - and write the input table restricted to them, with the configured trim
applied. That file is what a downstream analysis consumes, and it is the only stage whose
output is data rather than evidence.

Two scripts in the original pipeline did this, identical but for which list they read and
where they wrote; one stage with an explicit list argument replaces both.

## It also reports what a downstream filter would remove

Trimming removes values rather than rows, so the subset is not a rectangle of data: a feature
can be missing for most objects and still be present as a column. Anything consuming the
subset has to decide what to do about that, and the published pipeline did - it dropped
features with 30% or more missing values, at the point of building its dimension-reduction
subsets, not here.

That threshold is a property of the analysis downstream, not of the selection, so this stage
does not apply one. It writes the numbers the decision needs instead: per feature, how many
values are present, what fraction that is, and whether the feature is constant or empty. A
feature that is constant or entirely missing is worth acting on whatever the downstream
analysis is, so the report names those outright.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Sequence

if TYPE_CHECKING:  # pragma: no cover
    import pandas as pd

__all__ = ["QUALITY_COLUMNS", "SubsetReport", "read_feature_list", "subset_table"]

#: Column order of the per-feature quality table.
QUALITY_COLUMNS = ("feature", "n_present", "nonmissing_fraction", "n_distinct", "constant")


@dataclass(frozen=True)
class SubsetReport:
    """What the subset contains, and what is thin inside it."""

    rows: int
    requested: int
    matched: int
    #: Named in the feature list but absent from the table.
    absent: tuple[str, ...]
    values_present: int
    values_total: int
    #: Features that have data and never vary. A feature with one surviving value is thin,
    #: not constant, and is reported by its fraction instead - which is also how the
    #: published downstream filter read it, since a variance over one value is undefined.
    constant: tuple[str, ...]
    all_missing: tuple[str, ...]
    #: (feature, non-missing fraction) for the thinnest features, thinnest first.
    thinnest: tuple[tuple[str, float], ...]

    @property
    def fraction_present(self) -> float:
        return self.values_present / self.values_total if self.values_total else 0.0

    def summary(self) -> str:
        lines = [
            f"subset {self.rows:,} row(s) x {self.matched} feature(s)",
            f"  {self.values_present:,} of {self.values_total:,} feature values present "
            f"({self.fraction_present:.2%})",
        ]
        if self.absent:
            lines.append(
                f"  {len(self.absent)} feature(s) in the list are not columns of the table "
                f"and were skipped: {_listing(self.absent)}"
            )
        if self.thinnest:
            shown = ", ".join(
                f"{name} {fraction:.1%}" for name, fraction in self.thinnest
            )
            lines.append(f"  thinnest feature(s) by data present: {shown}")
        if self.all_missing:
            lines.append(
                f"  {len(self.all_missing)} feature(s) have no data at all: "
                f"{_listing(self.all_missing)}"
            )
        if self.constant:
            lines.append(
                f"  {len(self.constant)} feature(s) take a single value: "
                f"{_listing(self.constant)}"
            )
        return "\n".join(lines)


def _listing(names: Sequence[str], limit: int = 3) -> str:
    shown = ", ".join(names[:limit])
    return shown if len(names) <= limit else f"{shown}, +{len(names) - limit} more"


def read_feature_list(path: str | Path) -> list[str]:
    """Read a one-name-per-line feature list, dropping blanks and duplicates."""
    text = Path(path).read_text(encoding="utf-8")
    names = [line.strip() for line in text.splitlines() if line.strip()]
    return list(dict.fromkeys(names))


def subset_table(
    frame: pd.DataFrame,
    metadata: Sequence[str],
    features: Sequence[str],
) -> tuple[pd.DataFrame, pd.DataFrame, SubsetReport]:
    """Restrict `frame` to `metadata` plus whichever of `features` it actually has.

    Columns come out in the table's own order rather than the list's, so that two lists over
    the same table give column-comparable files. Returns the subset, a per-feature quality
    table, and a report.
    """
    import numpy as np
    import pandas as pd

    requested = list(dict.fromkeys(features))
    wanted = set(requested)
    matched = [column for column in frame.columns if column in wanted]
    absent = [name for name in requested if name not in frame.columns]
    kept = [column for column in frame.columns if column in set(metadata)] + matched
    subset = frame[kept]

    values = frame[matched].to_numpy(dtype=np.float64) if matched else np.empty((len(frame), 0))
    present = ~np.isnan(values)
    counts = present.sum(axis=0)
    distinct = np.array(
        [np.unique(values[present[:, index], index]).size for index in range(len(matched))],
        dtype=np.int64,
    )
    fraction = counts / len(frame) if len(frame) else np.zeros_like(counts, dtype=float)

    quality = pd.DataFrame(
        {
            "feature": matched,
            "n_present": counts,
            "nonmissing_fraction": fraction,
            "n_distinct": distinct,
            "constant": (distinct == 1) & (counts >= 2),
        },
        columns=list(QUALITY_COLUMNS),
    )
    order = np.argsort(fraction, kind="stable")
    report = SubsetReport(
        rows=len(frame),
        requested=len(requested),
        matched=len(matched),
        absent=tuple(absent),
        values_present=int(counts.sum()),
        values_total=len(frame) * len(matched),
        constant=tuple(quality.loc[quality["constant"], "feature"]),
        all_missing=tuple(quality.loc[quality["n_present"] == 0, "feature"]),
        thinnest=tuple(
            (matched[index], float(fraction[index]))
            for index in order[:3]
            if fraction[index] < 1.0
        ),
    )
    return subset, quality, report
