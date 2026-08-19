"""Apply a feature list back to the object-level table.

The last step of a run: take the features that survived - from `select`, or from `prune` when
it is enabled - and write the input table restricted to them, with the configured trim
applied. That file is what a downstream analysis consumes, and it is the only stage whose
output is data rather than evidence.

Two scripts in the original pipeline did this, identical but for which list they read and
where they wrote; one stage with an explicit list argument replaces both.

## The missing-data filter

Trimming removes values rather than rows, so the subset is not a rectangle of data: a feature
can be missing for most objects and still be present as a column. `[subset] drop_missing`
removes the columns that are too empty to be worth carrying - a feature missing
`[subset] max_missing` percent of the table **or more** is dropped, 30% by default.

Where it is applied matters as much as the threshold. It runs here, over the **whole** table,
before anything downstream subsamples: one decision for the whole experiment, so a feature is
either in the analysis or out of it. Deciding per subsample instead lets the same feature be
present in one day's file and absent from another, which makes the resulting sets
incomparable - and it is what the original pipeline did, in its dimension-reduction step.

Everything else is left alone and reported rather than filtered. A constant feature or one
with a single surviving value is named in the report and flagged in the per-feature table;
what to do about it depends on the analysis, and guessing is worse than saying so.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Sequence

if TYPE_CHECKING:  # pragma: no cover
    import pandas as pd

__all__ = ["QUALITY_COLUMNS", "SubsetReport", "read_feature_list", "subset_table"]

#: Column order of the per-feature quality table.
QUALITY_COLUMNS = (
    "feature",
    "n_present",
    "nonmissing_fraction",
    "n_distinct",
    "constant",
    "dropped",
)


@dataclass(frozen=True)
class SubsetReport:
    """What the subset contains, what was filtered out of it, and what is thin inside it."""

    rows: int
    requested: int
    #: Requested features that are columns of the table, before the missing-data filter.
    matched: int
    #: Named in the feature list but absent from the table.
    absent: tuple[str, ...]
    #: Removed by the missing-data filter, with the fraction of the table each was missing.
    dropped: tuple[tuple[str, float], ...]
    drop_missing: bool
    max_missing: float
    values_present: int
    values_total: int
    #: Features that have data and never vary. A feature with one surviving value is thin,
    #: not constant, and is reported by its fraction instead.
    constant: tuple[str, ...]
    all_missing: tuple[str, ...]
    #: (feature, non-missing fraction) for the thinnest features kept, thinnest first.
    thinnest: tuple[tuple[str, float], ...]

    @property
    def kept(self) -> int:
        return self.matched - len(self.dropped)

    @property
    def fraction_present(self) -> float:
        return self.values_present / self.values_total if self.values_total else 0.0

    def summary(self) -> str:
        lines = [
            f"subset {self.rows:,} row(s) x {self.kept} feature(s)",
            f"  {self.values_present:,} of {self.values_total:,} feature values present "
            f"({self.fraction_present:.2%})",
        ]
        if self.absent:
            lines.append(
                f"  {len(self.absent)} feature(s) in the list are not columns of the table "
                f"and were skipped: {_listing(self.absent)}"
            )
        if self.drop_missing:
            lines.append(
                f"  dropped {len(self.dropped)} of {self.matched} feature(s) missing "
                f"{self.max_missing:g}% of the table or more"
                + (":" if self.dropped else "")
            )
            lines.extend(
                f"    {name}  {fraction:.1%} missing" for name, fraction in self.dropped
            )
        if self.thinnest:
            shown = ", ".join(f"{name} {fraction:.1%}" for name, fraction in self.thinnest)
            lines.append(f"  thinnest feature(s) kept, by data present: {shown}")
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
    *,
    drop_missing: bool = False,
    max_missing: float = 30.0,
) -> tuple[pd.DataFrame, pd.DataFrame, SubsetReport]:
    """Restrict `frame` to `metadata` plus whichever of `features` it actually has.

    With `drop_missing`, a feature missing `max_missing` percent of the rows or more is left
    out. Columns come out in the table's own order rather than the list's, so that two lists
    over the same table give column-comparable files. Returns the subset, a per-feature
    quality table covering every matched feature including the dropped ones, and a report.
    """
    import numpy as np
    import pandas as pd

    requested = list(dict.fromkeys(features))
    wanted = set(requested)
    matched = [column for column in frame.columns if column in wanted]
    absent = [name for name in requested if name not in frame.columns]

    values = (
        frame[matched].to_numpy(dtype=np.float64)
        if matched
        else np.empty((len(frame), 0), dtype=np.float64)
    )
    present = ~np.isnan(values)
    counts = present.sum(axis=0)
    distinct = np.array(
        [np.unique(values[present[:, index], index]).size for index in range(len(matched))],
        dtype=np.int64,
    )
    fraction = counts / len(frame) if len(frame) else np.zeros(len(matched))
    # Compared as counts rather than as percentages, so that a feature missing exactly the
    # threshold falls on the documented side of it whatever the rounding.
    dropped = (
        (len(frame) - counts) * 100 >= max_missing * len(frame)
        if drop_missing
        else np.zeros(len(matched), dtype=bool)
    )

    quality = pd.DataFrame(
        {
            "feature": matched,
            "n_present": counts,
            "nonmissing_fraction": fraction,
            "n_distinct": distinct,
            "constant": (distinct == 1) & (counts >= 2),
            "dropped": dropped,
        },
        columns=list(QUALITY_COLUMNS),
    )
    kept = [name for name, gone in zip(matched, dropped) if not gone]
    subset = frame[[column for column in frame.columns if column in set(metadata)] + kept]

    surviving = quality[~quality["dropped"]]
    order = surviving["nonmissing_fraction"].to_numpy().argsort(kind="stable")
    report = SubsetReport(
        rows=len(frame),
        requested=len(requested),
        matched=len(matched),
        absent=tuple(absent),
        dropped=tuple(
            (name, 1.0 - value)
            for name, value in zip(
                quality.loc[quality["dropped"], "feature"],
                quality.loc[quality["dropped"], "nonmissing_fraction"],
            )
        ),
        drop_missing=drop_missing,
        max_missing=max_missing,
        values_present=int(surviving["n_present"].sum()),
        values_total=len(frame) * len(kept),
        constant=tuple(surviving.loc[surviving["constant"], "feature"]),
        all_missing=tuple(surviving.loc[surviving["n_present"] == 0, "feature"]),
        thinnest=tuple(
            (surviving["feature"].iloc[index], float(surviving["nonmissing_fraction"].iloc[index]))
            for index in order[:3]
            if surviving["nonmissing_fraction"].iloc[index] < 1.0
        ),
    )
    return subset, quality, report
