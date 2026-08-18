"""Optional extreme-value trimming.

Data-level quality control: within each group of `[trim] scope`, values below the lower
percentile or above the upper percentile are discarded. In the published run the scope is
one well of one replicate on one day, and the interval kept is `[p2.5, p97.5]` inclusive,
which removes badly segmented objects without reference to the exposure level.

Two properties of this step are surprising often enough to be worth stating:

1. **Trimming removes values, not rows.** An object that is an artifact in one feature is
   usually sound in the others, so each feature is trimmed independently.
2. **Trimmed values become NaN, so N differs per feature.** That is correct for
   per-feature QC, but it propagates: a feature with a missing EMD value in a stratum is
   never fitted, hence never selected. `TrimReport` therefore reports the (feature, group)
   cells that hold no data, because those are the ones that will quietly drop out later.

Because the kept interval is closed, trimming can never empty a group that had data -
the values sitting exactly on the percentile bounds always survive. Any empty cell the
report names was already empty on the way in.

The vectorized implementation here was checked against the per-group `np.nanpercentile`
loop of the original scripts and is identical value for value, NaNs included.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Sequence

if TYPE_CHECKING:  # pragma: no cover
    import pandas as pd

__all__ = ["TrimReport", "trim_extremes"]

#: Features processed per pass. The published table is 503,920 rows x 471 features, so a
#: whole-table temporary would be ~1.9 GB and three of them are live at once; blocking
#: bounds the peak at a few hundred MB without measurably changing the run time.
BLOCK_SIZE = 32


@dataclass(frozen=True)
class TrimReport:
    """What trimming did, in enough detail to put in a run manifest."""

    scope: tuple[str, ...]
    lower_percentile: float
    upper_percentile: float
    groups: int
    features: int
    rows: int
    #: Non-NaN feature values before trimming.
    values_present: int
    values_trimmed: int
    #: (feature, group) cells holding no data. These were already empty: see the module
    #: docstring. They are the cells that will drop out of the fit later.
    empty_cells: int
    features_with_empty_cells: tuple[str, ...]
    #: Rows whose scope columns are themselves missing; they are trimmed as one group.
    rows_with_missing_scope: int
    #: Values trimmed per feature, in feature order.
    per_feature: pd.Series

    @property
    def fraction_trimmed(self) -> float:
        return self.values_trimmed / self.values_present if self.values_present else 0.0

    def summary(self) -> str:
        lines = [
            f"trimmed to [p{self.lower_percentile:g}, p{self.upper_percentile:g}] within "
            f"{', '.join(self.scope)}",
            f"  {self.groups} group(s) x {self.features} feature(s) over {self.rows} row(s)",
            f"  {self.values_trimmed:,} of {self.values_present:,} values removed "
            f"({self.fraction_trimmed:.2%})",
        ]
        if self.empty_cells:
            shown = ", ".join(self.features_with_empty_cells[:3])
            more = (
                f", +{len(self.features_with_empty_cells) - 3} more"
                if len(self.features_with_empty_cells) > 3
                else ""
            )
            lines.append(
                f"  {self.empty_cells} (feature, group) cell(s) hold no data, across "
                f"{len(self.features_with_empty_cells)} feature(s): {shown}{more}"
            )
            lines.append(
                "  those features cannot yield a complete concentration series and will "
                "not be fitted"
            )
        if self.rows_with_missing_scope:
            lines.append(
                f"  {self.rows_with_missing_scope} row(s) have a missing value in a scope "
                f"column and were trimmed together as one group"
            )
        return "\n".join(lines)


def trim_extremes(
    frame: pd.DataFrame,
    features: Sequence[str],
    scope: Sequence[str],
    lower_percentile: float,
    upper_percentile: float,
    *,
    inplace: bool = False,
    block_size: int = BLOCK_SIZE,
) -> tuple[pd.DataFrame, TrimReport]:
    """Replace out-of-interval feature values with NaN, group by group.

    Percentiles are computed per (group, feature) over that cell's non-missing values,
    matching `numpy.nanpercentile` with linear interpolation. Values outside
    `[lower, upper]` become NaN; everything else, including the metadata columns, is
    untouched.

    Pass `inplace=True` to write into `frame` rather than a copy - worth it on the full
    per-object table, where the copy alone is gigabytes.
    """
    import numpy as np
    import pandas as pd

    features = list(features)
    scope = list(scope)
    missing = [column for column in (*features, *scope) if column not in frame.columns]
    if missing:
        raise KeyError(f"columns absent from the table: {', '.join(missing)}")
    if not 0.0 <= lower_percentile < upper_percentile <= 100.0:
        raise ValueError(
            f"need 0 <= lower < upper <= 100, got lower={lower_percentile!r} "
            f"upper={upper_percentile!r}"
        )

    if not inplace:
        frame = frame.copy()

    # dropna=False keeps rows whose scope value is missing: with the default they would
    # be excluded from the grouped quantiles and so silently escape trimming altogether.
    codes = frame.groupby(scope, sort=False, observed=True, dropna=False).ngroup()
    n_groups = int(codes.max()) + 1 if len(codes) else 0
    code_values = codes.to_numpy()

    trimmed_per_feature = np.zeros(len(features), dtype=np.int64)
    present_total = 0
    empty_cells = 0
    features_with_empty: list[str] = []

    for start in range(0, len(features), block_size):
        block = features[start : start + block_size]
        subset = frame[block]
        grouped = subset.groupby(codes, sort=False)
        low = grouped.transform("quantile", lower_percentile / 100.0).to_numpy()
        high = grouped.transform("quantile", upper_percentile / 100.0).to_numpy()
        values = subset.to_numpy(dtype=np.float64)

        present = ~np.isnan(values)
        # NaN never satisfies either comparison, so missing values are left as they are.
        drop = (values < low) | (values > high)

        present_total += int(present.sum())
        trimmed_per_feature[start : start + len(block)] = drop.sum(axis=0)

        for offset, name in enumerate(block):
            per_group = np.bincount(
                code_values, weights=present[:, offset], minlength=n_groups
            )
            empty = int((per_group == 0).sum())
            if empty:
                empty_cells += empty
                features_with_empty.append(name)

        frame[block] = np.where(drop, np.nan, values)

    report = TrimReport(
        scope=tuple(scope),
        lower_percentile=lower_percentile,
        upper_percentile=upper_percentile,
        groups=n_groups,
        features=len(features),
        rows=len(frame),
        values_present=present_total,
        values_trimmed=int(trimmed_per_feature.sum()),
        empty_cells=empty_cells,
        features_with_empty_cells=tuple(features_with_empty),
        rows_with_missing_scope=int(frame[scope].isna().any(axis=1).sum()),
        per_feature=pd.Series(trimmed_per_feature, index=features, name="values_trimmed"),
    )
    return frame, report
