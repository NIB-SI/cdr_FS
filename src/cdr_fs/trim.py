"""Optional extreme-value trimming.

Data-level quality control: within each group of `[trim] scope`, values below the lower
percentile or above the upper percentile are discarded. In the published run the scope is
one well of one replicate on one day, and the interval kept is `[p2.5, p97.5]` inclusive,
which removes badly segmented objects without reference to the exposure level.

Three properties of this step are surprising often enough to be worth stating:

1. **Trimming removes values, not rows.** An object that is an artifact in one feature is
   usually sound in the others, so each feature is trimmed independently.
2. **Trimmed values become NaN, so N differs per feature.** That is correct for
   per-feature QC, but it propagates: a feature with a missing EMD value in a stratum is
   never fitted, hence never selected. `TrimReport` therefore reports the (feature, group)
   cells left holding no data, because those are the ones that will quietly drop out.
3. **A group containing an infinity is deleted outright** - see below. It is the reason
   `TrimReport` counts infinities at all.

## Infinities

The rule is stated as *keep* `[lower, upper]`, not as *drop* outside it, and with
infinities in the data those are not each other's complement. A percentile whose
interpolation falls between two infinite order statistics is NaN, and NaN propagates:

    np.nanpercentile([1, 2, 3, 4, 5, 6, 7, 8, 9, inf], 97.5)  ->  nan

No value satisfies `value <= nan`, so the keep mask is empty and every value in that
(group, feature) cell becomes NaN - one infinity is enough. This matters here because
`AreaShape_FormFactor` on small organelle objects divides by a perimeter that can be zero,
so those columns carry `inf` and lose whole wells. It is visible in the published EMD
tables: `rp_norm_AreaShape_FormFactor_RelateMitoCell` is absent from every D7 and D9
contrast involving the top dose, because trimming had emptied its populations.

This is inherited behaviour, reproduced deliberately so that defaults reproduce the
published run. `TrimReport` reports the infinities it found so the consequence is visible
rather than mysterious.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Sequence

if TYPE_CHECKING:  # pragma: no cover
    import pandas as pd

__all__ = ["TrimReport", "trim_extremes"]

#: Features processed per pass. The published table is 503,920 rows x 471 features, so a
#: whole-table temporary would be ~1.9 GB and several are live at once; blocking bounds the
#: peak at a few hundred MB without measurably changing the run time.
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
    #: Non-NaN feature values before trimming. Infinities count as present, matching the
    #: `dropna()` of the original scripts.
    values_present: int
    values_trimmed: int
    #: (feature, group) cells holding no data *after* trimming - what will drop out of the
    #: fit. With infinities present these are usually cells trimming emptied, not cells
    #: that arrived empty; see the module docstring.
    empty_cells: int
    features_with_empty_cells: tuple[str, ...]
    #: Infinite values in the input. Any cell containing one is emptied outright.
    values_nonfinite: int
    features_with_nonfinite: tuple[str, ...]
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
        if self.values_nonfinite:
            lines.append(
                f"  {self.values_nonfinite:,} infinite value(s) in "
                f"{len(self.features_with_nonfinite)} feature(s): "
                f"{_listing(self.features_with_nonfinite)}"
            )
            lines.append(
                "  a percentile of a tail containing an infinity is NaN, so every "
                "(feature, group) cell holding one is emptied outright"
            )
        if self.empty_cells:
            lines.append(
                f"  {self.empty_cells} (feature, group) cell(s) left with no data, across "
                f"{len(self.features_with_empty_cells)} feature(s): "
                f"{_listing(self.features_with_empty_cells)}"
            )
            lines.append(
                "  those cells cannot contribute a distance, so any curve fitted to those "
                "features is short of points, or is not fitted at all"
            )
        if self.rows_with_missing_scope:
            lines.append(
                f"  {self.rows_with_missing_scope} row(s) have a missing value in a scope "
                f"column and were trimmed together as one group"
            )
        return "\n".join(lines)


def _listing(names: Sequence[str], limit: int = 3) -> str:
    shown = ", ".join(names[:limit])
    return shown if len(names) <= limit else f"{shown}, +{len(names) - limit} more"


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
    """Keep feature values inside `[lower, upper]` percentile of their group, NaN the rest.

    Percentiles are computed per (group, feature) over that cell's non-missing values,
    matching `numpy.nanpercentile` with linear interpolation. Metadata columns are
    untouched.

    The mask is built as the complement of the *keep* condition rather than as a
    drop condition, because the two differ once a percentile is NaN - see the module
    docstring on infinities. This is what makes the result identical to the per-group loop
    of the original scripts on real data.

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

    # Row indices of each group, so the exact fallback below can address one cell without
    # rescanning the table. Built once.
    order = np.argsort(code_values, kind="stable")
    starts = np.searchsorted(code_values[order], np.arange(n_groups + 1))

    trimmed_per_feature = np.zeros(len(features), dtype=np.int64)
    present_total = 0
    nonfinite_total = 0
    empty_cells = 0
    features_with_empty: list[str] = []
    features_with_nonfinite: list[str] = []

    for start in range(0, len(features), block_size):
        block = features[start : start + block_size]
        subset = frame[block]
        grouped = subset.groupby(codes, sort=False)
        # copy=True: pandas hands back a read-only view, and the fallback below writes.
        low = grouped.transform("quantile", lower_percentile / 100.0).to_numpy(copy=True)
        high = grouped.transform("quantile", upper_percentile / 100.0).to_numpy(copy=True)
        values = subset.to_numpy(dtype=np.float64)
        nonfinite = np.isinf(values)

        # pandas' grouped quantile and numpy's percentile agree on finite data but not on
        # infinities: numpy interpolates as b - (b-a)*(1-t) once t >= 0.5, so an infinite
        # order statistic yields inf - inf = NaN, where pandas yields +inf. The original
        # scripts used np.nanpercentile, and the difference decides whether a whole cell is
        # discarded, so any column carrying an infinity has its bounds recomputed exactly.
        # Only a handful of columns do, which is why the fast path is still worth having.
        for offset in np.flatnonzero(nonfinite.any(axis=0)):
            column = values[:, offset]
            for group in range(n_groups):
                rows = order[starts[group] : starts[group + 1]]
                cell = column[rows]
                cell = cell[~np.isnan(cell)]
                if cell.size:
                    with np.errstate(invalid="ignore"):
                        low[rows, offset] = np.percentile(cell, lower_percentile)
                        high[rows, offset] = np.percentile(cell, upper_percentile)
                else:
                    low[rows, offset] = high[rows, offset] = np.nan

        # Infinities count as present, matching the dropna() of the original scripts.
        present = ~np.isnan(values)
        # The complement of the *keep* condition, not a drop condition: with a NaN bound
        # the two differ, and keep is what the original scripts expressed.
        keep = (values >= low) & (values <= high)
        drop = present & ~keep
        remaining = present & keep

        present_total += int(present.sum())
        nonfinite_total += int(nonfinite.sum())
        trimmed_per_feature[start : start + len(block)] = drop.sum(axis=0)

        for offset, name in enumerate(block):
            if nonfinite[:, offset].any():
                features_with_nonfinite.append(name)
            per_group = np.bincount(
                code_values, weights=remaining[:, offset], minlength=n_groups
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
        values_nonfinite=nonfinite_total,
        features_with_nonfinite=tuple(features_with_nonfinite),
        rows_with_missing_scope=int(frame[scope].isna().any(axis=1).sum()),
        per_feature=pd.Series(trimmed_per_feature, index=features, name="values_trimmed"),
    )
    return frame, report
