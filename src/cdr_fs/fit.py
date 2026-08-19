"""Fitting the distance series, producing a table.

This is the half of the original `plots_emd_model_drc.py` that was entangled with
plotting: one function fitted, appended to a shared `results` list and drew into global
matplotlib state, so neither the fitting nor the drawing could be tested. Here fitting
produces a table and `plots` reads it.

For each feature and stratum the distances are laid out along the exposure axis in the order
`[design] levels` declares, the levels in `exclude_from_fit` are dropped, and the six models
are fitted to what remains. The x-axis follows `[fit] x_scale`:

* `rank` - the position of the level in the series, 0, 1, 2, ... This is what the published
  run used, and it is the default.
* `dose` - the values from `[design] dose`. This is what yields a standard log-logistic in
  dose, because the four sigmoid models evaluate `log(x)` themselves.

**These two are not equivalent, not even for a geometric series.** It is tempting to think
they are: rank is an exact affine transform of log-dose when the dilution is constant, and
the models are log-logistic, so surely the axis only gets reparameterised. But the models
take `log(x)`, and `log(affine(rank))` is not `affine(log(rank))`. Fitted to a response that
is a perfect logistic in log-dose, `x_scale = dose` recovers it exactly while `rank` does
not - AIC -298 against -41 on an eight-point series.

Only `Lin` and `Con` are genuinely invariant, being affine in x, and the sign of the linear
slope - which is what the retention rule tests - is preserved under any increasing
reparameterisation. So the published rank axis leaves the slope half of the rule untouched
and makes the sigmoid models fit less well than a dose axis would, which makes the
"constant model is not best" half marginally harder to pass. It is the conservative choice
rather than a wrong one, and it is the default because it is what was published.

## Incomplete series are not fitted

A series missing any point is skipped entirely rather than fitted to what is left. That is
the published behaviour and it is the right one: comparing information criteria across
models only makes sense at a fixed number of points, and a series with a hole has usually
lost it for a reason that matters - a dose that killed every cell, or a feature whose
values were infinite. `FitReport` names the features this happens to, because in the
published run it silently removed two of them from consideration.

## What counts as a point

With pooled replicates each contrast contributes one point, so the published series has
eight. With `emd.per_replicate` each contrast contributes one point per replicate, so the
same series has thirty-two and the same six models are fitted to all of them - which is the
statistical reason the option exists. `n_points` is recorded on every row, so a fit is never
read without knowing how much data stood behind it.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Sequence

if TYPE_CHECKING:  # pragma: no cover
    import pandas as pd

    from cdr_fs.config import Config

__all__ = ["COLUMNS", "FitReport", "fit_series", "series_from_emd"]

#: Column order of the fit table.
COLUMNS = (
    "feature",
    "stratum",
    "model",
    "n_points",
    "n_parameters",
    "aic",
    "bic",
    "aic_plus_bic",
    "slope",
    "parameters",
)


@dataclass(frozen=True)
class FitReport:
    """What was fitted, what was not, and why."""

    features: int
    strata: tuple[str, ...]
    points_per_series: int
    series_fitted: int
    series_incomplete: int
    #: Features with at least one incomplete series. These cannot be selected.
    features_incomplete: tuple[str, ...]
    #: Features whose every stratum was incomplete, so they were never fitted at all.
    features_never_fitted: tuple[str, ...]
    #: (feature, stratum) pairs with no distance at all, so not even an incomplete series.
    series_absent: int
    #: (series, model) combinations where the optimiser did not converge.
    fits_failed: int
    models_failed: dict[str, int]

    def summary(self) -> str:
        lines = [
            f"fitted {self.series_fitted:,} series of {self.points_per_series} point(s) "
            f"over {self.features} feature(s) x {len(self.strata)} stratum/strata",
        ]
        if self.series_incomplete:
            lines.append(
                f"  {self.series_incomplete} series not fitted for want of a complete "
                f"exposure series, across {len(self.features_incomplete)} feature(s): "
                f"{_listing(self.features_incomplete)}"
            )
        if self.series_absent:
            lines.append(
                f"  {self.series_absent} (feature, stratum) pair(s) had no distance at all"
            )
        if self.features_never_fitted:
            lines.append(
                f"  {len(self.features_never_fitted)} feature(s) were never fitted on any "
                f"stratum and so cannot be selected: "
                f"{_listing(self.features_never_fitted)}"
            )
        if self.fits_failed:
            failures = ", ".join(
                f"{model} {count}" for model, count in sorted(self.models_failed.items())
            )
            lines.append(f"  {self.fits_failed} fit(s) did not converge: {failures}")
        return "\n".join(lines)


def _listing(names: Sequence[str], limit: int = 3) -> str:
    shown = ", ".join(names[:limit])
    return shown if len(names) <= limit else f"{shown}, +{len(names) - limit} more"


def _axis(config: Config) -> dict[str, float]:
    """Level label -> x position, following `[fit] x_scale`."""
    levels = config.design.fitted_levels
    if config.fit.x_scale == "rank":
        return {level: float(index) for index, level in enumerate(levels)}
    doses = config.design.dose_of
    return {level: doses[level] for level in levels}


def series_from_emd(config: Config, emd_table: pd.DataFrame):
    """Group the distance table into one (x, y) series per feature and stratum.

    Yields `(feature, stratum, x, y, complete)`. `complete` is False when any fitted level
    contributed no usable distance, in which case the series must not be fitted.
    """
    import numpy as np

    axis = _axis(config)
    wanted = {
        f"{config.design.control}v{level}": level for level in config.design.fitted_levels
    }
    subset = emd_table[emd_table["contrast"].isin(wanted)]
    if config.select.strata is not None:
        subset = subset[subset["stratum"].isin(config.select.strata)]

    levels_needed = set(config.design.fitted_levels)
    for (feature, stratum), group in subset.groupby(["feature", "stratum"], sort=False):
        levels = group["contrast"].map(wanted)
        distances = group["emd"].to_numpy(dtype=np.float64)
        usable = ~np.isnan(distances)
        complete = levels_needed <= set(levels[usable])
        x = np.array([axis[level] for level in levels], dtype=np.float64)
        # Sorted by exposure, not by whatever order the rows arrived in. The fit is
        # indifferent - each y keeps its own x either way - but a series that is not
        # monotone in x draws as a scribble, and a deterministic order is worth having.
        # Stable, so replicates at one level keep their relative order.
        order = np.argsort(x[usable], kind="stable")
        yield feature, stratum, x[usable][order], distances[usable][order], complete


def fit_series(
    config: Config, emd_table: pd.DataFrame
) -> tuple[pd.DataFrame, FitReport]:
    """Fit every configured model to every complete series."""
    import pandas as pd

    from cdr_fs.models import MODEL_PARAMETER_NAMES, fit_model

    records: list[tuple] = []
    strata: dict[str, None] = {}
    features: dict[str, None] = {}
    fitted_by_feature: dict[str, int] = {}
    incomplete_by_feature: dict[str, int] = {}
    series_fitted = 0
    series_incomplete = 0
    points = 0
    failures: dict[str, int] = {}

    for feature, stratum, x, y, complete in series_from_emd(config, emd_table):
        features.setdefault(feature)
        strata.setdefault(stratum)
        if not complete:
            series_incomplete += 1
            incomplete_by_feature[feature] = incomplete_by_feature.get(feature, 0) + 1
            continue
        series_fitted += 1
        fitted_by_feature[feature] = fitted_by_feature.get(feature, 0) + 1
        points = max(points, len(y))

        for model in config.fit.models:
            result = fit_model(model, x, y)
            if result is None:
                failures[model] = failures.get(model, 0) + 1
                continue
            records.append(
                (
                    feature,
                    stratum,
                    model,
                    len(y),
                    len(result.parameters),
                    result.aic,
                    result.bic,
                    result.aic_plus_bic,
                    result.slope,
                    ",".join(
                        f"{name}={value:.6g}"
                        for name, value in zip(
                            MODEL_PARAMETER_NAMES[model], result.parameters
                        )
                    ),
                )
            )

    # Pairs absent from the distance table entirely would otherwise be invisible: they
    # are neither fitted nor incomplete, because nothing was there to be incomplete.
    considered = [
        stratum
        for stratum in dict.fromkeys(emd_table["stratum"])
        if config.select.strata is None or stratum in config.select.strata
    ]
    expected = len(dict.fromkeys(emd_table["feature"])) * len(considered)

    table = pd.DataFrame.from_records(records, columns=COLUMNS)
    report = FitReport(
        features=len(features),
        strata=tuple(strata),
        points_per_series=points,
        series_fitted=series_fitted,
        series_incomplete=series_incomplete,
        features_incomplete=tuple(incomplete_by_feature),
        features_never_fitted=tuple(
            feature for feature in incomplete_by_feature if feature not in fitted_by_feature
        ),
        series_absent=max(0, expected - series_fitted - series_incomplete),
        fits_failed=sum(failures.values()),
        models_failed=failures,
    )
    return table, report
