"""The six concentration-response models, and their information criteria.

**The formulas in this module are lifted verbatim from the published pipeline**
(`plots_emd_model_drc.py`, which this file is the direct descendant of - `git log --follow`
reaches it). They are the core of the method and are deliberately not rewritten,
reformulated or "improved": every parameterisation, every `1e-10` guard, the `p0` starting
values, the `maxfev` budget and the AIC/BIC expressions are exactly as they were when the
article's numbers were produced. What changed is only the surroundings: fitting now returns
a table instead of drawing into global matplotlib state.

The set spans the shapes a distance-versus-exposure series can take:

| Name | Shape | Parameters |
|---|---|---|
| `BC4` | Brain-Cousens hormesis, lower asymptote fixed at 0 | b, d, e, f |
| `BC5` | Brain-Cousens hormesis, lower asymptote free | b, c, d, e, f |
| `LL4` | four-parameter log-logistic | b, c, d, e |
| `WB1.4` | four-parameter Weibull | b, c, d, e |
| `Lin` | straight line | m, b |
| `Con` | constant | c |

`Lin` and `Con` are what the retention rule actually consults - the sign of the slope, and
whether a flat line already explains the series - so they are not optional extras. The
`f * x` term is what makes the Brain-Cousens pair hormesis models and distinguishes BC4
from the log-logistic; BC4 is BC5 with the lower asymptote fixed at 0.

## On information criteria at this sample size

With the published design the series has eight points and BC5 has five parameters. That
ratio is at the edge of what AIC/BIC comparison supports, which is exactly why
`emd.per_replicate` exists as an option - see that module. `fit` records the number of
points behind every fit so the ratio is never invisible.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

# Unlike the other modules, numpy is imported here rather than inside each function: the
# model functions below are handed to `curve_fit` and evaluated in its inner loop, so they
# have to be plain module-level functions closing over `np`, exactly as in the original.
# Nothing imports this module unless it is about to fit something.
import numpy as np

__all__ = [
    "FitResult",
    "MODEL_FUNCTIONS",
    "MODEL_PARAMETER_NAMES",
    "brain_cousens_bc4",
    "brain_cousens_bc5",
    "constant",
    "fit_model",
    "four_param_log_logistic",
    "four_param_weibull",
    "information_criteria",
    "initial_guess",
    "linear",
]

#: Optimiser budget, as in the published run.
MAX_EVALUATIONS = 50_000


# --------------------------------------------------------------- the models, verbatim
#
# Do not reformulate these. See the module docstring.


def brain_cousens_bc4(x, b, d, e, f):
    return (d + f * x) / (1 + np.exp(b * (np.log(x + 1e-10) - np.log(e + 1e-10))))


def brain_cousens_bc5(x, b, c, d, e, f):
    return c + (d - c + f * x) / (1 + np.exp(b * (np.log(x + 1e-10) - np.log(e + 1e-10))))


def four_param_log_logistic(x, b, c, d, e):
    return c + (d - c) / (1 + np.exp(b * (np.log(x + 1e-10) - np.log(e + 1e-10))))


def four_param_weibull(x, b, c, d, e):
    return c + (d - c) * np.exp(-np.exp(b * (np.log(x + 1e-10) - np.log(e + 1e-10))))


def linear(x, m, b):
    return m * x + b


def constant(x, c):
    return np.full_like(x, c, dtype=np.float64)


MODEL_FUNCTIONS: dict[str, Callable] = {
    "BC4": brain_cousens_bc4,
    "BC5": brain_cousens_bc5,
    "LL4": four_param_log_logistic,
    "WB1.4": four_param_weibull,
    "Lin": linear,
    "Con": constant,
}

MODEL_PARAMETER_NAMES: dict[str, tuple[str, ...]] = {
    "BC4": ("b", "d", "e", "f"),
    "BC5": ("b", "c", "d", "e", "f"),
    "LL4": ("b", "c", "d", "e"),
    "WB1.4": ("b", "c", "d", "e"),
    "Lin": ("m", "b"),
    "Con": ("c",),
}


# ------------------------------------------------------------------- starting values


def initial_guess(model: str, x, y) -> list[float]:
    """Starting parameters, exactly as the published run chose them.

    They are derived from the data rather than fixed, so they follow `x_scale`: `e`, the
    inflection, starts at the median of `x`, which is the middle rank under `rank` and the
    median dose under `dose`.
    """
    guesses = {
        "BC4": [1, np.max(y), np.median(x), 0.1],
        "BC5": [1, np.min(y), np.max(y), np.median(x), 0.1],
        "LL4": [1, np.min(y), np.max(y), np.median(x)],
        "WB1.4": [1, np.min(y), np.max(y), np.median(x)],
        "Lin": [1, np.mean(y)],
        "Con": [np.mean(y)],
    }
    if model not in guesses:
        raise KeyError(f"unknown model {model!r}; known: {', '.join(MODEL_FUNCTIONS)}")
    return guesses[model]


def information_criteria(residuals, n_parameters: int) -> tuple[float, float]:
    """AIC and BIC from the residual sum of squares, as the published run computed them.

    Both are the Gaussian-likelihood forms `n*log(RSS/n)` plus a parameter penalty, so they
    are comparable across the six models but are not on the same scale as the AIC of any
    other software - only differences within one series mean anything.
    """
    residual_sum = np.sum(residuals**2)
    n = len(residuals)
    aic = n * np.log(residual_sum / n) + 2 * n_parameters
    bic = n * np.log(residual_sum / n) + n_parameters * np.log(n)
    return float(aic), float(bic)


@dataclass(frozen=True)
class FitResult:
    model: str
    parameters: tuple[float, ...]
    aic: float
    bic: float
    #: Only the linear model has one; it is what the retention rule tests the sign of.
    slope: float | None

    @property
    def aic_plus_bic(self) -> float:
        return self.aic + self.bic


def fit_model(model: str, x, y, *, maxfev: int = MAX_EVALUATIONS) -> FitResult | None:
    """Fit one model to one series, or return None when the optimiser cannot.

    Missing values are dropped before fitting, matching the original. A failure to converge
    is not an error - some shapes simply do not fit some series - so it returns None and the
    caller records the absence.
    """
    from scipy.optimize import curve_fit

    function = MODEL_FUNCTIONS.get(model)
    if function is None:
        raise KeyError(f"unknown model {model!r}; known: {', '.join(MODEL_FUNCTIONS)}")

    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    present = ~np.isnan(y)
    x, y = x[present], y[present]
    if x.size == 0:
        return None

    try:
        # The optimiser probes negative `e`, and log of a negative is NaN. That is how it
        # learns to back off, and the original run relied on exactly this, so the arithmetic
        # is untouched; errstate only stops numpy narrating it a few thousand times per run.
        # "ignore" changes reporting, not results: NaN and inf still propagate as before.
        with np.errstate(invalid="ignore", over="ignore", divide="ignore"):
            parameters, _ = curve_fit(
                function, x, y, p0=initial_guess(model, x, y), maxfev=maxfev
            )
            residuals = y - function(x, *parameters)
        aic, bic = information_criteria(residuals, len(parameters))
    except (RuntimeError, OverflowError):
        # RuntimeError: least squares did not converge within maxfev.
        # OverflowError: the exponential in the sigmoid models overflowed.
        return None

    return FitResult(
        model=model,
        parameters=tuple(float(value) for value in parameters),
        aic=aic,
        bic=bic,
        slope=float(parameters[0]) if model == "Lin" else None,
    )
