"""The six models and their information criteria.

These are the inherited formulas, so the tests exist to stop a future edit "tidying" one of
them into a different curve. That only works if the expected value is derived *independently*
of the expression under test - a test that re-types the source line passes for any rewrite
that the same hand makes to both. So each model is pinned at the one argument where its
logistic term collapses to something exact:

    at x = e, `b * (log(x) - log(e))` is 0, and `exp(0)` is 1

which leaves a closed form per model, each different from the others:

    LL4     c + (d - c)/2          the midpoint, symmetric
    WB1.4   c + (d - c)/e          lower, because the Weibull is not symmetric
    BC4     (d + f*e)/2            the hormesis term survives into the numerator
    BC5     c + (d - c + f*e)/2

Those four numbers are what separate the four sigmoids from each other. Swap the Weibull's
expression for a log-logistic and it lands on the midpoint instead of at `exp(-1)` of the
way, and the test says so.
"""

from __future__ import annotations

import numpy as np
import pytest

from cdr_fs.models import (
    brain_cousens_bc4,
    brain_cousens_bc5,
    fit_model,
    four_param_log_logistic,
    four_param_weibull,
    information_criteria,
)

X = np.arange(1, 9, dtype=float)


def test_the_hormesis_term_is_what_distinguishes_bc4_from_log_logistic():
    """The post-article correction, pinned.

    Before it, BC4 was `c + (d-c)/(1+exp(...))` - algebraically LL4, so the six-model
    comparison was really five. With `f = 0` the corrected BC4 still coincides with an LL4
    whose lower asymptote is zero; it is the `f * x` term that makes it a hormesis model.
    """
    b, d, e = 3.0, 4.0, 5.0
    assert brain_cousens_bc4(X, b, d, e, 0.0) == pytest.approx(
        four_param_log_logistic(X, b, 0.0, d, e)
    )
    with_hormesis = brain_cousens_bc4(X, b, d, e, 0.5)
    assert not np.allclose(with_hormesis, four_param_log_logistic(X, b, 0.0, d, e))
    # The hormetic rise: below the inflection the f*x term lifts the curve above the
    # asymptote it would otherwise decay from. No log-logistic can do that.
    assert with_hormesis.max() > d


def test_information_criteria_match_the_published_formula():
    residuals = np.array([0.1, -0.2, 0.3, -0.1, 0.05, -0.05, 0.2, -0.3])
    n, k = len(residuals), 4
    rss = float(np.sum(residuals**2))
    aic, bic = information_criteria(residuals, k)
    assert aic == pytest.approx(n * np.log(rss / n) + 2 * k)
    assert bic == pytest.approx(n * np.log(rss / n) + k * np.log(n))
    # BIC penalises harder than AIC once n > e^2, i.e. from 8 points on.
    assert bic > aic


def test_each_sigmoid_is_pinned_where_its_logistic_term_is_exact():
    """One closed form per model at `x = e`, derived from the shape rather than the source.

    The four sigmoids differ only in what surrounds the same logistic term, so this is the
    argument at which those differences are visible as four distinct numbers.
    """
    b, c, d, e, f = 1.5, 0.5, 4.0, 3.0, 0.25
    at_e = np.array([e])

    assert four_param_log_logistic(at_e, b, c, d, e)[0] == pytest.approx((c + d) / 2)
    assert four_param_weibull(at_e, b, c, d, e)[0] == pytest.approx(c + (d - c) * np.exp(-1))
    assert brain_cousens_bc4(at_e, b, d, e, f)[0] == pytest.approx((d + f * e) / 2)
    assert brain_cousens_bc5(at_e, b, c, d, e, f)[0] == pytest.approx(c + (d - c + f * e) / 2)

    # The Weibull's asymmetry is the whole difference between it and the log-logistic.
    assert four_param_weibull(at_e, b, c, d, e)[0] != pytest.approx((c + d) / 2)

    # Both run from d at negligible exposure to c at overwhelming exposure, for b > 0.
    for model in (four_param_log_logistic, four_param_weibull):
        assert model(np.array([1e-9]), b, c, d, e)[0] == pytest.approx(d)
        assert model(np.array([1e9]), b, c, d, e)[0] == pytest.approx(c)

    # BC4 is BC5 with the lower asymptote pinned at zero - a relation between two of the
    # expressions, so neither can be rewritten without the other following.
    grid = np.array([0.5, 1.0, 3.0, 9.0])
    assert brain_cousens_bc4(grid, b, d, e, f) == pytest.approx(
        brain_cousens_bc5(grid, b, 0.0, d, e, f)
    )


def test_a_fit_that_runs_out_of_budget_returns_none():
    """`fit_model` returns None rather than raising when the optimiser gives up.

    The published run has 541 such fits, and `FitReport` counts them, so the branch is part
    of the method's output rather than an edge case. Starving `maxfev` is the deterministic
    way to reach it: a series that defeats the optimiser on its own budget would make this
    test a hostage to the scipy version.
    """
    y = np.array([0.0, 0.0, 0.0, 0.0, 1e6, 1e6, 1e6, 1e6])
    for model in ("BC4", "BC5", "LL4", "WB1.4", "Lin"):
        assert fit_model(model, X, y, maxfev=3) is None, model
    # `Con` needs no iteration, so it survives any budget - which is what shows the None
    # above comes from the optimiser giving up rather than from the call failing.
    assert fit_model("Con", X, y, maxfev=3) is not None
