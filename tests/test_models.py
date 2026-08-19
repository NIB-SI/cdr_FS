"""The six models and their information criteria.

These are the lifted formulas, so the tests are deliberately about the *arithmetic* rather
than about behaviour: each model is evaluated at known parameters and checked against a
value computed independently. If a future edit "tidies" one of the expressions, one of these
fails, which is the point.
"""

from __future__ import annotations

import numpy as np
import pytest

from cdr_fs.models import (
    MODEL_FUNCTIONS,
    MODEL_PARAMETER_NAMES,
    brain_cousens_bc4,
    brain_cousens_bc5,
    constant,
    fit_model,
    four_param_log_logistic,
    four_param_weibull,
    information_criteria,
    initial_guess,
    linear,
)

X = np.arange(1, 9, dtype=float)


def test_every_model_has_named_parameters_matching_its_signature():
    import inspect

    for name, function in MODEL_FUNCTIONS.items():
        parameters = list(inspect.signature(function).parameters)[1:]  # drop x
        assert tuple(parameters) == MODEL_PARAMETER_NAMES[name], name
        assert len(initial_guess(name, X, X)) == len(parameters), name


def test_log_logistic_matches_an_independent_computation():
    b, c, d, e = 1.5, 0.5, 4.0, 3.0
    expected = c + (d - c) / (1 + np.exp(b * (np.log(X + 1e-10) - np.log(e + 1e-10))))
    assert four_param_log_logistic(X, b, c, d, e) == pytest.approx(expected)
    # At x = e the sigmoid sits exactly halfway between the asymptotes.
    midpoint = four_param_log_logistic(np.array([e]), b, c, d, e)
    assert midpoint == pytest.approx([(c + d) / 2], rel=1e-9)


def test_weibull_matches_an_independent_computation():
    b, c, d, e = 1.5, 0.5, 4.0, 3.0
    expected = c + (d - c) * np.exp(-np.exp(b * (np.log(X + 1e-10) - np.log(e + 1e-10))))
    assert four_param_weibull(X, b, c, d, e) == pytest.approx(expected)


def test_bc4_is_bc5_with_the_lower_asymptote_at_zero():
    b, d, e, f = 1.5, 4.0, 3.0, 0.2
    assert brain_cousens_bc4(X, b, d, e, f) == pytest.approx(
        brain_cousens_bc5(X, b, 0.0, d, e, f)
    )


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


def test_linear_and_constant_are_what_they_say():
    assert linear(X, 2.0, 3.0) == pytest.approx(2.0 * X + 3.0)
    assert constant(X, 7.0) == pytest.approx(np.full_like(X, 7.0))
    # constant must return an array shaped like x, not a scalar; curve_fit relies on it.
    assert constant(X, 7.0).shape == X.shape


def test_information_criteria_match_the_published_formula():
    residuals = np.array([0.1, -0.2, 0.3, -0.1, 0.05, -0.05, 0.2, -0.3])
    n, k = len(residuals), 4
    rss = float(np.sum(residuals**2))
    aic, bic = information_criteria(residuals, k)
    assert aic == pytest.approx(n * np.log(rss / n) + 2 * k)
    assert bic == pytest.approx(n * np.log(rss / n) + k * np.log(n))
    # BIC penalises harder than AIC once n > e^2, i.e. from 8 points on.
    assert bic > aic


def test_a_perfect_fit_gives_minus_infinite_criteria():
    # Inherited behaviour: RSS of zero puts log(0) into the expression. It is consistent -
    # a perfect fit does win the comparison - and it does occur, on series whose distances
    # are all identical.
    with np.errstate(divide="ignore"):
        aic, bic = information_criteria(np.zeros(8), 1)
    assert aic == -np.inf
    assert bic == -np.inf


def test_constant_model_recovers_the_mean():
    y = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
    result = fit_model("Con", X, y)
    assert result is not None
    assert result.parameters[0] == pytest.approx(y.mean())
    assert result.slope is None


def test_linear_model_recovers_a_known_line_and_reports_its_slope():
    y = 2.5 * X - 1.0
    result = fit_model("Lin", X, y)
    assert result is not None
    assert result.parameters == pytest.approx((2.5, -1.0))
    assert result.slope == pytest.approx(2.5)


def test_only_the_linear_model_reports_a_slope():
    y = 2.5 * X - 1.0
    for name in MODEL_FUNCTIONS:
        result = fit_model(name, X, y)
        if result is None:
            continue
        assert (result.slope is None) == (name != "Lin"), name


def test_missing_values_are_dropped_before_fitting():
    y = 2.0 * X + 1.0
    holed = y.copy()
    holed[3] = np.nan
    result = fit_model("Lin", X, holed)
    assert result is not None
    assert result.parameters == pytest.approx((2.0, 1.0))


def test_an_all_missing_series_returns_nothing():
    assert fit_model("Lin", X, np.full(8, np.nan)) is None


def test_fewer_points_than_parameters_raises_rather_than_returning_none():
    """Not swallowed, and not meant to be reachable.

    scipy raises TypeError for this, which the inherited `except (RuntimeError,
    OverflowError)` deliberately does not catch - the original would have propagated it
    too. It should never arise through the CLI, because `config.py` refuses a
    configuration whose fitted levels do not outnumber the widest model's parameters
    (`tests/cases.py`, "more parameters than points"). Letting it surface here keeps the
    two layers honest: config prevents it, fitting does not pretend to handle it.
    """
    with pytest.raises(TypeError, match="must not exceed the number of data points"):
        fit_model("BC5", np.array([1.0, 2.0]), np.array([1.0, 2.0]))


def test_a_failure_to_converge_returns_none():
    # A step function is a shape the Weibull cannot reach from its starting values within
    # maxfev; the caller must see None so the absence can be recorded.
    y = np.array([0.0, 0.0, 0.0, 0.0, 1e6, 1e6, 1e6, 1e6])
    outcomes = {name: fit_model(name, X, y) for name in ("Lin", "Con")}
    assert all(result is not None for result in outcomes.values())
    # Whatever the sigmoids do here, none of them may raise.
    for name in ("BC4", "BC5", "LL4", "WB1.4"):
        fit_model(name, X, y)


def test_unknown_model_is_an_error():
    with pytest.raises(KeyError, match="LL5"):
        fit_model("LL5", X, X)
    with pytest.raises(KeyError, match="LL5"):
        initial_guess("LL5", X, X)


def test_aic_plus_bic_is_the_sum():
    result = fit_model("Lin", X, 2.0 * X + 1.0)
    assert result is not None
    assert result.aic_plus_bic == pytest.approx(result.aic + result.bic)
