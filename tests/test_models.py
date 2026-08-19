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
    brain_cousens_bc4,
    four_param_log_logistic,
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
