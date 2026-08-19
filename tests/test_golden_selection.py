"""The retention rule against the published feature lists.

The sharpest test in the suite, and the cheapest: given the fit table the article's run
produced, the rule either returns the published set of features or it does not. No
tolerance, no numerical drift, no large input - everything it needs is committed, so it runs
everywhere.

The fit table **predates the BC4 correction**, which is why its filename says so: back then
BC4 was implemented as the four-parameter log-logistic, so its AIC is bit-identical to LL4's
in all 1,760 series where both converged. That is exactly right here - it is the article's
fit, and these are the article's lists, and together they are what yields 182. Swapping in a
fit table from the current pipeline would give 199 and this test would fail, correctly;
`test_the_fixture_is_the_pre_correction_fit` guards against exactly that mix-up.

The lists were taken from the headers of the published `*_trimmed_features.txt` tables, which
are those lists applied back to the data. See `fixtures/README.md`.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from cases import write_config
from tests_support import as_fit_table, feature_list

from cdr_fs.config import load_config
from cdr_fs.select import select_features

FIXTURES = Path(__file__).resolve().parent / "fixtures"
PUBLISHED_FITS = FIXTURES / "published_model_fit_results_pre_bc4_fix.txt"

#: (strata, fixture, expected count). 182 is the figure the article reports.
GATES = [
    pytest.param("D1,D5,D7,D9", "published_selected_all_days.txt", 182, id="all-days"),
    pytest.param("D5", "published_selected_D5.txt", 374, id="D5-only"),
]


@pytest.mark.parametrize(("strata", "fixture", "expected"), GATES)
def test_retention_rule_reproduces_the_published_list(tmp_path, strata, fixture, expected):
    config = load_config(write_config(tmp_path, {"select.strata": strata}))
    retained, _, report = select_features(config, as_fit_table(PUBLISHED_FITS))
    published = feature_list(FIXTURES / fixture)

    assert len(published) == expected
    assert set(retained) == set(published), (
        f"{len(set(retained) - set(published))} extra, "
        f"{len(set(published) - set(retained))} missing"
    )
    assert report.retained == expected


def test_the_fixture_is_the_pre_correction_fit():
    """Guard the provenance of the fixture, because it is easy to get backwards.

    Before the correction, BC4 was `c + (d-c)/(1+exp(...))` - algebraically LL4 - so the two
    columns agree exactly wherever both converged. If someone replaces this file with a fit
    from the current pipeline, the 182 above stops being reproducible and this says why.
    """
    table = as_fit_table(PUBLISHED_FITS)
    wide = table.pivot_table(
        index=["feature", "stratum"], columns="model", values="aic", aggfunc="first"
    )
    both = wide[["BC4", "LL4"]].dropna()
    assert len(both) == 1760
    assert (both["BC4"] == both["LL4"]).all()
