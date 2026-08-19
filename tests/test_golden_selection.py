"""The retention rule against the published feature lists.

The sharpest test in the suite, and the cheapest: given the fit table the article's run
produced, the rule either returns the published set of features or it does not. No
tolerance, no numerical drift, no 3.7 GB input - just the 1 MB `model_fit_results.txt` in
`data/` and two committed lists of names.

The lists were taken from the headers of the published `*_trimmed_features.txt` tables,
which are the selected feature lists applied back to the data. See `fixtures/README.md`.

Note the fit table predates the BC4 correction, so its BC4 column is the log-logistic
duplicate BC4 used to be. That is exactly right for this test: it is the article's fit, and
these are the article's lists.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from cases import write_config
from tests_support import as_fit_table, feature_list

from cdr_fs.config import load_config
from cdr_fs.select import select_features

FIXTURES = Path(__file__).resolve().parent / "fixtures"
PUBLISHED_FITS = Path(__file__).resolve().parents[1] / "data" / "model_fit_results.txt"

pytestmark = pytest.mark.skipif(
    not PUBLISHED_FITS.exists(),
    reason=(
        f"{PUBLISHED_FITS.name} absent from data/ - it is the published fit table, "
        "1 MB, and cannot be committed"
    ),
)

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


def test_the_two_gates_are_not_nested():
    """Why the hybrid rule has to be a hybrid.

    The all-days gate is the stricter one for the non-constant test and the *looser* one for
    the slope test, because a positive slope is required on ANY stratum and all-days offers
    four chances where D5 offers one. So neither published list contains the other, and the
    organelle counts land on opposite sides: `counts_RelateLysoCell` is retained across all
    days but not on D5 alone, `counts_RelateMitoCell` the other way round.

    A uniform quantifier could not produce that. If both tests were `all`, the D5 set would
    necessarily contain the all-days set. So this asymmetry is evidence from the published
    output itself that the rule is `any` on slope and `all` on non-constant - it is not just
    a reading of the code.
    """
    all_days = set(feature_list(FIXTURES / "published_selected_all_days.txt"))
    d5_only = set(feature_list(FIXTURES / "published_selected_D5.txt"))

    assert not all_days <= d5_only
    assert not d5_only <= all_days
    assert "counts_RelateLysoCell" in all_days - d5_only
    assert "counts_RelateMitoCell" in d5_only - all_days


def test_the_published_selection_kept_an_organelle_count():
    """The `counts_` trap, seen from the far end of the pipeline.

    `counts_RelateLysoCell` is the first entry of the published retained list. A `^counts_`
    metadata pattern would have classified it as metadata and it could never have got here.
    """
    all_days = feature_list(FIXTURES / "published_selected_all_days.txt")
    assert all_days[0] == "counts_RelateLysoCell"
    assert sum(1 for name in all_days if name.startswith("rp_norm_")) == 181
    assert len(all_days) == 182


@pytest.mark.parametrize(
    ("slope", "nonconstant", "expected"),
    [("any", "any", 443), ("any", "all", 182), ("all", "any", 272), ("all", "all", 151)],
)
def test_only_the_published_quantifiers_give_the_published_count(
    tmp_path, slope, nonconstant, expected
):
    """The published pair is `any`/`all`, and nothing else lands on 182.

    Recorded so that a future change to either default has to confront the number it breaks.
    """
    config = load_config(
        write_config(
            tmp_path,
            {"select.slope_positive": slope, "select.nonconstant": nonconstant},
        )
    )
    retained, _, _ = select_features(config, as_fit_table(PUBLISHED_FITS))
    assert len(retained) == expected
