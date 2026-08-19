"""The retention rule.

Fit tables are written by hand here, so each case states exactly the situation it is about:
a feature that rises on one stratum and not another, a constant fit that ties for best, a
slope that points the wrong way. That is the only way to pin the two quantifiers, which
differ from each other and were previously documented wrongly in both the pipeline README
and its config comments.
"""

from __future__ import annotations

import pandas as pd
import pytest
from cases import write_config

from cdr_fs.config import load_config
from cdr_fs.fit import COLUMNS as FIT_COLUMNS
from cdr_fs.select import EVIDENCE_COLUMNS, select_features

STRATA = ("D1", "D5", "D7", "D9")


def fit_table(entries):
    """`entries[(feature, stratum)] = {model: score}`; the score is used for aic and bic.

    A `Lin` entry may be given as `(score, slope)` to set the slope explicitly.
    """
    rows = []
    for (feature, stratum), models in entries.items():
        for model, value in models.items():
            score, slope = value if isinstance(value, tuple) else (value, None)
            if model == "Lin" and slope is None:
                slope = 1.0
            rows.append(
                (feature, stratum, model, 8, 2, score / 2, score / 2, score, slope, "")
            )
    return pd.DataFrame.from_records(rows, columns=list(FIT_COLUMNS))


def rising(slope: float = 1.0, best: float = -50.0, con: float = -10.0):
    """A stratum where a non-constant model wins and the slope points up."""
    return {"LL4": best, "Lin": (best + 5, slope), "Con": con}


def flatish(slope: float = 1.0, con: float = -50.0):
    """A stratum where the constant model is the best fit."""
    return {"LL4": -10.0, "Lin": (-9.0, slope), "Con": con}


def config_for(tmp_path, **overrides):
    return load_config(write_config(tmp_path, overrides))


def test_a_feature_rising_everywhere_is_retained(tmp_path):
    table = fit_table({("f", stratum): rising() for stratum in STRATA})
    retained, evidence, report = select_features(config_for(tmp_path), table)
    assert retained == ["f"]
    assert list(evidence.columns) == list(EVIDENCE_COLUMNS)
    assert report.retained == 1
    assert report.candidates == 1


def test_the_published_rule_is_a_hybrid(tmp_path):
    """Positive slope on ANY stratum, constant model beaten on ALL of them.

    `mixed_slope` rises on D1 and falls on the rest, and is kept, because the slope test
    only needs one stratum. `mixed_constant` rises everywhere but a flat line explains D7
    best, and is dropped, because that test needs all of them. A rule with one quantifier
    could not produce both outcomes, which is why they are two config keys.
    """
    entries = {}
    for stratum in STRATA:
        entries[("mixed_slope", stratum)] = rising(slope=1.0 if stratum == "D1" else -1.0)
        entries[("mixed_constant", stratum)] = (
            flatish() if stratum == "D7" else rising()
        )
    retained, _, report = select_features(config_for(tmp_path), fit_table(entries))
    assert retained == ["mixed_slope"]
    assert report.rejected_constant == 1


def test_all_all_is_the_stricter_rule(tmp_path):
    entries = {}
    for stratum in STRATA:
        entries[("mixed_slope", stratum)] = rising(slope=1.0 if stratum == "D1" else -1.0)
    strict = config_for(tmp_path, **{"select.slope_positive": "all"})
    retained, _, report = select_features(strict, fit_table(entries))
    assert retained == []
    assert report.rejected_slope == 1
    assert report.slope_quantifier == "all"


def test_a_tie_counts_against_the_feature(tmp_path):
    # The original compared with `==`, so a constant fit that merely ties is enough to
    # drop the feature: the extra parameters bought nothing.
    tie = {"LL4": -50.0, "Lin": (-40.0, 1.0), "Con": -50.0}
    table = fit_table({("f", stratum): tie for stratum in STRATA})
    retained, evidence, _ = select_features(config_for(tmp_path), table)
    assert retained == []
    assert not evidence["nonconstant"].any()


def test_a_negative_slope_is_rejected(tmp_path):
    table = fit_table({("f", stratum): rising(slope=-2.0) for stratum in STRATA})
    retained, _, report = select_features(config_for(tmp_path), table)
    assert retained == []
    assert report.rejected_slope == 1
    assert report.rejected_constant == 0


def test_a_zero_slope_is_rejected(tmp_path):
    # Strictly greater than zero, as published: a flat line is not a response.
    table = fit_table({("f", stratum): rising(slope=0.0) for stratum in STRATA})
    assert select_features(config_for(tmp_path), table)[0] == []


def test_all_means_the_strata_the_feature_was_fitted_on(tmp_path):
    """Published behaviour, and worth stating: `all` does not mean all declared strata.

    `partial` was only fitted on D1, and passes on D1, so it is retained even though three
    strata have no evidence for it. The original looped over each feature's own rows.
    """
    entries = {("partial", "D1"): rising()}
    entries.update({("whole", stratum): rising() for stratum in STRATA})
    retained, _, _ = select_features(config_for(tmp_path), fit_table(entries))
    assert set(retained) == {"partial", "whole"}


def test_strata_outside_the_selection_do_not_count(tmp_path):
    # D9 would sink the feature, but the configuration does not ask about D9.
    entries = {("f", stratum): rising() for stratum in ("D1", "D5", "D7")}
    entries[("f", "D9")] = flatish()
    table = fit_table(entries)
    assert select_features(config_for(tmp_path), table)[0] == []
    narrowed = config_for(tmp_path, **{"select.strata": "D1,D5,D7"})
    retained, _, report = select_features(narrowed, table)
    assert retained == ["f"]
    assert report.strata == ("D1", "D5", "D7")


def test_rank_by_changes_the_verdict(tmp_path):
    # AIC and BIC penalise parameters differently, so which model wins can differ. Here
    # LL4 wins on AIC and Con wins on BIC.
    rows = [
        ("f", "D1", "LL4", 8, 4, -30.0, -5.0, -35.0, None, ""),
        ("f", "D1", "Lin", 8, 2, -20.0, -15.0, -35.0, 1.0, ""),
        ("f", "D1", "Con", 8, 1, -10.0, -25.0, -35.0, None, ""),
    ]
    table = pd.DataFrame.from_records(rows, columns=list(FIT_COLUMNS))
    assert select_features(config_for(tmp_path, **{"fit.rank_by": "aic"}), table)[0] == ["f"]
    assert select_features(config_for(tmp_path, **{"fit.rank_by": "bic"}), table)[0] == []
    # On aic_plus_bic all three tie, and a tie counts against the feature.
    assert select_features(config_for(tmp_path), table)[0] == []


def test_evidence_records_the_winning_model(tmp_path):
    entries = {("f", "D1"): rising(), ("f", "D5"): flatish()}
    _, evidence, _ = select_features(config_for(tmp_path), fit_table(entries))
    by_stratum = evidence.set_index("stratum")
    assert by_stratum.loc["D1", "best_model"] == "LL4"
    assert by_stratum.loc["D5", "best_model"] == "Con"
    assert bool(by_stratum.loc["D1", "nonconstant"]) is True
    assert bool(by_stratum.loc["D5", "nonconstant"]) is False


def test_output_order_follows_the_input_not_the_alphabet(tmp_path):
    # `prune` keeps the alphabetically first member of each cluster, so it must do that
    # sorting itself rather than inherit an accident of ordering from here.
    entries = {}
    for feature in ("zeta", "alpha", "mu"):
        entries.update({(feature, stratum): rising() for stratum in STRATA})
    retained, _, _ = select_features(config_for(tmp_path), fit_table(entries))
    assert retained == ["zeta", "alpha", "mu"]


def test_report_counts_add_up(tmp_path):
    entries = {}
    for stratum in STRATA:
        entries[("kept", stratum)] = rising()
        entries[("no_slope", stratum)] = rising(slope=-1.0)
        entries[("flat", stratum)] = flatish()
    _, _, report = select_features(config_for(tmp_path), fit_table(entries))
    assert report.candidates == 3
    assert report.retained + report.rejected_slope + report.rejected_constant == 3
    assert (report.retained, report.rejected_slope, report.rejected_constant) == (1, 1, 1)
    assert "retained 1 of 3" in report.summary()
