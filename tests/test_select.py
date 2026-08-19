"""The retention rule.

Fit tables are written by hand here, so each case states exactly the situation it is about:
a feature that rises on one stratum and not another, a constant fit that ties for best, a
slope that points the wrong way. That is the only way to pin the two quantifiers, which
differ from each other and were previously documented wrongly in both the pipeline README
and its config comments.
"""

from __future__ import annotations

import pandas as pd
from cases import write_config

from cdr_fs.config import load_config
from cdr_fs.fit import COLUMNS as FIT_COLUMNS
from cdr_fs.select import select_features

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


def test_a_tie_counts_against_the_feature(tmp_path):
    # The original compared with `==`, so a constant fit that merely ties is enough to
    # drop the feature: the extra parameters bought nothing.
    tie = {"LL4": -50.0, "Lin": (-40.0, 1.0), "Con": -50.0}
    table = fit_table({("f", stratum): tie for stratum in STRATA})
    retained, evidence, _ = select_features(config_for(tmp_path), table)
    assert retained == []
    assert not evidence["nonconstant"].any()


def test_all_means_the_strata_the_feature_was_fitted_on(tmp_path):
    """Published behaviour, and worth stating: `all` does not mean all declared strata.

    `partial` was only fitted on D1, and passes on D1, so it is retained even though three
    strata have no evidence for it. The original looped over each feature's own rows.
    """
    entries = {("partial", "D1"): rising()}
    entries.update({("whole", stratum): rising() for stratum in STRATA})
    retained, _, _ = select_features(config_for(tmp_path), fit_table(entries))
    assert set(retained) == {"partial", "whole"}
