"""The retention rule.

A feature is kept when its distance from the control **rises with exposure** and when a
**flat line is not already the best explanation** of that rise. Both halves are needed: the
slope test alone would keep noisy features that happen to trend upward, and the
model-comparison test alone would keep features whose distance shrinks with dose.

    slope_positive   the linear model's slope is > 0
    nonconstant      the constant model is not the best fit by [fit] rank_by

## The two quantifiers differ, and that was invisible before

Across strata the published rule is a **hybrid**: a positive slope on **any** stratum, and
the constant model beaten on **all** strata. Both the original `config.ini` comment and the
pipeline README described it as uniformly strict, which it is not - the asymmetry was an
artifact of how the two tests happened to be written, one as an inclusion list built from
all rows at once and the other as a per-day exclusion loop.

Splitting them into two keys fixes the documentation by construction and lets a user choose
`all`/`all` for the stricter rule the article's Methods section describes. The defaults
reproduce the published run.

## Ties count against a feature

The original compared the constant model's score to the best score with `==`, so a feature
whose constant fit merely *ties* for best is dropped. That is preserved: a tie means the
extra parameters bought nothing.

## Features that were never fitted

They cannot pass, because they have no linear slope to test, and they never reach this
module at all - `fit` names them instead. On the reference dataset that is what removes the
two `AreaShape_FormFactor` organelle features from contention on the strata where
infinities had emptied their populations.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Sequence

if TYPE_CHECKING:  # pragma: no cover
    import pandas as pd

    from cdr_fs.config import Config

__all__ = ["EVIDENCE_COLUMNS", "SelectReport", "select_features"]

#: Column order of the per-stratum evidence table.
EVIDENCE_COLUMNS = (
    "feature",
    "stratum",
    "slope",
    "slope_positive",
    "best_model",
    "best_score",
    "constant_score",
    "nonconstant",
)

#: `[fit] rank_by` -> the column of the fit table it names.
SCORE_COLUMNS = {"aic_plus_bic": "aic_plus_bic", "aic": "aic", "bic": "bic"}


@dataclass(frozen=True)
class SelectReport:
    """Why each feature was kept or dropped."""

    strata: tuple[str, ...]
    score: str
    slope_quantifier: str
    nonconstant_quantifier: str
    candidates: int
    retained: int
    #: Failed the slope test: distance flat or falling with exposure.
    rejected_slope: int
    #: Passed the slope test but a constant fit was as good or better.
    rejected_constant: int

    def summary(self) -> str:
        return "\n".join(
            [
                f"retained {self.retained:,} of {self.candidates:,} feature(s)",
                f"  rule: positive slope on {self.slope_quantifier} stratum/strata, "
                f"constant model beaten on {self.nonconstant_quantifier}, "
                f"ranked by {self.score}",
                # An unstratified experiment has one stratum whose label is "", which
                # would otherwise print as though the line had lost its value.
                f"  strata: {', '.join(self.strata) or '(one, unstratified)'}",
                f"  rejected: {self.rejected_slope:,} for slope, "
                f"{self.rejected_constant:,} for being no better than constant",
            ]
        )


def _quantify(flags: pd.Series, quantifier: str) -> bool:
    """Apply `any`/`all` to one feature's per-stratum verdicts.

    `all` means every stratum **this feature was fitted on**, not every stratum in
    `[select] strata`. That is the published behaviour - the original looped over the days
    present in each feature's own rows - and it is worth being explicit about, because it
    makes `all` slightly *weaker* for a feature with an incomplete series: fewer strata
    means fewer opportunities to fail. A feature fitted on no stratum at all has no rows
    here and so never reaches this function.
    """
    if quantifier == "any":
        return bool(flags.any())
    return bool(flags.all())


def select_features(
    config: Config, fit_table: pd.DataFrame
) -> tuple[list[str], pd.DataFrame, SelectReport]:
    """Apply the retention rule to a fit table.

    Returns the retained feature names in table order, a per-stratum evidence table, and a
    report.
    """
    import numpy as np
    import pandas as pd

    score_column = SCORE_COLUMNS[config.fit.rank_by]
    table = fit_table
    if config.select.strata is not None:
        table = table[table["stratum"].isin(config.select.strata)]
    strata = tuple(dict.fromkeys(table["stratum"]))

    # Per (feature, stratum): the best score of any model, and the constant model's score.
    grouped = table.groupby(["feature", "stratum"], sort=False)
    best = grouped[score_column].min().rename("best_score")
    best_model = table.loc[grouped[score_column].idxmin(), "model"]
    best_model.index = best.index
    constant = (
        table[table["model"] == "Con"]
        .set_index(["feature", "stratum"])[score_column]
        .rename("constant_score")
    )
    slope = (
        table[table["model"] == "Lin"]
        .set_index(["feature", "stratum"])["slope"]
        .rename("slope")
    )

    evidence = pd.concat([best, constant, slope], axis=1)
    evidence["best_model"] = best_model
    # A tie counts as the constant model winning: the extra parameters bought nothing.
    evidence["nonconstant"] = ~(evidence["constant_score"] == evidence["best_score"])
    evidence["slope_positive"] = pd.to_numeric(evidence["slope"], errors="coerce") > 0
    evidence = evidence.reset_index()[list(EVIDENCE_COLUMNS)]

    retained: list[str] = []
    rejected_slope = 0
    rejected_constant = 0
    for feature, group in evidence.groupby("feature", sort=False):
        if not _quantify(group["slope_positive"], config.select.slope_positive):
            rejected_slope += 1
            continue
        if not _quantify(group["nonconstant"], config.select.nonconstant):
            rejected_constant += 1
            continue
        retained.append(feature)

    # Order the output by the input's feature order rather than alphabetically: it is what
    # `prune` breaks ties on, so it must not depend on how the fit table happened to sort.
    order = {feature: index for index, feature in enumerate(dict.fromkeys(fit_table["feature"]))}
    retained.sort(key=lambda feature: order[feature])

    report = SelectReport(
        strata=strata,
        score=config.fit.rank_by,
        slope_quantifier=config.select.slope_positive,
        nonconstant_quantifier=config.select.nonconstant,
        candidates=int(evidence["feature"].nunique()),
        retained=len(retained),
        rejected_slope=rejected_slope,
        rejected_constant=rejected_constant,
    )
    return retained, evidence, report
