"""Shared readers for the published intermediates.

Kept apart from `cases.py` because these are about the *published outputs*, not about
building configurations, and two test modules need them.
"""

from __future__ import annotations

from pathlib import Path

METADATA = {
    "Concentration",
    "counts_Cells",
    "counts_Cytoplasm",
    "counts_FilteredNuclei",
    "Metadata_Well",
    "Metadata_Day",
    "Metadata_Biorep",
    "Tech_replica",
    "Day_Well_BR",
    "cell_ID",
}


def as_fit_table(path: Path):
    """The published `model_fit_results.txt`, relabelled into this package's layout.

    Its BC4 rows come from before the correction, when BC4 was the log-logistic duplicate,
    so this table is the article's fit and not the current pipeline's.
    """
    import pandas as pd

    from cdr_fs.fit import COLUMNS

    published = pd.read_csv(
        path, sep="\t", dtype={"Feature": str, "Day": str, "Model": str}
    )
    table = pd.DataFrame(
        {
            "feature": published["Feature"],
            "stratum": published["Day"],
            "model": published["Model"],
            "n_points": 8,
            "n_parameters": 0,
            "aic": published["AIC"],
            "bic": published["BIC"],
            "aic_plus_bic": published["AIC"] + published["BIC"],
            "slope": published["Slope"],
            "parameters": "",
        }
    )
    return table[list(COLUMNS)]


def feature_list(path: Path) -> list[str]:
    """A committed feature-list fixture, one name per line."""
    return Path(path).read_text(encoding="utf-8").split()
