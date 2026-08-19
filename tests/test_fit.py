"""Laying the distances out along the exposure axis, and fitting them.

The distances are synthesised here rather than computed, so a series can be given an exact
shape - flat, rising, hormetic - and the fit checked against what that shape ought to
produce. `test_golden.py` covers the real numbers.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from cases import write_config

from cdr_fs.config import load_config
from cdr_fs.emd import COLUMNS as EMD_COLUMNS
from cdr_fs.fit import COLUMNS, fit_series, series_from_emd

STRATA = ("D1", "D5", "D7", "D9")
FITTED = ("10", "9", "8", "7", "6", "5", "4", "3")
#: The reference dose vector, so the `dose` axis can be exercised.
DOSE = ",".join(f"{1000 / 1.75**k:.6f}" for k in range(8, -1, -1))


def emd_table(shapes: dict[str, dict[str, np.ndarray]], include_top: bool = True):
    """Build a distance table. `shapes[feature][stratum]` is the 8-point series."""
    rows = []
    for feature, per_stratum in shapes.items():
        for stratum, values in per_stratum.items():
            for level, value in zip(FITTED, values):
                rows.append((feature, stratum, f"11v{level}", "11", level, "", 100, 100, value))
            if include_top:
                # The withheld top dose is measured but must never be fitted.
                rows.append((feature, stratum, "11v2", "11", "2", "", 100, 100, 999.0))
    return pd.DataFrame.from_records(rows, columns=list(EMD_COLUMNS))


def rising(scale: float = 1.0) -> np.ndarray:
    return np.arange(8, dtype=float) * scale


def flat(value: float = 2.0) -> np.ndarray:
    return np.full(8, value)


def config_for(tmp_path, **overrides):
    return load_config(write_config(tmp_path, overrides))


def test_series_are_laid_out_in_declared_level_order(tmp_path):
    config = config_for(tmp_path)
    table = emd_table({"f": {"D1": rising()}})
    # Shuffle the rows: the layout must come from [design] levels, not from row order.
    table = table.sample(frac=1.0, random_state=0).reset_index(drop=True)
    (feature, stratum, x, y, complete), = list(series_from_emd(config, table))
    assert (feature, stratum, complete) == ("f", "D1", True)
    assert list(x) == list(range(8))
    assert list(y) == list(rising())


def test_the_withheld_top_dose_is_not_fitted(tmp_path):
    config = config_for(tmp_path)
    (_, _, x, y, _), = list(series_from_emd(config, emd_table({"f": {"D1": rising()}})))
    assert len(y) == 8
    assert 999.0 not in set(y)  # the 11v2 distance is present but excluded from the fit


def test_rank_and_dose_axes_are_not_interchangeable(tmp_path):
    """The plan assumed they were, for a geometric series. They are not.

    The four sigmoid models evaluate log(x) internally, so `log(affine(rank))` is not
    `affine(log(rank))`. Only Lin and Con are invariant - and the *sign* of the linear
    slope, which is what the retention rule actually tests, survives either way.
    """
    series = 5.0 / (1 + np.exp(-1.2 * (np.arange(8) - 4.0)))
    table = emd_table({"f": {"D1": series}})

    by_rank, _ = fit_series(config_for(tmp_path), table)
    by_dose, _ = fit_series(
        config_for(tmp_path, **{"design.dose": DOSE, "fit.x_scale": "dose"}), table
    )

    rank_scores = by_rank.set_index("model")["aic"]
    dose_scores = by_dose.set_index("model")["aic"]
    assert rank_scores["LL4"] != pytest.approx(dose_scores["LL4"])
    # Lin and Con are affine-invariant, so their scores are identical.
    assert rank_scores["Con"] == pytest.approx(dose_scores["Con"])
    # The slope differs in magnitude, because the axis is rescaled, but not in sign.
    assert rank_scores["Lin"] != pytest.approx(dose_scores["Lin"])
    assert by_rank.set_index("model")["slope"]["Lin"] > 0
    assert by_dose.set_index("model")["slope"]["Lin"] > 0


def test_incomplete_series_is_not_fitted_and_is_named(tmp_path):
    values = rising()
    values[3] = np.nan
    table = emd_table({"whole": {"D1": rising()}, "holed": {"D1": values}})
    fits, report = fit_series(config_for(tmp_path), table)
    assert set(fits["feature"]) == {"whole"}
    assert report.series_incomplete == 1
    assert report.features_incomplete == ("holed",)
    assert report.features_never_fitted == ("holed",)
    assert "holed" in report.summary()


def test_a_missing_level_makes_the_series_incomplete(tmp_path):
    # Not a NaN distance: the row is absent altogether, as when a population was empty.
    table = emd_table({"f": {"D1": rising()}})
    table = table[table["contrast"] != "11v6"]
    (_, _, _, _, complete), = list(series_from_emd(config_for(tmp_path), table))
    assert complete is False


def test_series_absent_entirely_is_counted(tmp_path):
    # `f` has D1 and D5; `g` has only D1. The (g, D5) pair is neither fitted nor
    # incomplete - nothing was there to be incomplete - so it needs its own tally.
    table = emd_table({"f": {"D1": rising(), "D5": rising()}, "g": {"D1": rising()}})
    _, report = fit_series(config_for(tmp_path, **{"select.strata": "D1,D5"}), table)
    assert report.series_fitted == 3
    assert report.series_incomplete == 0
    assert report.series_absent == 1


def test_points_per_series_is_recorded(tmp_path):
    fits, report = fit_series(config_for(tmp_path), emd_table({"f": {"D1": rising()}}))
    assert report.points_per_series == 8
    assert set(fits["n_points"]) == {8}
    # Eight points against BC5's five parameters: the ratio must stay visible.
    assert fits.set_index("model")["n_parameters"]["BC5"] == 5


def test_per_replicate_gives_the_curve_more_points(tmp_path):
    rows = []
    for replicate in ("BR1", "BR2", "BR3", "BR4"):
        for level, value in zip(FITTED, rising()):
            rows.append(
                ("f", "D1", f"11v{level}", "11", level, replicate, 25, 25, value + 0.1)
            )
    table = pd.DataFrame.from_records(rows, columns=list(EMD_COLUMNS))
    _, report = fit_series(
        config_for(tmp_path, **{"emd.per_replicate": "true"}), table
    )
    assert report.points_per_series == 32
    assert report.series_fitted == 1


def test_only_the_configured_models_are_fitted(tmp_path):
    fits, _ = fit_series(
        config_for(tmp_path, **{"fit.models": "Lin,Con"}),
        emd_table({"f": {"D1": rising()}}),
    )
    assert set(fits["model"]) == {"Lin", "Con"}
    assert list(fits.columns) == list(COLUMNS)


def test_strata_outside_the_selection_are_ignored(tmp_path):
    table = emd_table({"f": {stratum: rising() for stratum in STRATA}})
    fits, report = fit_series(config_for(tmp_path, **{"select.strata": "D1,D9"}), table)
    assert set(fits["stratum"]) == {"D1", "D9"}
    assert report.strata == ("D1", "D9")


def test_a_flat_series_is_fitted_and_the_constant_model_wins(tmp_path):
    fits, _ = fit_series(config_for(tmp_path), emd_table({"f": {"D1": flat()}}))
    scores = fits.set_index("model")["aic_plus_bic"]
    assert scores["Con"] == scores.min()


def test_a_rising_series_beats_the_constant_model(tmp_path):
    fits, _ = fit_series(config_for(tmp_path), emd_table({"f": {"D1": rising()}}))
    scores = fits.set_index("model")["aic_plus_bic"]
    assert scores["Con"] > scores.min()
    assert fits.set_index("model")["slope"]["Lin"] > 0


def test_a_falling_series_has_a_negative_slope(tmp_path):
    fits, _ = fit_series(config_for(tmp_path), emd_table({"f": {"D1": rising()[::-1]}}))
    assert fits.set_index("model")["slope"]["Lin"] < 0


def test_parameters_are_recorded_by_name(tmp_path):
    fits, _ = fit_series(config_for(tmp_path), emd_table({"f": {"D1": rising()}}))
    linear = fits.set_index("model")["parameters"]["Lin"]
    assert linear.startswith("m=")
    assert "b=" in linear
