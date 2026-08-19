"""The shipped example configuration must stay valid.

`examples/published.ini` is what users copy, and it is also the closest thing to a
specification of the published run. Nothing else would notice if a rename in `config.py`
left it behind, so these tests point it at the subset fixture and validate it for real.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from cases import SUBSET, validate

from cdr_fs.config import MODELS, load_config

PUBLISHED_INI = Path(__file__).resolve().parents[1] / "examples" / "published.ini"

#: A 1.75-fold serial dilution from 1000 mg/L: 1000 / 1.75**k for k = 8..0. The article's
#: SI prints this series truncated (1000, 571.42, 326.530, 186.588 ... 11.36), which is how
#: the dilution factor was pinned to exactly 1.75 rather than the 1.7502 a back-calculation
#: from those truncated endpoints suggests.
DOSE = tuple(sorted(1000 / 1.75**k for k in range(9)))


@pytest.fixture
def localised(tmp_path):
    """The example config with its /PATH/TO placeholders pointed at the fixture."""
    text = PUBLISHED_INI.read_text(encoding="utf-8")
    text = re.sub(r"(?m)^table = .*$", f"table = {SUBSET.as_posix()}", text)
    text = re.sub(r"(?m)^dir = .*$", f"dir = {(tmp_path / 'results').as_posix()}", text)
    path = tmp_path / "published.ini"
    path.write_text(text, encoding="utf-8")
    return path


def test_example_is_valid_against_real_data(localised):
    config, warnings = validate(localised, observed=True)
    assert warnings == []
    assert config.trim.enabled is True
    assert config.select.strata == ("D1", "D5", "D7", "D9")


def test_example_still_carries_the_published_defaults(localised):
    config = load_config(localised)
    assert config.fit.models == MODELS
    assert config.fit.x_scale == "rank"
    assert (config.select.slope_positive, config.select.nonconstant) == ("any", "all")
    assert config.emd.per_replicate is False
    assert config.design.control == "11"
    assert config.design.fitted_levels == ("10", "9", "8", "7", "6", "5", "4", "3")


def test_example_declares_no_real_paths():
    # House rule: every user-editable path in a committed file stays /PATH/TO/...
    text = PUBLISHED_INI.read_text(encoding="utf-8")
    for line in text.splitlines():
        if re.match(r"^(table|dir) = ", line):
            assert "/PATH/TO/" in line, line


def test_dose_vector_is_the_exact_dilution_series(localised):
    config = load_config(localised)
    # The config carries six decimal places, so compare against the series rounded the
    # same way rather than against a tolerance.
    assert config.design.dose == tuple(round(dose, 6) for dose in DOSE)
    # The lowest exposure pairs with level 10 and the highest with level 2.
    assert config.design.dose_of["10"] == 11.368302
    assert config.design.dose_of["2"] == 1000.0


def test_log_dose_is_available_now_that_doses_are_declared(localised, tmp_path):
    # Every dose is positive, so the alternative x-axis validates. Rank remains the
    # default because rank is what the published run fitted.
    text = localised.read_text(encoding="utf-8").replace("x_scale = rank", "x_scale = log_dose")
    path = tmp_path / "log_dose.ini"
    path.write_text(text, encoding="utf-8")
    assert load_config(path).fit.x_scale == "log_dose"
