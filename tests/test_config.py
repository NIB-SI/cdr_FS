"""Configuration validation.

The bulk of this file is a sweep over `cases.INVALID_CASES`: `.ini` supplies no types, so
the guarantee `config.py` offers is that a wrong configuration is refused with a message
naming the section, the key and the fix. Each case asserts the message actually says that.

The second half locks the invariant the whole extraction rests on: **the defaults
reproduce the published run.** If one of those assertions has to change, the golden test
downstream is no longer testing the paper.
"""

from __future__ import annotations

import pytest
from cases import INVALID_CASES, validate, write_config

from cdr_fs.config import ConfigError, MODELS, load_config


@pytest.mark.parametrize(
    ("overrides", "fragment"),
    [pytest.param(o, f, id=label) for label, o, f in INVALID_CASES],
)
def test_rejects_misconfiguration(tmp_path, overrides, fragment):
    with pytest.raises(ConfigError) as raised:
        validate(write_config(tmp_path, overrides))
    assert fragment in str(raised.value)


def test_accepts_the_reference_configuration(tmp_path):
    """The counterweight to fifty rejections: validation must still accept the paper.

    Without it every new rule can only make `config.py` stricter, and the run this package
    exists to reproduce is the one that eventually stops loading.
    """
    config, warnings = validate(write_config(tmp_path))
    assert config.input.sep == "\t"
    assert warnings == []


def test_rejects_a_design_the_data_does_not_contain(tmp_path):
    """The three checks that need the data rather than only the header.

    They cannot run at load time - the values exist only once a column has been read - so
    they are called by whichever stage reads the table first, and by `check --scan`.
    """
    config = load_config(write_config(tmp_path))

    with pytest.raises(ConfigError, match="do not occur in column Concentration: 4"):
        config.validate_observed(levels={"11", "10", "9", "8", "7", "6", "5", "3", "2"})
    with pytest.raises(ConfigError, match="do not occur in column Metadata_Day: D5"):
        config.validate_observed(strata={"D1", "D7", "D9"})
    # A reproducibility floor is every pair of replicates, so one replicate gives no pairs.
    with pytest.raises(ConfigError, match="at least two replicates"):
        config.validate_observed(replicates={"BR1"})
    assert config.validate_observed(replicates={"BR1", "BR2"}) == []


# ------------------------------------------------------- defaults reproduce the paper


@pytest.fixture
def published(tmp_path):
    """The reference configuration with every optional key left to its default."""
    return load_config(
        write_config(
            tmp_path,
            {"select.strata": None, "design.dose": None},
        )
    )


def test_defaults_fit_eight_contrasts_of_nine(published):
    # EMD is computed for the whole series; only the fit withholds the top dose, which
    # is what keeps selection within sub-cytotoxic exposure levels.
    assert [f"{a}v{b}" for a, b in published.emd.contrasts] == [
        "11v10", "11v9", "11v8", "11v7", "11v6", "11v5", "11v4", "11v3", "11v2",
    ]
    assert published.design.fitted_levels == ("10", "9", "8", "7", "6", "5", "4", "3")


def test_defaults_are_the_published_choices(published):
    assert published.fit.models == MODELS
    assert published.fit.x_scale == "rank"
    assert published.fit.rank_by == "aic_plus_bic"
    # The published retention rule is a hybrid, and the two halves differ.
    assert (published.select.slope_positive, published.select.nonconstant) == ("any", "all")
    assert published.emd.baseline == "control_across_replicates"
    assert published.emd.per_replicate is False  # pooled replicates: 8 points per curve
    assert published.prune.threshold == 0.9
    assert published.prune.linkage == "average"
    assert published.prune.representative == "alphabetical"
    assert published.prune.fill_missing == "column_mean"
    assert published.trim.enabled is True
    assert (published.trim.lower_percentile, published.trim.upper_percentile) == (2.5, 97.5)
