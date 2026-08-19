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
from cases import (
    INVALID_CASES,
    METADATA_PATTERNS,
    OBSERVED_CASES,
    validate,
    write_config,
)

from cdr_fs.config import ConfigError, MODELS, load_config


@pytest.mark.parametrize(
    ("overrides", "fragment"),
    [pytest.param(o, f, id=label) for label, o, f in INVALID_CASES],
)
def test_rejects_misconfiguration(tmp_path, overrides, fragment):
    with pytest.raises(ConfigError) as raised:
        validate(write_config(tmp_path, overrides))
    assert fragment in str(raised.value)


@pytest.mark.parametrize(
    ("overrides", "fragment"),
    [pytest.param(o, f, id=label) for label, o, f in OBSERVED_CASES],
)
def test_rejects_design_absent_from_the_data(tmp_path, overrides, fragment):
    with pytest.raises(ConfigError) as raised:
        validate(write_config(tmp_path, overrides), observed=True)
    assert fragment in str(raised.value)


def test_accepts_the_reference_configuration(tmp_path):
    config, warnings = validate(write_config(tmp_path), observed=True)
    assert config.input.sep == "\t"
    assert warnings == []


def test_unused_pattern_is_a_warning_not_an_error(tmp_path):
    # A pattern that matches nothing is how a metadata column silently becomes a feature,
    # but it is legitimate when one config serves several tables, so it must not be fatal.
    _, warnings = validate(
        write_config(
            tmp_path,
            {"schema.metadata_patterns": "\n".join([*METADATA_PATTERNS, "^counts_Cell$"])},
        )
    )
    assert any("^counts_Cell$" in warning for warning in warnings)


def test_undeclared_level_in_the_data_is_a_warning(tmp_path):
    # Dropping the top dose from [design] is a legitimate way to ignore it entirely.
    _, warnings = validate(
        write_config(
            tmp_path,
            {"design.levels": "10,9,8,7,6,5,4,3", "design.exclude_from_fit": "3"},
        ),
        observed=True,
    )
    assert any("will be ignored: 2" in warning for warning in warnings)


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


def test_omitted_keys_are_reported_as_defaulted(published):
    # Nothing may be silently defaulted: `cdr-fs check` prints this list.
    assert "fit.models" in published.defaulted
    assert "emd.per_replicate" in published.defaulted
