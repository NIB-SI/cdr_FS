"""The shipped example configuration must stay valid.

`examples/published.ini` is what users copy, and it is also the closest thing to a
specification of the published run. Nothing else would notice if a rename in `config.py`
left it behind, so these tests point it at the subset fixture and validate it for real.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from cases import COLUMNS_PUBLISHED, SUBSET, validate

from cdr_fs.config import load_config

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
    import pandas as pd

    config, warnings = validate(localised)
    # And against the values in the table, not only its header: the levels, days and
    # replicates the example declares have to be the ones the data actually holds.
    frame = pd.read_csv(
        config.input.table,
        sep=config.input.sep,
        usecols=["Concentration", "Metadata_Day", "Metadata_Biorep"],
        dtype=str,
    )
    warnings += config.validate_observed(
        levels=set(frame["Concentration"]),
        strata=set(frame["Metadata_Day"]),
        replicates=set(frame["Metadata_Biorep"]),
    )
    # One warning, and it is about the fixture rather than the configuration: `subset.tsv` is
    # a 30-column slice that happens to contain none of the eight object-index columns, so the
    # pattern naming them matches nothing here. On the full header it matches all eight - see
    # `test_the_object_index_columns_are_metadata`.
    assert len(warnings) == 1
    assert "Number_Object_Number" in warnings[0]
    assert config.trim.enabled is True
    assert config.select.strata == ("D1", "D5", "D7", "D9")


def test_the_object_index_columns_are_metadata(localised):
    """`Number_Object_Number` is CellProfiler's within-image object label, not a measurement.

    The published run carried all eight of those columns as features - seven passed the
    concentration-response gate and two survived correlation pruning - and removed them by name
    downstream. Declaring them metadata makes the same decision at the start.
    """
    from cdr_fs.schema import resolve_schema

    columns = COLUMNS_PUBLISHED.read_text(encoding="utf-8").split()
    resolved = resolve_schema(columns, load_config(localised).schema.compiled)
    assert len(columns) == 481
    assert (len(resolved.metadata), len(resolved.features)) == (18, 463)
    assert resolved.prefix_breakdown() == [("rp_", 461), ("counts_", 2)]
    assert not any("Number_Object_Number" in name for name in resolved.features)
    # Both organelle counts stay features: they are appended measurements, not labels.
    assert set(resolved.features_with_prefix("counts_")) == {
        "counts_RelateLysoCell",
        "counts_RelateMitoCell",
    }


def test_the_exclusion_list_ships_empty(localised):
    """The mechanism is there; the reference run needs no hand exclusions.

    With the object indices classified as metadata the rules reach the published feature count
    on their own, and the two features the published run struck out by hand survive on their
    own merits. Naming them here would take the final list from 95 to 93.
    """
    assert load_config(localised).subset.exclude == ()
