"""The shipped example configurations must stay valid.

Three files ship. `examples/template.ini` is the annotated file a new dataset starts from, so
it has to carry every section and every switch. `examples/published.ini` is the closest thing
to a specification of the published run. `examples/quickstart.ini` is what the README tells a
first-time reader to run, unedited, against the committed fixture - so it is the one shipped
artefact whose paths must resolve as they stand.

Nothing else in the suite would notice if a rename in `config.py` left one of them behind.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from cases import COLUMNS_PUBLISHED, SUBSET, validate

from cdr_fs.config import load_config

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE_INI = ROOT / "examples" / "template.ini"
PUBLISHED_INI = ROOT / "examples" / "published.ini"
QUICKSTART_INI = ROOT / "examples" / "quickstart.ini"

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
    concentration-response gate and two survived the correlation stage - and removed them by
    name downstream. Declaring them metadata makes the same decision at the start.
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
    assert load_config(localised).drop_missing.exclude == ()


def test_the_quickstart_runs_exactly_as_shipped(monkeypatch):
    """No edits, no substitution: the README tells a reader to run this file as it is.

    Its paths are relative to the repository root, which is the only place a shipped example
    can point without knowing where someone keeps their data.
    """
    from cdr_fs.schema import read_header, resolve_schema

    monkeypatch.chdir(ROOT)
    config = load_config(QUICKSTART_INI)
    columns = read_header(config.input.table, config.input.sep)
    assert config.validate_columns(columns) == []

    resolved = resolve_schema(columns, config.schema.compiled)
    # The number the quickstart tells the reader to check on the `[columns]` line.
    assert (len(resolved.metadata), len(resolved.features)) == (10, 20)
    assert {"counts_RelateLysoCell", "counts_RelateMitoCell"} <= resolved.feature_set


def test_the_template_carries_every_section_and_switch(tmp_path):
    """It is the file a new dataset starts from, so its structure has to be the whole structure.

    A switch it leaves out is a switch nobody discovers. The keys a user must decide are
    written plain and the rest are written out commented at their defaults, which is what makes
    the file complete and keeps `cdr-fs check`'s list of defaults in effect worth reading. Its
    column names are generic on purpose, so this stops at `load_config`: there is no table the
    template's schema is supposed to match.
    """
    import configparser

    from cdr_fs.config import _ALLOWED

    text = TEMPLATE_INI.read_text(encoding="utf-8")
    text = re.sub(r"(?m)^table = .*$", f"table = {SUBSET.as_posix()}", text)
    text = re.sub(r"(?m)^dir = .*$", f"dir = {(tmp_path / 'results').as_posix()}", text)
    path = tmp_path / "template.ini"
    path.write_text(text, encoding="utf-8")

    config = load_config(path)
    assert config.design.levels == ("L1", "L2", "L3", "L4", "L5", "L6")

    parser = configparser.ConfigParser(interpolation=None, inline_comment_prefixes=(";",))
    parser.read_string(text)
    assert set(parser.sections()) == set(_ALLOWED)
    for section, (required, optional) in _ALLOWED.items():
        present = set(parser.options(section))
        assert required <= present, f"[{section}] omits required {sorted(required - present)}"
        for key in sorted(optional - present):
            assert re.search(rf"(?m)^;\s*{re.escape(key)} =", text), (
                f"[{section}] {key} is neither set nor shown commented at its default"
            )


def test_no_shipped_configuration_carries_a_local_path():
    """A committed absolute path is the one mistake a DOI'd release cannot take back.

    `template.ini` and `published.ini` use `/PATH/TO/...` placeholders; `quickstart.ini` uses
    repository-relative paths. None may name a drive, a home directory, or an absolute location.
    """
    for path in (TEMPLATE_INI, PUBLISHED_INI, QUICKSTART_INI):
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            setting = line.split(";")[0].strip()
            if not setting.startswith(("table", "dir")):
                continue
            value = setting.partition("=")[2].strip()
            assert not re.match(r"^[A-Za-z]:", value), f"{path.name}:{number} names a drive"
            assert "~" not in value, f"{path.name}:{number} names a home directory"
            assert value.startswith("/PATH/TO/") or not value.startswith("/"), (
                f"{path.name}:{number} is an absolute path: {value}"
            )
