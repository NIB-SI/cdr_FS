"""Configuration fixtures for the reproduction gates.

`write_config` renders the published design against the committed fixture, so a gate can
state its own departures from it as overrides rather than carrying a whole `.ini`.
"""

from __future__ import annotations

from pathlib import Path

FIXTURES = Path(__file__).parent / "fixtures"
#: 1,272 cells x 30 columns: the 10 metadata columns, both organelle-count features and
#: 18 rp_norm_* features, drawn from the published table. See fixtures/README.md.
SUBSET = FIXTURES / "subset.tsv"
#: The 481 column names of the full published table, one per line.
COLUMNS_PUBLISHED = FIXTURES / "columns_published.txt"

METADATA_PATTERNS = [
    "^Metadata_",
    "^Concentration$",
    "^Tech_replica$",
    "^Day_Well_BR$",
    "^cell_ID$",
    "^counts_Cells$",
    "^counts_Cytoplasm$",
    "^counts_FilteredNuclei$",
]

#: A valid configuration for the subset fixture, mirroring examples/published.ini.
BASE: dict[str, dict[str, str]] = {
    "input": {"table": "<table>", "sep": "tab"},
    "schema": {
        "metadata_patterns": "\n".join(METADATA_PATTERNS),
        "condition": "Concentration",
        "group_by": "Metadata_Day",
        "pool_over": "Metadata_Biorep",
    },
    "design": {
        "control": "11",
        "levels": "10,9,8,7,6,5,4,3,2",
        "dose": "",
        "exclude_from_fit": "2",
    },
    "trim": {
        "enabled": "true",
        "lower_percentile": "2.5",
        "upper_percentile": "97.5",
        "scope": "Metadata_Day,Metadata_Biorep,Metadata_Well",
    },
    "emd": {},
    "fit": {},
    "select": {"strata": "D1,D5,D7,D9"},
    "correlation": {"enabled": "true"},
    "missing_data": {"drop_missing": "false"},
    "output": {"dir": "<output>"},
}


def config_text(
    overrides: dict[str, str | None] | None = None,
    *,
    table: Path | str = SUBSET,
    output: Path | str = "results",
    base: dict[str, dict[str, str]] | None = None,
) -> str:
    """Render a configuration, applying dotted-key `overrides`.

    An override value of `None` removes the key, or the whole section when the key is a
    bare section name. `"section.key": "value"` adds or replaces.
    """
    sections = {name: dict(keys) for name, keys in (base or BASE).items()}
    sections["input"]["table"] = str(table)
    sections["output"]["dir"] = str(output)

    for target, value in (overrides or {}).items():
        section, _, key = target.partition(".")
        if not key:
            if value is None:
                sections.pop(section, None)
            else:
                sections.setdefault(section, {})
            continue
        if value is None:
            sections.get(section, {}).pop(key, None)
        else:
            sections.setdefault(section, {})[key] = value

    lines = []
    for name, keys in sections.items():
        lines.append(f"[{name}]")
        for key, value in keys.items():
            if "\n" in value:
                lines.append(f"{key} =")
                lines.extend(f"    {part}" for part in value.splitlines())
            else:
                lines.append(f"{key} = {value}")
        lines.append("")
    return "\n".join(lines)


def write_config(
    directory: Path,
    overrides: dict[str, str | None] | None = None,
    *,
    name: str = "config.ini",
    table: Path | str = SUBSET,
    base: dict[str, dict[str, str]] | None = None,
) -> Path:
    path = Path(directory) / name
    path.write_text(
        config_text(overrides, table=table, output=Path(directory) / "results", base=base),
        encoding="utf-8",
    )
    return path


def validate(path: Path):
    """Everything `cdr-fs check` validates from the file and the header, as one call.

    Raises `ConfigError` for the first problem found. The checks that need the data
    itself live in `Config.validate_observed` and are exercised directly.
    """
    from cdr_fs.config import load_config
    from cdr_fs.schema import read_header

    config = load_config(path)
    warnings = config.validate_columns(read_header(config.input.table, config.input.sep))
    return config, warnings
