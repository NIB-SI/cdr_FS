"""Configuration fixtures and the table of misconfigurations that must be rejected.

Kept free of pytest imports so the same cases can be replayed by a plain script when
pytest is not installed.

`INVALID_CASES` is the specification of `config.py`'s job: every entry is a way to get a
run wrong, and the expected fragment is what the user should be told about it. Adding a
validation rule means adding a case here.
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
    "prune": {"enabled": "true"},
    "output": {"dir": "<output>"},
}


def config_text(
    overrides: dict[str, str | None] | None = None,
    *,
    table: Path | str = SUBSET,
    output: Path | str = "results",
) -> str:
    """Render a configuration, applying dotted-key `overrides`.

    An override value of `None` removes the key, or the whole section when the key is a
    bare section name. `"section.key": "value"` adds or replaces.
    """
    sections = {name: dict(keys) for name, keys in BASE.items()}
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
) -> Path:
    path = Path(directory) / name
    path.write_text(
        config_text(overrides, table=table, output=Path(directory) / "results"),
        encoding="utf-8",
    )
    return path


def validate(path: Path, *, observed: bool = False):
    """Everything `cdr-fs check` validates, in the same order, as one call.

    Raises `ConfigError` for the first problem found. With `observed=True` it also reads
    the condition and group_by columns and checks the declared design against them.
    """
    from cdr_fs.config import load_config
    from cdr_fs.schema import read_header

    config = load_config(path)
    warnings = config.validate_columns(read_header(config.input.table, config.input.sep))
    if observed:
        import pandas as pd

        usecols = [config.schema.condition]
        if config.schema.group_by:
            usecols.append(config.schema.group_by)
        frame = pd.read_csv(
            config.input.table,
            sep=config.input.sep,
            usecols=usecols,
            dtype={name: str for name in usecols},
        )
        warnings += config.validate_observed(
            levels=set(frame[config.schema.condition]),
            strata=set(frame[config.schema.group_by]) if config.schema.group_by else None,
        )
    return config, warnings


# (label, overrides, fragment the message must contain)
INVALID_CASES: list[tuple[str, dict[str, str | None], str]] = [
    # --- structure -----------------------------------------------------------
    ("unknown section", {"selection": ""}, "unknown section(s): [selection]"),
    ("missing section", {"prune": None}, "missing section(s): [prune]"),
    (
        "unknown key",
        {"trim.lower_percentil": "2.5"},
        "[trim] has unknown key(s): lower_percentil",
    ),
    # --- input ---------------------------------------------------------------
    (
        "input table absent",
        {"input.table": "no/such/table.tsv"},
        "[input] table - does not exist",
    ),
    (
        "sep given as an escape",
        {"input.sep": "\\t"},
        "expected one of: tab, comma, semicolon",
    ),
    (
        "wrong separator for the file",
        {"input.sep": "comma"},
        "check [input] sep",
    ),
    # --- schema --------------------------------------------------------------
    (
        "uncompilable regex",
        {"schema.metadata_patterns": "^counts_[Cells$"},
        "is not a valid regex",
    ),
    (
        "no metadata patterns",
        {"schema.metadata_patterns": ""},
        "every column would be treated as a feature",
    ),
    (
        "condition column not declared metadata",
        {"schema.metadata_patterns": "\n".join(METADATA_PATTERNS[:1])},
        "would be treated as features",
    ),
    (
        "condition column misspelled",
        {"schema.condition": "Concentraton"},
        "are not in subset.tsv: Concentraton",
    ),
    (
        "every column matched",
        {"schema.metadata_patterns": "."},
        "leaving no features to select from",
    ),
    # --- design --------------------------------------------------------------
    (
        "control listed among the levels",
        {"design.levels": "11,10,9,8,7,6,5,4,3,2"},
        "contains the control level",
    ),
    ("one level only", {"design.levels": "10"}, "needs at least 2"),
    ("repeated level", {"design.levels": "10,9,9,8,7,6,5,4,3"}, "repeats: 9"),
    (
        "dose length mismatch",
        {"design.dose": "11.36,19.88,34.8"},
        "index-matched",
    ),
    (
        "dose not numeric",
        {"design.dose": "11.36,19.88,34.8,60.9,106.58,186.54,326.47,571.38,low"},
        "'low' is not a number",
    ),
    (
        "excluded level unknown",
        {"design.exclude_from_fit": "1"},
        "absent from [design] levels: 1",
    ),
    (
        "every level excluded",
        {"design.exclude_from_fit": "10,9,8,7,6,5,4,3,2"},
        "withholds every level",
    ),
    # --- trim ----------------------------------------------------------------
    (
        "trim enabled without percentiles",
        {"trim.lower_percentile": None},
        "[trim] lower_percentile - is required but missing",
    ),
    ("trim enabled without scope", {"trim.scope": None}, "is required when [trim] enabled"),
    (
        "inverted percentiles",
        {"trim.lower_percentile": "97.5", "trim.upper_percentile": "2.5"},
        "would be empty",
    ),
    (
        "percentile out of range",
        {"trim.upper_percentile": "197.5"},
        "expected a number in [0, 100]",
    ),
    ("non-boolean flag", {"trim.enabled": "maybe"}, "expected a boolean"),
    (
        "trim scope column absent",
        {"trim.scope": "Metadata_Day,Metadata_Plate"},
        "are not in subset.tsv: Metadata_Plate",
    ),
    # --- emd -----------------------------------------------------------------
    (
        "baseline without a replicate column",
        {"schema.pool_over": ""},
        "[schema] pool_over names no replicate column",
    ),
    (
        "per_replicate without a replicate column",
        {"schema.pool_over": "", "emd.baseline": "none", "emd.per_replicate": "true"},
        "[emd] per_replicate - is true",
    ),
    (
        "contrast naming an unknown level",
        {"emd.contrasts": "11v10,11v99"},
        "neither the keyword",
    ),
    (
        "keyword mixed with explicit contrasts",
        {"emd.contrasts": "control_vs_each,11v10"},
        "use one or the other",
    ),
    # --- fit -----------------------------------------------------------------
    ("unknown model", {"fit.models": "BC4,LL5"}, "names unknown model(s): LL5"),
    (
        "log_dose without a dose vector",
        {"fit.x_scale": "log_dose"},
        "[design] dose is empty",
    ),
    (
        "log_dose with a non-positive dose",
        {
            "fit.x_scale": "log_dose",
            "design.dose": "0,19.88,34.8,60.9,106.58,186.54,326.47,571.38,1000",
        },
        "non-positive value(s): 10=0",
    ),
    (
        "more parameters than points",
        {"design.levels": "10,9,8,7,6", "design.exclude_from_fit": "6"},
        "more points than parameters",
    ),
    ("unknown x_scale", {"fit.x_scale": "linear"}, "expected one of: rank, dose, log_dose"),
    # --- select --------------------------------------------------------------
    (
        "slope rule without the linear model",
        {"fit.models": "BC4,BC5,LL4,WB1.4,Con"},
        "does not include Lin",
    ),
    (
        "nonconstant rule without the constant model",
        {"fit.models": "BC4,BC5,LL4,WB1.4,Lin"},
        "does not include Con",
    ),
    ("unknown quantifier", {"select.slope_positive": "some"}, "expected one of: any, all"),
    (
        "strata without a grouping column",
        {"schema.group_by": ""},
        "[schema] group_by is empty",
    ),
    # --- prune ---------------------------------------------------------------
    (
        "threshold above one",
        {"prune.threshold": "1.5"},
        "expected a number in (0, 1]",
    ),
    ("threshold of zero", {"prune.threshold": "0"}, "expected a number in (0, 1]"),
    ("ward linkage", {"prune.linkage": "ward"}, "expected one of: average, complete, single"),
]

# Cases that need the data, not just the header: (label, overrides, fragment).
OBSERVED_CASES: list[tuple[str, dict[str, str | None], str]] = [
    (
        "level absent from the data",
        {"design.levels": "10,9,8,7,6,5,4,3,1", "design.exclude_from_fit": "1"},
        "do not occur in column Concentration: 1",
    ),
    (
        "stratum absent from the data",
        {"select.strata": "D1,D3"},
        "do not occur in column Metadata_Day: D3",
    ),
]
