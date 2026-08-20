"""Configuration loading and validation.

`.ini` gives neither types nor structure, so every value a run depends on is parsed,
range-checked and cross-checked here, before any data is touched. The contract is that a
misconfigured run fails in the first second with a message the user can act on, rather
than producing plausible-looking wrong numbers.

Validation happens in three places, by how much it costs to check:

* `load_config` - everything decidable from the file alone: types, ranges, vocabularies,
  parallel-list alignment, cross-field consistency. Instant.
* `Config.validate_columns` - the columns named in `[schema]` and `[trim]` must exist in
  the table *and* be classified as metadata. Needs the header only, so also instant.
* `Config.validate_observed` - the declared levels and strata must occur in the data.
  Needs a pass over two columns, so it is called by whichever stage first reads the
  table rather than at startup.
"""

from __future__ import annotations

import configparser
import re
from dataclasses import dataclass
from pathlib import Path

__all__ = ["Config", "ConfigError", "load_config"]


class ConfigError(Exception):
    """A configuration problem the user must fix.

    The message is written to be printed as-is: the CLI reports it without a traceback,
    because a traceback tells the user nothing they can act on.
    """


SEPARATORS = {"tab": "\t", "comma": ",", "semicolon": ";"}
MODELS = ("BC4", "BC5", "LL4", "WB1.4", "Lin", "Con")
#: Free parameters per model, used to check the fit is identifiable at all.
MODEL_PARAMS = {"BC4": 4, "BC5": 5, "LL4": 4, "WB1.4": 4, "Lin": 2, "Con": 1}
# No "log_dose": the four sigmoid models already evaluate log(x) internally, so
# handing them log(dose) would log it twice - and would be NaN for any dose below 1.
# "dose" is what yields a standard log-logistic in dose.
X_SCALES = ("rank", "dose")
RANK_BY = ("aic_plus_bic", "aic", "bic")
QUANTIFIERS = ("any", "all")
BASELINES = ("control_across_replicates", "none")
LINKAGES = ("average", "complete", "single")
REPRESENTATIVES = ("alphabetical", "first")
FILL_MISSING = ("column_mean", "none")
CONTRASTS_KEYWORD = "control_vs_each"

# Keys absent from this table are rejected, so a typo cannot silently fall back to a
# default. Per section: (required keys, optional keys).
_ALLOWED: dict[str, tuple[frozenset[str], frozenset[str]]] = {
    "input": (frozenset({"table"}), frozenset({"sep"})),
    "schema": (
        frozenset({"metadata_patterns", "condition"}),
        frozenset({"group_by", "pool_over"}),
    ),
    "design": (
        frozenset({"control", "levels"}),
        frozenset({"dose", "exclude_from_fit"}),
    ),
    "trim": (
        frozenset({"enabled"}),
        frozenset({"lower_percentile", "upper_percentile", "scope"}),
    ),
    "emd": (frozenset(), frozenset({"contrasts", "baseline", "per_replicate"})),
    "fit": (frozenset(), frozenset({"models", "x_scale", "rank_by"})),
    "select": (
        frozenset({"enabled"}),
        frozenset({"slope_positive", "nonconstant", "strata"}),
    ),
    "correlation": (
        frozenset({"enabled"}),
        frozenset(
            {"threshold", "linkage", "representative", "aggregate_by", "fill_missing"}
        ),
    ),
    "drop_missing": (frozenset({"enabled"}), frozenset({"max_missing", "exclude"})),
    "output": (frozenset({"dir"}), frozenset()),
}


# ---------------------------------------------------------------------------- specs


@dataclass(frozen=True)
class InputSpec:
    table: Path
    sep_keyword: str

    @property
    def sep(self) -> str:
        """The separator character that `sep_keyword` names."""
        return SEPARATORS[self.sep_keyword]


@dataclass(frozen=True)
class SchemaSpec:
    metadata_patterns: tuple[str, ...]
    condition: str
    group_by: str | None
    pool_over: str | None

    @property
    def compiled(self) -> tuple[re.Pattern[str], ...]:
        return tuple(re.compile(pattern) for pattern in self.metadata_patterns)


@dataclass(frozen=True)
class DesignSpec:
    control: str
    #: Ordered low -> high. This is the concentration-response axis.
    levels: tuple[str, ...]
    dose: tuple[float, ...] | None
    exclude_from_fit: frozenset[str]

    @property
    def fitted_levels(self) -> tuple[str, ...]:
        """`levels`, in order, minus the ones withheld from fitting."""
        return tuple(level for level in self.levels if level not in self.exclude_from_fit)

    @property
    def dose_of(self) -> dict[str, float]:
        """Level label -> dose, empty when no dose vector was given."""
        return {} if self.dose is None else dict(zip(self.levels, self.dose))


@dataclass(frozen=True)
class TrimSpec:
    enabled: bool
    lower_percentile: float
    upper_percentile: float
    scope: tuple[str, ...]


@dataclass(frozen=True)
class EmdSpec:
    #: (reference level, compared level) pairs, in dose order.
    contrasts: tuple[tuple[str, str], ...]
    contrasts_keyword: str | None
    baseline: str
    per_replicate: bool


@dataclass(frozen=True)
class FitSpec:
    models: tuple[str, ...]
    x_scale: str
    rank_by: str

    @property
    def max_params(self) -> int:
        return max(MODEL_PARAMS[model] for model in self.models)


@dataclass(frozen=True)
class SelectSpec:
    #: False runs no concentration-response selection at all: the filtering stages start
    #: from every feature instead of from a retained list.
    enabled: bool
    slope_positive: str
    nonconstant: str
    #: None means every stratum present in the data.
    strata: tuple[str, ...] | None


@dataclass(frozen=True)
class CorrelationSpec:
    enabled: bool
    threshold: float
    linkage: str
    representative: str
    #: Columns whose combinations define the unit correlations are computed across. Empty
    #: correlates the rows as they are.
    aggregate_by: tuple[str, ...]
    fill_missing: str


@dataclass(frozen=True)
class DropMissingSpec:
    enabled: bool
    #: Percent. A feature missing this much of the table or more is dropped.
    max_missing: float
    #: Exact feature names to leave out of the final table, whatever else says otherwise.
    exclude: tuple[str, ...]


@dataclass(frozen=True)
class OutputSpec:
    dir: Path


@dataclass(frozen=True)
class Config:
    path: Path
    input: InputSpec
    schema: SchemaSpec
    design: DesignSpec
    trim: TrimSpec
    emd: EmdSpec
    fit: FitSpec
    select: SelectSpec
    correlation: CorrelationSpec
    drop_missing: DropMissingSpec
    output: OutputSpec
    #: "section.key" for every value that came from a default rather than the file.
    defaulted: tuple[str, ...]

    # ------------------------------------------------------------------- roles

    @property
    def role_columns(self) -> dict[str, str]:
        """Column name -> the role the configuration gives it.

        Each of these must exist in the table and be classified as metadata: a feature
        column doubling as a design column would be both trimmed and fitted.
        """
        roles = {self.schema.condition: "[schema] condition"}
        if self.schema.group_by:
            roles.setdefault(self.schema.group_by, "[schema] group_by")
        if self.schema.pool_over:
            roles.setdefault(self.schema.pool_over, "[schema] pool_over")
        if self.trim.enabled:
            for column in self.trim.scope:
                roles.setdefault(column, "[trim] scope")
        if self.correlation.enabled:
            for column in self.correlation.aggregate_by:
                roles.setdefault(column, "[correlation] aggregate_by")
        return roles

    # ------------------------------------------------------------- data checks

    def validate_columns(self, columns: list[str]) -> list[str]:
        """Check the configuration against the table's header.

        Returns warnings; raises `ConfigError` for anything that would corrupt a run.
        """
        # Imported here, not at module level: schema.py imports ConfigError from this
        # module, so a module-level import in this direction would be a cycle.
        from cdr_fs.schema import resolve_schema

        missing = [column for column in self.role_columns if column not in columns]
        if missing:
            # Truncate the preview: with the wrong separator the whole header arrives as
            # a single "column", and printing it in full buries the message.
            preview = ", ".join(
                column if len(column) <= 40 else f"{column[:40]}..."
                for column in columns[:4]
            )
            raise ConfigError(
                f"{self.path}: columns named in the configuration are not in "
                f"{self.input.table.name}: {', '.join(missing)}\n"
                f"  the table has {len(columns)} column(s), starting: {preview}\n"
                f"  if that looks like one long column, check [input] sep "
                f"(currently {self.input.sep_keyword})"
            )

        resolved = resolve_schema(columns, self.schema.compiled)
        misclassified = {
            column: role
            for column, role in self.role_columns.items()
            if column in resolved.feature_set
        }
        if misclassified:
            detail = "\n".join(
                f"    {column}  ({role})" for column, role in sorted(misclassified.items())
            )
            example = sorted(misclassified)[0]
            raise ConfigError(
                f"{self.path}: [schema] metadata_patterns does not match columns that "
                f"the configuration uses as design columns, so they would be treated as "
                f"features:\n{detail}\n"
                f"  add one pattern per column, e.g. ^{re.escape(example)}$"
            )

        if not resolved.features:
            raise ConfigError(
                f"{self.path}: [schema] metadata_patterns matched every one of the "
                f"{len(columns)} columns in {self.input.table.name}, leaving no features "
                f"to select from"
            )

        warnings = []
        if resolved.unused_patterns:
            warnings.append(
                "[schema] metadata_patterns entries matching no column: "
                f"{', '.join(resolved.unused_patterns)}\n"
                "  a typo here silently turns a metadata column into a feature"
            )
        # A mistyped exclusion fails in the dangerous direction: the feature stays in the
        # output and nothing says so. Names are matched exactly, so say which ones missed.
        unmatched = [
            name for name in self.drop_missing.exclude if name not in resolved.feature_set
        ]
        if unmatched:
            warnings.append(
                f"[drop_missing] exclude names {len(unmatched)} entry/entries that are not "
                f"features of {self.input.table.name}: {', '.join(unmatched)}\n"
                "  exclusions are matched exactly, so a mistyped name excludes nothing"
            )
        return warnings

    def validate_observed(
        self,
        levels: set[str] | None = None,
        strata: set[str] | None = None,
        replicates: set[str] | None = None,
    ) -> list[str]:
        """Check the declared design against the values the data actually contains.

        Called by whichever stage first reads the table, where the values come for free.
        """
        warnings = []
        if levels is not None:
            declared = {self.design.control, *self.design.levels}
            missing = sorted(declared - levels)
            if missing:
                raise ConfigError(
                    f"{self.path}: [design] declares levels that do not occur in column "
                    f"{self.schema.condition}: {', '.join(missing)}\n"
                    f"  values present: {', '.join(sorted(levels))}"
                )
            undeclared = sorted(levels - declared)
            if undeclared:
                warnings.append(
                    f"column {self.schema.condition} holds values that [design] does not "
                    f"declare, which will be ignored: {', '.join(undeclared)}"
                )
        if strata is not None and self.select.strata is not None:
            missing = sorted(set(self.select.strata) - strata)
            if missing:
                raise ConfigError(
                    f"{self.path}: [select] strata names values that do not occur in "
                    f"column {self.schema.group_by}: {', '.join(missing)}\n"
                    f"  values present: {', '.join(sorted(strata))}"
                )
        # The baseline set is every pair of replicates, so a single replicate yields no
        # pairs at all. Caught here rather than left to the file: an empty baseline table
        # is not a reproducibility floor, and this is only decidable once the column has
        # been read, which is why `[emd] baseline` cannot check it on its own.
        if (
            replicates is not None
            and self.emd.baseline == "control_across_replicates"
            and len(replicates) < 2
        ):
            held = f"only {sorted(replicates)[0]!r}" if replicates else "no values"
            raise ConfigError(
                f"{self.path}: [emd] baseline is 'control_across_replicates', which "
                f"compares replicates with each other, but column "
                f"{self.schema.pool_over} holds {held}\n"
                f"  a reproducibility floor needs at least two replicates; set [emd] "
                f"baseline = none if this experiment has none"
            )
        return warnings


# --------------------------------------------------------------------------- reader


class _Reader:
    """Typed access to one `.ini` file, with errors naming the file, section and key."""

    def __init__(self, path: Path, parser: configparser.ConfigParser) -> None:
        self.path = path
        self.parser = parser
        self.defaulted: list[str] = []

    def fail(self, section: str, key: str | None, message: str) -> ConfigError:
        where = f"[{section}] {key}" if key else f"[{section}]"
        return ConfigError(f"{self.path}: {where} - {message}")

    def _raw(self, section: str, key: str, default: str | None) -> str:
        if not self.parser.has_option(section, key):
            if default is None:
                raise self.fail(section, key, "is required but missing")
            self.defaulted.append(f"{section}.{key}")
            return default
        return self.parser.get(section, key).strip()

    def text(self, section: str, key: str, default: str | None = None) -> str:
        value = self._raw(section, key, default)
        if not value:
            raise self.fail(section, key, "is required but empty")
        return value

    def optional_text(self, section: str, key: str) -> str | None:
        return self._raw(section, key, "") or None

    def choice(
        self,
        section: str,
        key: str,
        allowed: tuple[str, ...],
        default: str | None = None,
    ) -> str:
        value = self.text(section, key, default)
        if value not in allowed:
            raise self.fail(
                section, key, f"is {value!r}; expected one of: {', '.join(allowed)}"
            )
        return value

    def boolean(self, section: str, key: str, default: str | None = None) -> bool:
        value = self.text(section, key, default).lower()
        if value in ("true", "yes", "on", "1"):
            return True
        if value in ("false", "no", "off", "0"):
            return False
        raise self.fail(
            section, key, f"is {value!r}; expected a boolean (true/false, yes/no, 1/0)"
        )

    def number(
        self,
        section: str,
        key: str,
        default: str | None = None,
        *,
        minimum: float | None = None,
        maximum: float | None = None,
        exclude_minimum: bool = False,
        exclude_maximum: bool = False,
    ) -> float:
        raw = self.text(section, key, default)
        try:
            value = float(raw)
        except ValueError:
            raise self.fail(section, key, f"is {raw!r}; expected a number") from None
        too_low = minimum is not None and (
            value <= minimum if exclude_minimum else value < minimum
        )
        too_high = maximum is not None and (
            value >= maximum if exclude_maximum else value > maximum
        )
        if too_low or too_high:
            low = "-inf" if minimum is None else f"{minimum:g}"
            high = "inf" if maximum is None else f"{maximum:g}"
            span = (
                f"{'(' if exclude_minimum else '['}{low}, "
                f"{high}{')' if exclude_maximum else ']'}"
            )
            raise self.fail(section, key, f"is {value:g}; expected a number in {span}")
        return value

    def items(
        self, section: str, key: str, default: str | None = None
    ) -> tuple[str, ...]:
        """A comma- or newline-separated list. Empty entries are dropped."""
        raw = self._raw(section, key, default)
        parts = [part.strip() for line in raw.splitlines() for part in line.split(",")]
        return tuple(part for part in parts if part)

    def lines(
        self, section: str, key: str, default: str | None = None
    ) -> tuple[str, ...]:
        """One entry per line, so an entry may itself contain commas - used for regexes."""
        raw = self._raw(section, key, default)
        return tuple(line.strip() for line in raw.splitlines() if line.strip())

    def unique(
        self, section: str, key: str, values: tuple[str, ...]
    ) -> tuple[str, ...]:
        seen: dict[str, int] = {}
        for value in values:
            seen[value] = seen.get(value, 0) + 1
        repeated = sorted(value for value, count in seen.items() if count > 1)
        if repeated:
            raise self.fail(section, key, f"repeats: {', '.join(repeated)}")
        return values


# ---------------------------------------------------------------------------- load


def load_config(path: str | Path) -> Config:
    """Read and fully validate a `cdr_FS` configuration file.

    Raises `ConfigError` with an actionable message for any problem decidable without
    reading the data table.
    """
    path = Path(path)
    if not path.is_file():
        raise ConfigError(f"configuration file not found: {path}")

    # interpolation=None: a value may contain '%' (a percentile, a regex) without
    # configparser trying to expand it.
    # inline_comment_prefixes=(';',): the ';' trailing-comment style the example config
    # uses would otherwise become part of the value.
    parser = configparser.ConfigParser(
        interpolation=None, inline_comment_prefixes=(";",)
    )
    try:
        with path.open(encoding="utf-8") as handle:
            parser.read_file(handle)
    except UnicodeDecodeError as error:
        raise ConfigError(f"{path}: is not valid UTF-8 text ({error})") from None
    except configparser.Error as error:
        raise ConfigError(f"{path}: could not be parsed - {error}") from None

    _check_sections_and_keys(path, parser)
    reader = _Reader(path, parser)

    input_spec = _read_input(reader)
    schema_spec = _read_schema(reader)
    design_spec = _read_design(reader)
    trim_spec = _read_trim(reader)
    emd_spec = _read_emd(reader, design_spec, schema_spec)
    # Read before both, because it decides which of their cross-checks apply: a run that
    # fits nothing does not need its models to be identifiable.
    selecting = reader.boolean("select", "enabled")
    fit_spec = _read_fit(reader, design_spec, selecting=selecting)
    select_spec = _read_select(reader, schema_spec, fit_spec, enabled=selecting)
    correlation_spec = _read_correlation(reader, trim_spec)
    drop_missing_spec = _read_drop_missing(reader)
    output_spec = _read_output(reader)

    return Config(
        path=path,
        input=input_spec,
        schema=schema_spec,
        design=design_spec,
        trim=trim_spec,
        emd=emd_spec,
        fit=fit_spec,
        select=select_spec,
        correlation=correlation_spec,
        drop_missing=drop_missing_spec,
        output=output_spec,
        defaulted=tuple(reader.defaulted),
    )


def _check_sections_and_keys(path: Path, parser: configparser.ConfigParser) -> None:
    """Reject unknown sections and keys, and require every known section to be present.

    Requiring the sections is what makes a section-name typo an error instead of a
    silent fall-back to defaults.
    """
    if parser.defaults():
        raise ConfigError(
            f"{path}: a [DEFAULT] section is not supported - its keys would be injected "
            f"into every other section. Found: {', '.join(sorted(parser.defaults()))}"
        )
    present = set(parser.sections())
    unknown = sorted(present - set(_ALLOWED))
    absent = sorted(set(_ALLOWED) - present)
    # Reported together rather than one at a time: a file written against an earlier version
    # has both, and "unknown [prune]" beside "missing [correlation]" lets the mapping read
    # itself. Reporting only the unknown half leaves the user to infer the rename.
    if unknown:
        lines = [
            f"{path}: unknown section(s): {', '.join('[' + s + ']' for s in unknown)}",
            f"  known sections: {', '.join('[' + s + ']' for s in _ALLOWED)}",
        ]
        if absent:
            lines.insert(1, f"  missing: {', '.join('[' + s + ']' for s in absent)}")
        raise ConfigError("\n".join(lines))
    if absent:
        raise ConfigError(
            f"{path}: missing section(s): {', '.join('[' + s + ']' for s in absent)}\n"
            f"  every section must be present, even when all its keys are left at their "
            f"defaults, so that a mistyped section name cannot pass unnoticed"
        )
    for section, (required, optional) in _ALLOWED.items():
        known = required | optional
        stray = sorted(set(parser.options(section)) - known)
        if stray:
            raise ConfigError(
                f"{path}: [{section}] has unknown key(s): {', '.join(stray)}\n"
                f"  known keys in [{section}]: {', '.join(sorted(known))}"
            )


def _read_input(reader: _Reader) -> InputSpec:
    table = Path(reader.text("input", "table"))
    sep_keyword = reader.choice(
        "input", "sep", tuple(SEPARATORS), default="tab"
    )
    if not table.exists():
        raise reader.fail("input", "table", f"does not exist: {table}")
    if not table.is_file():
        raise reader.fail("input", "table", f"is not a file: {table}")
    return InputSpec(table=table, sep_keyword=sep_keyword)


def _read_output(reader: _Reader) -> OutputSpec:
    directory = Path(reader.text("output", "dir"))
    # Decidable from the file alone, so it belongs here. Every stage writes into this
    # directory, and a path already occupied by a file fails on the first `mkdir` - which
    # is after that stage has done its work, with a traceback instead of a message.
    if directory.exists() and not directory.is_dir():
        raise reader.fail("output", "dir", f"exists but is not a directory: {directory}")
    return OutputSpec(dir=directory)


def _read_schema(reader: _Reader) -> SchemaSpec:
    patterns = reader.lines("schema", "metadata_patterns")
    if not patterns:
        raise reader.fail(
            "schema",
            "metadata_patterns",
            "is empty; every column would be treated as a feature. List one regex per "
            "line naming the metadata columns",
        )
    for pattern in patterns:
        try:
            re.compile(pattern)
        except re.error as error:
            raise reader.fail(
                "schema", "metadata_patterns", f"{pattern!r} is not a valid regex - {error}"
            ) from None
    return SchemaSpec(
        metadata_patterns=patterns,
        condition=reader.text("schema", "condition"),
        group_by=reader.optional_text("schema", "group_by"),
        pool_over=reader.optional_text("schema", "pool_over"),
    )


def _read_design(reader: _Reader) -> DesignSpec:
    control = reader.text("design", "control")
    levels = reader.unique("design", "levels", reader.items("design", "levels"))
    if len(levels) < 2:
        raise reader.fail(
            "design",
            "levels",
            f"lists {len(levels)} level(s); a concentration-response axis needs at least 2",
        )
    if control in levels:
        raise reader.fail(
            "design",
            "levels",
            f"contains the control level {control!r}; list only the exposure levels, "
            f"ordered low to high",
        )

    dose_raw = reader.items("design", "dose", default="")
    dose: tuple[float, ...] | None = None
    if dose_raw:
        if len(dose_raw) != len(levels):
            raise reader.fail(
                "design",
                "dose",
                f"has {len(dose_raw)} value(s) but [design] levels has {len(levels)}; "
                f"the two lists are index-matched, so a mismatch would shift every dose "
                f"by one position",
            )
        values = []
        for entry in dose_raw:
            try:
                values.append(float(entry))
            except ValueError:
                raise reader.fail(
                    "design", "dose", f"{entry!r} is not a number"
                ) from None
        # `levels` is declared low to high, so an index-matched dose vector has to rise
        # with it. When it does not, one of the two lists is in the wrong order - and a
        # reversed `levels` is otherwise undetectable: it produces a full, plausible run
        # in which every linear slope has the wrong sign and the retention rule inverts.
        # This is the only automatic guard against that, and it is why supplying `dose` is
        # worth doing even when the fit is on `rank`.
        falling = [
            f"{levels[index]}={values[index]:g} then {levels[index + 1]}={values[index + 1]:g}"
            for index in range(len(values) - 1)
            if values[index + 1] <= values[index]
        ]
        if falling:
            raise reader.fail(
                "design",
                "dose",
                f"does not rise with [design] levels, which are declared low to high: "
                f"{'; '.join(falling[:3])}\n"
                f"  either the doses or the levels are in the wrong order; a reversed "
                f"levels list runs to completion and silently inverts the retention rule",
            )
        dose = tuple(values)

    excluded = frozenset(reader.items("design", "exclude_from_fit", default=""))
    unknown = sorted(excluded - set(levels))
    if unknown:
        raise reader.fail(
            "design",
            "exclude_from_fit",
            f"names level(s) absent from [design] levels: {', '.join(unknown)}",
        )
    if len(excluded) == len(levels):
        raise reader.fail(
            "design", "exclude_from_fit", "withholds every level, leaving nothing to fit"
        )
    return DesignSpec(
        control=control, levels=levels, dose=dose, exclude_from_fit=excluded
    )


def _read_trim(reader: _Reader) -> TrimSpec:
    enabled = reader.boolean("trim", "enabled")
    # When trimming is on the percentiles are required: enabling the step and then
    # falling back to 0/100 would silently trim nothing.
    bounds_default = None if enabled else "0"
    lower = reader.number(
        "trim", "lower_percentile", default=bounds_default, minimum=0.0, maximum=100.0
    )
    upper = reader.number(
        "trim",
        "upper_percentile",
        default=None if enabled else "100",
        minimum=0.0,
        maximum=100.0,
    )
    if lower >= upper:
        raise reader.fail(
            "trim",
            "lower_percentile",
            f"is {lower:g} and upper_percentile is {upper:g}; the kept interval "
            f"[lower, upper] would be empty",
        )
    scope = reader.unique("trim", "scope", reader.items("trim", "scope", default=""))
    if enabled and not scope:
        raise reader.fail(
            "trim",
            "scope",
            "is required when [trim] enabled is true: percentiles are computed within "
            "each group of these columns",
        )
    return TrimSpec(
        enabled=enabled, lower_percentile=lower, upper_percentile=upper, scope=scope
    )


def _read_emd(
    reader: _Reader, design: DesignSpec, schema: SchemaSpec
) -> EmdSpec:
    raw = reader.items("emd", "contrasts", default=CONTRASTS_KEYWORD)
    # Every level, including any withheld from fitting: EMD is computed for the whole
    # series and the fitter is what skips the withheld levels.
    known = {f"{design.control}v{level}": (design.control, level) for level in design.levels}
    if raw == (CONTRASTS_KEYWORD,):
        keyword: str | None = CONTRASTS_KEYWORD
        contrasts = tuple(known[f"{design.control}v{lv}"] for lv in design.levels)
    elif CONTRASTS_KEYWORD in raw:
        raise reader.fail(
            "emd",
            "contrasts",
            f"mixes the keyword {CONTRASTS_KEYWORD!r} with explicit contrasts; use one "
            f"or the other",
        )
    else:
        keyword = None
        unknown = [entry for entry in raw if entry not in known]
        if unknown:
            raise reader.fail(
                "emd",
                "contrasts",
                f"has entries that are neither the keyword {CONTRASTS_KEYWORD!r} nor a "
                f"control-vs-level pair: {', '.join(unknown)}\n"
                f"  valid pairs: {', '.join(known)}",
            )
        order = {level: index for index, level in enumerate(design.levels)}
        contrasts = tuple(
            sorted((known[entry] for entry in dict.fromkeys(raw)), key=lambda p: order[p[1]])
        )

    baseline = reader.choice("emd", "baseline", BASELINES, default=BASELINES[0])
    per_replicate = reader.boolean("emd", "per_replicate", default="false")
    if baseline == "control_across_replicates" and not schema.pool_over:
        raise reader.fail(
            "emd",
            "baseline",
            "is 'control_across_replicates', which pairs replicates against each other, "
            "but [schema] pool_over names no replicate column",
        )
    if per_replicate and not schema.pool_over:
        raise reader.fail(
            "emd",
            "per_replicate",
            "is true, which computes one EMD per replicate, but [schema] pool_over names "
            "no replicate column",
        )
    return EmdSpec(
        contrasts=contrasts,
        contrasts_keyword=keyword,
        baseline=baseline,
        per_replicate=per_replicate,
    )


def _read_fit(reader: _Reader, design: DesignSpec, *, selecting: bool) -> FitSpec:
    models = reader.unique(
        "fit", "models", reader.items("fit", "models", default=",".join(MODELS))
    )
    unknown = [model for model in models if model not in MODELS]
    if unknown:
        raise reader.fail(
            "fit",
            "models",
            f"names unknown model(s): {', '.join(unknown)}\n"
            f"  available: {', '.join(MODELS)}",
        )
    x_scale = reader.choice("fit", "x_scale", X_SCALES, default="rank")
    if x_scale == "dose" and design.dose is None:
        raise reader.fail(
            "fit",
            "x_scale",
            "is 'dose', which needs the actual doses, but [design] dose is empty",
        )
    if x_scale == "dose" and design.dose is not None:
        # The four sigmoid models evaluate log(x) internally, so a non-positive dose puts
        # a negative number into a logarithm and the fit silently becomes NaN.
        nonpositive = [
            f"{level}={value:g}"
            for level, value in zip(design.levels, design.dose)
            if value <= 0
        ]
        if nonpositive:
            raise reader.fail(
                "fit",
                "x_scale",
                f"is 'dose' but [design] dose has non-positive value(s): "
                f"{', '.join(nonpositive)}\n"
                f"  the concentration-response models take log(x), so a dose of zero or "
                f"less cannot be fitted",
            )

    spec = FitSpec(
        models=models, x_scale=x_scale, rank_by=reader.choice(
            "fit", "rank_by", RANK_BY, default="aic_plus_bic"
        )
    )
    # Only when something will actually be fitted. With [select] enabled false the models
    # are never called, and refusing a four-level design for a model nobody runs would be
    # demanding a fix to an irrelevant key.
    fitted = len(design.fitted_levels)
    if selecting and fitted <= spec.max_params:
        worst = max(spec.models, key=lambda model: MODEL_PARAMS[model])
        raise reader.fail(
            "fit",
            "models",
            f"includes {worst}, which has {MODEL_PARAMS[worst]} free parameters, but only "
            f"{fitted} level(s) are fitted after [design] exclude_from_fit; the fit needs "
            f"more points than parameters",
        )
    return spec


def _read_select(
    reader: _Reader, schema: SchemaSpec, fit: FitSpec, *, enabled: bool
) -> SelectSpec:
    slope_positive = reader.choice(
        "select", "slope_positive", QUANTIFIERS, default="any"
    )
    nonconstant = reader.choice("select", "nonconstant", QUANTIFIERS, default="all")
    if enabled and "Lin" not in fit.models:
        raise reader.fail(
            "select",
            "slope_positive",
            f"is {slope_positive!r}, which tests the sign of the linear slope, but "
            f"[fit] models does not include Lin",
        )
    if enabled and "Con" not in fit.models:
        raise reader.fail(
            "select",
            "nonconstant",
            f"is {nonconstant!r}, which tests whether the constant model is the best "
            f"fit, but [fit] models does not include Con",
        )
    strata = reader.items("select", "strata", default="")
    if strata and not schema.group_by:
        raise reader.fail(
            "select",
            "strata",
            f"names strata ({', '.join(strata)}) but [schema] group_by is empty, so the "
            f"experiment has a single stratum",
        )
    return SelectSpec(
        enabled=enabled,
        slope_positive=slope_positive,
        nonconstant=nonconstant,
        strata=reader.unique("select", "strata", strata) if strata else None,
    )


def _read_drop_missing(reader: _Reader) -> DropMissingSpec:
    # A filter that changes which columns leave the pipeline has to be stated, for the same
    # reason [trim] enabled and [correlation] enabled are required rather than defaulted.
    return DropMissingSpec(
        enabled=reader.boolean("drop_missing", "enabled"),
        # Exclusive at zero: "missing 0% or more" is true of every feature, so a threshold of
        # zero would empty the table rather than filter it.
        max_missing=reader.number(
            "drop_missing",
            "max_missing",
            default="30",
            minimum=0.0,
            maximum=100.0,
            exclude_minimum=True,
        ),
        # Exact names, not patterns: an exclusion is a judgement about one feature, and a
        # regex that quietly matches a second one is the wrong tool for that.
        exclude=reader.unique(
            "drop_missing", "exclude", reader.items("drop_missing", "exclude", default="")
        ),
    )


def _read_correlation(reader: _Reader, trim: TrimSpec) -> CorrelationSpec:
    enabled = reader.boolean("correlation", "enabled")
    # Absent means "the same unit trimming works within", which is the published run's well.
    # Present but empty is a deliberate choice to correlate the rows as they are, so the two
    # cases have to be told apart rather than both read as "no columns".
    if reader.parser.has_option("correlation", "aggregate_by"):
        aggregate_by = reader.unique(
            "correlation", "aggregate_by", reader.items("correlation", "aggregate_by")
        )
    else:
        aggregate_by = trim.scope
        reader.defaulted.append("correlation.aggregate_by")
        if enabled and not aggregate_by:
            raise reader.fail(
                "correlation",
                "aggregate_by",
                "is required when [correlation] enabled is true and [trim] scope is empty: "
                "correlations are computed between features aggregated within each group of "
                "these columns. Set it to the columns identifying one experimental unit, or "
                "set it to nothing to correlate the rows as they are",
            )
    return CorrelationSpec(
        enabled=enabled,
        threshold=reader.number(
            "correlation",
            "threshold",
            default="0.9",
            minimum=0.0,
            maximum=1.0,
            exclude_minimum=True,
        ),
        linkage=reader.choice("correlation", "linkage", LINKAGES, default="average"),
        representative=reader.choice(
            "correlation", "representative", REPRESENTATIVES, default="alphabetical"
        ),
        aggregate_by=aggregate_by,
        fill_missing=reader.choice(
            "correlation", "fill_missing", FILL_MISSING, default="column_mean"
        ),
    )
