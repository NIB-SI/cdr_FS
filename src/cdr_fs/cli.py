"""Command line entry point.

One subcommand per stage, plus `check`, which validates a configuration and reports the
metadata/feature split without running anything, and `run`, which executes the chain.
Every subcommand takes `-c/--config`, which is what removes both the `/PATH/TO/...`
editing ritual and the run-from-the-script's-own-directory constraint of the original
scripts.

Each stage is two functions: a handler that parses `argparse` and loads the
configuration, and a `_stage_*` that takes an already-loaded `Config` and explicit
options. `run` loads the configuration once and calls the second of each pair, so a
misconfigured run fails before the first stage rather than between two of them.

Stage modules are imported inside their handlers, so `--help` and `check` work on an
install without pandas, scipy or matplotlib present.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from cdr_fs import __version__
from cdr_fs.config import Config, ConfigError, load_config

__all__ = ["main"]

_EXTENSIONS = {"tab": ".tsv", "comma": ".csv", "semicolon": ".csv"}

#: `cdr-fs plot --only` names, and what each figure needs to exist.
FIGURES = ("fits", "emd", "baseline", "dendrogram")

#: The stages `cdr-fs run` executes, in order. `trim` is deliberately not one of them:
#: every stage trims from `[input] table` itself, so a trimmed copy is an inspection
#: artefact rather than an input - 3.9 GB of one on the reference dataset.
CHAIN = ("emd", "fit", "select", "correlation", "drop_missing", "plot")

#: The three that are the concentration-response selection. `[select] enabled` turns all
#: three off together, because `select` needs `fit` and `fit` needs `emd`, and with nothing
#: reading their tables the two expensive stages are work for its own sake.
SELECTION = ("emd", "fit", "select")

#: Width of the rules `run` prints between stages.
_RULE_WIDTH = 78


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help()
        return 0
    try:
        return args.handler(args)
    except ConfigError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:  # pragma: no cover
        print("interrupted", file=sys.stderr)
        return 130


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cdr-fs",
        # Raw, so the two epilog lines stay two lines: they point at different files for
        # different jobs, and reflowed into one paragraph that distinction is lost.
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "Feature selection for high-content screening by fitting\n"
            "concentration/dose-response models to earth mover's distance scores."
        ),
        epilog=(
            "Start from examples/template.ini for a new dataset.\n"
            "examples/published.ini is the published run:\n  "
            "https://doi.org/10.1021/acs.est.5c18316"
        ),
    )
    parser.add_argument("--version", action="version", version=f"cdr-fs {__version__}")
    subparsers = parser.add_subparsers(dest="command", metavar="<stage>")

    check = _add(subparsers, "check", "validate a configuration and report the schema")
    check.add_argument(
        "--scan",
        action="store_true",
        help=(
            "also read the condition and group_by columns to confirm the declared levels "
            "and strata occur in the data; one pass over the whole table"
        ),
    )
    check.set_defaults(handler=_check)

    run = _add(
        subparsers,
        "run",
        "run the whole chain: " + ", ".join(CHAIN),
        description=(
            "Run the whole chain in order: "
            + " -> ".join(CHAIN)
            + ". The configuration is read once, at the start, and each stage prints its "
            "own report. A stage a switch turns off is reported as skipped rather than as "
            "a failure: [select] enabled false skips "
            + ", ".join(SELECTION)
            + " together, carrying every feature into the filtering stages, and "
            "[correlation] enabled false skips that one. `drop_missing` always runs, "
            "because it writes the final table. Two things are separate commands on "
            "purpose and no run calls them: `trim`, which only writes a trimmed copy of "
            "the input for inspection, and `plot` over every feature - the figures drawn "
            "here cover the retained features only, and a figure that will not draw is "
            "reported without failing the run. Exit codes are 0 when every planned stage "
            "finished, 2 for a configuration error, and 3 when a stage refused and the "
            "chain stopped there."
        ),
    )
    run.set_defaults(handler=_run)

    trim = _add(subparsers, "trim", "remove extreme feature values (optional QC step)")
    trim.add_argument(
        "--dry-run",
        action="store_true",
        help="report what would be trimmed without writing the trimmed table",
    )
    trim.set_defaults(handler=_trim)

    emd = _add(
        subparsers,
        "emd",
        "earth mover's distance between control and each exposure level",
    )
    emd.set_defaults(handler=_emd)

    fit = _add(
        subparsers, "fit", "fit the concentration-response models to the distance series"
    )
    fit.set_defaults(handler=_fit)

    select = _add(subparsers, "select", "apply the retention rule to the fitted models")
    select.set_defaults(handler=_select)

    correlation = _add(
        subparsers,
        "correlation",
        "collapse near-redundant features by correlation (optional)",
    )
    correlation.set_defaults(handler=_correlation)

    # The help opens on the table because the stage name does not: `drop_missing` names the
    # filter this stage applies, but its product is the final data table, and a reader
    # scanning `--help` for "where does my output come from" has to land here.
    drop_missing = _add(
        subparsers,
        "drop_missing",
        "write the final table - the input restricted to the retained features, minus "
        "those missing too much data",
    )
    drop_missing.add_argument(
        "--features",
        metavar="FILE",
        help=(
            "feature list to apply, one name per line; the default is the representatives "
            "list when [correlation] is enabled, the selected list when only [select] is, "
            "and every feature when neither is"
        ),
    )
    drop_missing.set_defaults(handler=_drop_missing)

    plot = _add(subparsers, "plot", "draw the diagnostic figures from the tables written so far")
    plot.add_argument(
        "--only",
        metavar="FIGURE",
        help=(
            "comma-separated subset of: "
            + ", ".join(FIGURES)
            + "; the default draws whichever the available tables allow"
        ),
    )
    plot.add_argument(
        "--grid",
        type=int,
        default=3,
        metavar="N",
        help="fit panels per row and column of a page (default 3, so 9 per page)",
    )
    plot.add_argument(
        "--features",
        metavar="FILE",
        help=(
            "restrict the fit panels to this feature list; the default draws every feature "
            "in the distance table, which on a full run is hundreds of pages"
        ),
    )
    plot.set_defaults(handler=_plot)

    return parser


def _add(
    subparsers, name: str, help_text: str, description: str | None = None
) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        name, help=help_text, description=description or help_text
    )
    parser.add_argument(
        "-c",
        "--config",
        required=True,
        metavar="FILE",
        help="configuration file (.ini)",
    )
    return parser


# ------------------------------------------------------------------------ reporting


def _human_bytes(count: int) -> str:
    size = float(count)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if size < 1024 or unit == "TiB":
            return f"{size:.0f} {unit}" if unit == "B" else f"{size:.2f} {unit}"
        size /= 1024
    return f"{size:.2f} TiB"  # pragma: no cover


def _field(label: str, value: str) -> str:
    return f"  {label:<18}{value}"


def _abbreviate(values, limit: int = 12) -> str:
    values = list(values)
    if len(values) <= limit:
        return ", ".join(str(value) for value in values)
    head = ", ".join(str(value) for value in values[: limit - 1])
    return f"{head}, ... {values[-1]}  ({len(values)} in total)"


def _exposure_axis(config: Config) -> str:
    """Which end of `[design] levels` is the low exposure, spelled out.

    Not left to "(low to high)": the order of that list is the response axis, and a list
    written the wrong way round runs to completion, negates every slope and retains
    nothing. This line is where a reader can catch that, which is why both `check` and
    `run` print it.
    """
    design = config.design
    ends = f"{design.levels[0]} is the LOWEST exposure, {design.levels[-1]} the highest"
    if design.dose is not None:
        ends += f"   ({design.dose[0]:g} -> {design.dose[-1]:g})"
    return ends


def _describe(config: Config) -> list[str]:
    design, emd, fit = config.design, config.emd, config.fit
    lines = [
        "[input]",
        _field("table", f"{config.input.table}  ({_human_bytes(config.input.table.stat().st_size)})"),
        _field("sep", config.input.sep_keyword),
        "[schema]",
        _field("condition", config.schema.condition),
        _field("group_by", config.schema.group_by or "(none - a single stratum)"),
        _field("pool_over", config.schema.pool_over or "(none)"),
        _field("metadata_patterns", _abbreviate(config.schema.metadata_patterns)),
        "[design]",
        _field("control", design.control),
        _field("levels", _abbreviate(design.levels)),
    ]
    lines.append(_field("exposure axis", _exposure_axis(config)))
    if design.dose is not None:
        lines.append(_field("dose", _abbreviate(f"{value:g}" for value in design.dose)))
    else:
        lines.append(
            _field("dose", "(none - levels are used by rank, and nothing checks their order)")
        )
    lines += [
        _field(
            "exclude_from_fit",
            _abbreviate(sorted(design.exclude_from_fit)) if design.exclude_from_fit else "(none)",
        ),
        _field("fitted levels", f"{len(design.fitted_levels)}: {_abbreviate(design.fitted_levels)}"),
        "[trim]",
        _field("enabled", str(config.trim.enabled).lower()),
    ]
    if config.trim.enabled:
        lines += [
            _field(
                "interval kept",
                f"[p{config.trim.lower_percentile:g}, p{config.trim.upper_percentile:g}]",
            ),
            _field("scope", ", ".join(config.trim.scope)),
        ]
    lines += [
        "[emd]",
        _field(
            "contrasts",
            f"{len(emd.contrasts)}: "
            f"{_abbreviate(f'{ref}v{level}' for ref, level in emd.contrasts)}",
        ),
        _field("baseline", emd.baseline),
        _field("per_replicate", str(emd.per_replicate).lower()),
        "[fit]",
        _field("models", ", ".join(fit.models)),
        _field("x_scale", fit.x_scale),
        _field("rank_by", fit.rank_by),
        "[select]",
        _field("enabled", str(config.select.enabled).lower()),
    ]
    # As for [trim] and [correlation]: the rule is only described when it is applied, so a
    # reader cannot mistake a printed quantifier for a gate that ran.
    if config.select.enabled:
        lines += [
            _field(
                "slope_positive",
                f"{config.select.slope_positive}  (a positive linear slope on "
                f"{'at least one' if config.select.slope_positive == 'any' else 'every'} "
                f"stratum)",
            ),
            _field(
                "nonconstant",
                f"{config.select.nonconstant}  (the constant model is not the best fit on "
                f"{'at least one' if config.select.nonconstant == 'any' else 'every'} "
                f"stratum)",
            ),
            _field(
                "strata",
                ", ".join(config.select.strata)
                if config.select.strata
                else "(every stratum present in the data)",
            ),
        ]
    lines += [
        "[correlation]",
        _field("enabled", str(config.correlation.enabled).lower()),
    ]
    if config.correlation.enabled:
        lines += [
            _field("threshold", f"|r| >= {config.correlation.threshold:g}"),
            _field("linkage", config.correlation.linkage),
            _field("representative", config.correlation.representative),
            _field(
                "aggregate_by",
                ", ".join(config.correlation.aggregate_by)
                if config.correlation.aggregate_by
                else "(none - the rows are correlated as they are)",
            ),
            _field("fill_missing", config.correlation.fill_missing),
        ]
    lines += [
        "[drop_missing]",
        _field("enabled", str(config.drop_missing.enabled).lower()),
    ]
    if config.drop_missing.enabled:
        lines.append(
            _field(
                "max_missing",
                f"{config.drop_missing.max_missing:g}%  (a feature missing this much of the "
                f"table or more is dropped)",
            )
        )
    lines.append(
        _field(
            "exclude",
            f"{len(config.drop_missing.exclude)}: "
            f"{_abbreviate(config.drop_missing.exclude, limit=4)}"
            if config.drop_missing.exclude
            else "(none)",
        )
    )
    lines += ["[output]", _field("dir", str(config.output.dir))]
    return lines


def _describe_schema(config: Config, columns: list[str]) -> list[str]:
    from cdr_fs.schema import resolve_schema

    resolved = resolve_schema(columns, config.schema.compiled)
    breakdown = resolved.prefix_breakdown()
    lines = [
        f"[columns]  {len(columns)} = {len(resolved.metadata)} metadata "
        f"+ {len(resolved.features)} feature(s)",
        _field("metadata", _abbreviate(resolved.metadata, limit=12)),
        _field(
            "features",
            ", ".join(f"{head}* {count}" for head, count in breakdown),
        ),
    ]
    # Name the members of every small group: the two organelle counts among 461
    # rp_norm_* features are exactly the case where a lazy metadata pattern goes wrong.
    for head, count in breakdown:
        if count <= 5:
            lines.append(
                _field(f"  {head}*", ", ".join(resolved.features_with_prefix(head)))
            )
    return lines


# -------------------------------------------------------------------------- handlers


def _check(args: argparse.Namespace) -> int:
    return _stage_check(load_config(args.config), scan=args.scan)


def _stage_check(config: Config, *, scan: bool = False) -> int:
    from cdr_fs.schema import read_header

    print(f"cdr_FS {__version__} - configuration check")
    print(_field("config", str(config.path)))
    for line in _describe(config):
        print(line)

    columns = read_header(config.input.table, config.input.sep)
    warnings = config.validate_columns(columns)
    print()
    for line in _describe_schema(config, columns):
        print(line)

    if scan:
        warnings += _scan(config)
    else:
        print()
        print(
            "note: the declared levels and strata were not checked against the data; "
            "pass --scan to read those two columns"
        )

    if config.defaulted:
        print()
        print("defaults in effect (absent from the file):")
        for key in config.defaulted:
            print(f"  {key}")
    if warnings:
        print()
        for warning in warnings:
            print(f"warning: {warning}")
    print()
    print("configuration is valid")
    return 0


def _scan(config: Config) -> list[str]:
    import pandas as pd

    usecols = [config.schema.condition]
    for column in (config.schema.group_by, config.schema.pool_over):
        if column and column not in usecols:
            usecols.append(column)
    frame = pd.read_csv(
        config.input.table,
        sep=config.input.sep,
        usecols=usecols,
        dtype={name: str for name in usecols},
    )
    levels = set(frame[config.schema.condition].str.strip())
    strata = (
        set(frame[config.schema.group_by].str.strip()) if config.schema.group_by else None
    )
    replicates = (
        set(frame[config.schema.pool_over].str.strip()) if config.schema.pool_over else None
    )
    print()
    print(f"[data]  {len(frame):,} row(s)")
    print(_field(config.schema.condition, _abbreviate(sorted(levels), limit=12)))
    if strata is not None:
        print(_field(config.schema.group_by, _abbreviate(sorted(strata), limit=12)))
    if replicates is not None:
        print(_field(config.schema.pool_over, _abbreviate(sorted(replicates), limit=12)))
    return config.validate_observed(levels=levels, strata=strata, replicates=replicates)


def _prepare(config: Config, *, trim: bool = True):
    """Read the input table, validate it against the configuration, and trim it.

    Every stage that needs cell-level data starts here, so trimming is applied once from
    `[input] table` rather than being inherited from a materialised intermediate. That is
    why `cdr-fs trim` is optional: it exists to write the trimmed table out for inspection,
    not because later stages need the file.
    """
    from cdr_fs.schema import read_header, read_table, resolve_schema
    from cdr_fs.trim import trim_extremes

    columns = read_header(config.input.table, config.input.sep)
    for warning in config.validate_columns(columns):
        print(f"warning: {warning}", file=sys.stderr)
    resolved = resolve_schema(columns, config.schema.compiled)

    print(f"reading {config.input.table} ...", flush=True)
    frame = read_table(
        config.input.table,
        config.input.sep,
        metadata=resolved.metadata,
        features=resolved.features,
    )
    print(f"  {len(frame):,} row(s), {len(resolved.features)} feature(s)")
    for warning in config.validate_observed(
        levels=set(frame[config.schema.condition]),
        strata=set(frame[config.schema.group_by]) if config.schema.group_by else None,
        replicates=set(frame[config.schema.pool_over]) if config.schema.pool_over else None,
    ):
        print(f"warning: {warning}", file=sys.stderr)

    if trim and config.trim.enabled:
        frame, report = trim_extremes(
            frame,
            resolved.features,
            config.trim.scope,
            config.trim.lower_percentile,
            config.trim.upper_percentile,
            inplace=True,
        )
        print(report.summary())
    return frame, resolved


def _table_path(config: Config, name: str) -> Path:
    """Where a stage table lives. One convention, so the stages can find each other."""
    return Path(config.output.dir) / f"{name}{_EXTENSIONS[config.input.sep_keyword]}"


def _write(config: Config, table, name: str):
    destination = Path(config.output.dir)
    destination.mkdir(parents=True, exist_ok=True)
    target = _table_path(config, name)
    table.to_csv(target, sep=config.input.sep, index=False)
    print(f"wrote {target}  ({_human_bytes(target.stat().st_size)})")
    return target


def _list_path(config: Config, name: str) -> Path:
    """Where a feature list lives. One convention, so the stages can find each other."""
    return Path(config.output.dir) / f"{name}.txt"


def _collapsed_prefix(config: Config) -> str:
    """What `correlation`'s outputs are called: nothing extra, or `all_`.

    A correlation of every feature is not a correlation of the selected ones, and the two
    must not overwrite each other's answers in a shared output directory. The method's own
    chain keeps the unqualified names; a run with no selection qualifies its outputs with
    the input they were built from.
    """
    return "" if config.select.enabled else "all_"


def _features_stem(config: Config) -> str | None:
    """Which list the filtering stages apply, by name, or None for every feature.

    `correlation` narrows to its representatives and `select` to its retained list, so the
    later stage wins. With both switched off nothing has narrowed anything, and the whole
    feature set is what there is to filter.
    """
    if config.correlation.enabled:
        return f"{_collapsed_prefix(config)}representatives"
    if config.select.enabled:
        return "selected"
    return None


def _features_source(config: Config, override: str | None) -> Path | None:
    """The feature list `drop_missing` applies, or None when that is every feature."""
    if override:
        source = Path(override)
        if not source.is_file():
            raise ConfigError(f"--features file not found: {source}")
        return source
    stem = _features_stem(config)
    if stem is None:
        return None
    source = _list_path(config, stem)
    if not source.exists():
        producer = "correlation" if config.correlation.enabled else "select"
        raise ConfigError(
            f"{source} is missing - run `cdr-fs {producer} -c {config.path}` first, "
            f"or name a list with --features"
        )
    return source


def _write_list(config: Config, features: list[str], name: str) -> Path:
    destination = Path(config.output.dir)
    destination.mkdir(parents=True, exist_ok=True)
    target = _list_path(config, name)
    target.write_text("\n".join(features) + "\n", encoding="utf-8")
    print(f"wrote {target}  ({len(features)} feature(s))")
    return target


def _read_stage(config: Config, name: str, produced_by: str):
    """Read a previous stage's output, or say which command produces it."""
    from cdr_fs.schema import read_stage_table

    target = _table_path(config, name)
    if not target.exists():
        raise ConfigError(
            f"{target} is missing - run `cdr-fs {produced_by} -c {config.path}` first"
        )
    return read_stage_table(target, config.input.sep)


def _trim(args: argparse.Namespace) -> int:
    return _stage_trim(load_config(args.config), dry_run=args.dry_run)


def _stage_trim(config: Config, *, dry_run: bool = False) -> int:
    if not config.trim.enabled:
        print(
            "error: [trim] enabled is false, so there is nothing to do",
            file=sys.stderr,
        )
        return 3

    frame, _ = _prepare(config)
    if dry_run:
        print("dry run - nothing written")
        return 0
    _write(config, frame, "trimmed")
    return 0


def _emd(args: argparse.Namespace) -> int:
    return _stage_emd(load_config(args.config))


def _stage_emd(config: Config) -> int:
    from cdr_fs.emd import compute_baseline, compute_contrasts

    frame, resolved = _prepare(config)

    contrasts, report = compute_contrasts(config, frame, resolved.features)
    print(report.summary())
    _write(config, contrasts, "emd")

    if config.emd.baseline != "none":
        baseline, baseline_report = compute_baseline(config, frame, resolved.features)
        print(baseline_report.summary())
        _write(config, baseline, "emd_baseline")
    return 0


def _fit(args: argparse.Namespace) -> int:
    return _stage_fit(load_config(args.config))


def _stage_fit(config: Config) -> int:
    from cdr_fs.fit import fit_series

    # `load_config` raises this only when [select] enabled is true. This stage runs either
    # way, so it asks for itself: an over-parameterised `curve_fit` raises `TypeError`,
    # which reaches the user as a traceback and an exit code of 1.
    problem = config.fit_problem
    if problem:
        print(f"error: {config.path}: [fit] models - {problem}", file=sys.stderr)
        return 2

    distances = _read_stage(config, "emd", "emd")
    table, report = fit_series(config, distances)
    print(report.summary())
    if table.empty:
        # An empty fit table is not a result: `select` would read it and retain nothing,
        # which looks like "no feature responds" rather than "nothing was fitted".
        print(
            "error: no series was fitted, so there is nothing to select from - see the "
            "report above for which features were incomplete",
            file=sys.stderr,
        )
        return 3
    _write(config, table, "fit")
    return 0


def _select(args: argparse.Namespace) -> int:
    return _stage_select(load_config(args.config))


def _stage_select(config: Config) -> int:
    from cdr_fs.select import select_features

    # Same shape as `trim` and `correlation`: `run` reads the switch itself and never calls
    # the stage when it is off, so this is only ever reached by `cdr-fs select`.
    if not config.select.enabled:
        print(
            "error: [select] enabled is false, so there is nothing to do",
            file=sys.stderr,
        )
        return 3

    fits = _read_stage(config, "fit", "fit")
    retained, evidence, report = select_features(config, fits)
    print(report.summary())

    _write(config, evidence, "select_evidence")
    _write_list(config, retained, "selected")
    return 0


def _correlation(args: argparse.Namespace) -> int:
    return _stage_correlation(load_config(args.config))


def _stage_correlation(config: Config) -> int:
    from cdr_fs.correlation import collapse_correlated
    from cdr_fs.drop_missing import read_feature_list

    # `run` reads this switch itself and does not call the stage when it is off, so this
    # refusal is only ever reached by `cdr-fs correlation`, where the user asked for the
    # stage by name and an exit code is the honest answer.
    if not config.correlation.enabled:
        print(
            "error: [correlation] enabled is false, so there is nothing to do",
            file=sys.stderr,
        )
        return 3

    if config.select.enabled:
        selected = _list_path(config, "selected")
        if not selected.exists():
            raise ConfigError(
                f"{selected} is missing - run `cdr-fs select -c {config.path}` first"
            )
        features = read_feature_list(selected)
        if not features:
            print(
                f"error: {selected.name} is empty, so there is nothing to collapse",
                file=sys.stderr,
            )
            return 3
        frame, _ = _prepare(config)
    else:
        # No selection ran, so there is no retained list to collapse: the redundancy
        # question is asked of the whole feature set instead.
        print("collapsing every feature ([select] enabled is false)")
        frame, resolved = _prepare(config)
        features = resolved.features

    kept, clusters, tree, report = collapse_correlated(config, frame, features)
    print(report.summary())
    prefix = _collapsed_prefix(config)
    _write(config, clusters, f"{prefix}correlation_clusters")
    _write(config, tree, f"{prefix}correlation_linkage")
    _write_list(config, kept, f"{prefix}representatives")
    return 0


def _drop_missing(args: argparse.Namespace) -> int:
    return _stage_drop_missing(load_config(args.config), features=args.features)


def _stage_drop_missing(config: Config, *, features: str | None = None) -> int:
    from cdr_fs.drop_missing import read_feature_list, restrict_table

    source = _features_source(config, features)
    if source is None:
        # Neither stage that narrows the feature set ran, so there is no list to apply and
        # the whole set is the answer. `_prepare` is where it comes from.
        print("applying every feature (no selection, no correlation collapsing)")
        frame, resolved = _prepare(config)
        names, stem = resolved.features, "all"
    else:
        print(f"applying {source.name}")
        names = read_feature_list(source)
        if not names:
            # Otherwise the output is a metadata-only table, which looks like a result.
            print(
                f"error: {source.name} is empty, so there are no features to write",
                file=sys.stderr,
            )
            return 3
        frame, resolved = _prepare(config)
        stem = source.stem
    final, quality, report = restrict_table(
        frame,
        resolved.metadata,
        names,
        drop_missing=config.drop_missing.enabled,
        max_missing=config.drop_missing.max_missing,
        exclude=config.drop_missing.exclude,
    )
    print(report.summary())

    # Named for the list applied, not just for the stage: two lists over one table are a
    # normal thing to run, and they must not overwrite each other's output.
    name = f"final_{stem}"
    _write(config, final, name)
    _write(config, quality, f"{name}_features")
    _write_list(config, list(quality.loc[~quality["dropped"], "feature"]), f"{name}_retained")
    return 0


def _plot(args: argparse.Namespace) -> int:
    return _stage_plot(
        load_config(args.config), only=args.only, grid=args.grid, features=args.features
    )


def _stage_plot(
    config: Config,
    *,
    only: str | None = None,
    grid: int = 3,
    features: str | None = None,
) -> int:
    import pandas as pd

    from cdr_fs.drop_missing import read_feature_list
    from cdr_fs.plots import plot_dendrogram, plot_distribution, plot_fit_panels
    from cdr_fs.schema import read_stage_table

    wanted = FIGURES
    if only:
        wanted = tuple(part.strip() for part in only.split(",") if part.strip())
        unknown = [name for name in wanted if name not in FIGURES]
        if unknown:
            raise ConfigError(
                f"--only names unknown figure(s): {', '.join(unknown)}\n"
                f"  available: {', '.join(FIGURES)}"
            )

    results = Path(config.output.dir)
    extension = _EXTENSIONS[config.input.sep_keyword]

    def available(name: str):
        target = results / f"{name}{extension}"
        return target if target.exists() else None

    selection = None
    if features:
        source = Path(features)
        if not source.is_file():
            raise ConfigError(f"--features file not found: {source}")
        selection = read_feature_list(source)

    drawn: list[Path] = []
    skipped: list[str] = []
    read = lambda target: read_stage_table(target, config.input.sep)  # noqa: E731

    # A figure function refuses rather than draw something misleading - an empty
    # distribution, a truncated tree. That refusal is a ValueError, and it must reach the
    # user as a skipped figure with a reason, not as a traceback: `plot` draws several
    # figures and one bad table should not lose the others.
    def draw(name: str, call):
        try:
            drawn.extend(call())
        except ValueError as refusal:
            skipped.append(f"{name} - {refusal}")

    if "fits" in wanted:
        distances, fits = available("emd"), available("fit")
        if distances and fits:
            draw("fits", lambda: plot_fit_panels(
                config,
                read(distances),
                read(fits),
                results,
                grid=grid,
                features=selection,
            ))
        else:
            skipped.append("fits - needs emd and fit (run `cdr-fs fit`)")

    for name, source, title in (
        ("emd", "emd", "Distances from the control, per feature"),
        ("baseline", "emd_baseline", "Between-replicate control distances, per feature"),
    ):
        if name not in wanted:
            continue
        target = available(source)
        if target:
            draw(name, lambda target=target, source=source, title=title: [
                plot_distribution(read(target), results / f"{source}.png", title=title)
            ])
        elif name == "baseline" and config.emd.baseline == "none":
            # `emd` ran and will never write this file, so naming it as the fix is wrong.
            skipped.append(
                "baseline - [emd] baseline is none, so no between-replicate distances "
                "were computed"
            )
        else:
            skipped.append(f"{name} - needs {source}{extension} (run `cdr-fs emd`)")

    if "dendrogram" in wanted:
        # Only when the stage is on: a tree left in the output directory by an earlier
        # run would otherwise be drawn as though it belonged to this one.
        prefix = _collapsed_prefix(config)
        tree = (
            available(f"{prefix}correlation_linkage")
            if config.correlation.enabled
            else None
        )
        clusters = available(f"{prefix}correlation_clusters")
        if tree:
            draw("dendrogram", lambda: [
                plot_dendrogram(
                    pd.read_csv(tree, sep=config.input.sep),
                    results / f"{prefix}dendrogram.png",
                    cut=1.0 - config.correlation.threshold,
                    clusters=read(clusters) if clusters else None,
                )
            ])
        elif not config.correlation.enabled:
            skipped.append(
                "dendrogram - [correlation] enabled is false, so there is no tree to draw"
            )
        else:
            skipped.append(
                f"dendrogram - needs {prefix}correlation_linkage "
                f"(run `cdr-fs correlation`)"
            )

    for path in drawn:
        print(f"wrote {path}  ({_human_bytes(path.stat().st_size)})")
    for note in skipped:
        print(f"skipped {note}", file=sys.stderr)
    if not drawn:
        print("error: nothing to draw", file=sys.stderr)
        return 3
    return 0


# ------------------------------------------------------------------------------- run


def _rule(title: str) -> str:
    """A titled rule, so each stage's own report reads as a block within the run's."""
    opening = f"-- {title} "
    return opening + "-" * max(3, _RULE_WIDTH - len(opening))


def _run(args: argparse.Namespace) -> int:
    return _run_chain(load_config(args.config))


def _run_chain(config: Config) -> int:
    """Run the chain, printing a header, each stage's own output, and a summary.

    The configuration is loaded once, by the caller: a run is fifty seconds on the
    quickstart fixture and a good deal longer on a real dataset, so a misconfiguration
    has to fail before `emd` rather than between two stages.

    Exit codes are the stages' own. A stage the configuration turns off is a declared
    outcome and leaves the run at 0; an exit 3 is a stage refusing to write something
    that would read as a result, and it stops the chain where it happened.
    """
    from cdr_fs.drop_missing import read_feature_list

    selecting, correlated = config.select.enabled, config.correlation.enabled
    # The rule `_stage_drop_missing` applies to pick its input list, so that `plot` and
    # the closing line name the file this run will actually have written.
    stem = _features_stem(config) or "all"
    retained = _list_path(config, f"final_{stem}_retained")
    # Only the figures this run's own tables can support. Left to draw whatever it finds,
    # `plot` would put a previous run's distances in a report that says nothing was fitted.
    drawable = (["fits", "emd", "baseline"] if selecting else []) + (
        ["dendrogram"] if correlated else []
    )
    stages = {
        "emd": lambda: _stage_emd(config),
        "fit": lambda: _stage_fit(config),
        "select": lambda: _stage_select(config),
        "correlation": lambda: _stage_correlation(config),
        "drop_missing": lambda: _stage_drop_missing(config),
        # With the list `drop_missing` has just written. `plot` given no list draws every
        # feature in the distance table, which on a real run is hundreds of pages.
        "plot": lambda: _stage_plot(
            config, only=",".join(drawable), features=str(retained)
        ),
    }

    # The switches are read here rather than left to the stages, which refuse with exit 3
    # when they are off - and a configuration skip is not a failure. Reading them here is
    # also what makes any exit 3 the chain does see a genuine stop.
    off = set()
    if not selecting:
        off |= set(SELECTION)
    if not correlated:
        off.add("correlation")
    if not drawable:
        off.add("plot")
    planned = [stage for stage in CHAIN if stage not in off]
    status = dict.fromkeys(CHAIN, "not reached")
    status.update(dict.fromkeys(off, "off"))

    for line in _run_header(config, planned):
        print(line)

    code = 0
    for number, stage in enumerate(planned, 1):
        print()
        print(_rule(f"{number}/{len(planned)} {stage}"), flush=True)
        result = stages[stage]()
        if result and stage == "plot":
            # The run's product is the final table; the figures are diagnostics, and
            # `plot` has already named on stderr the ones it could not draw. A picture
            # that will not draw must not fail a run that wrote its table - reachable
            # with one feature, where the dendrogram has nothing to cut.
            status[stage] = "no figures"
            continue
        if result:
            status[stage] = "stopped"
            code = result
            break
        status[stage] = "ok"

    print()
    print(_rule("summary"))
    for stage in CHAIN:
        print(_field(stage, status[stage]))
    print(_field("trim", "not run"))  # last: it is not part of the chain
    if status["drop_missing"] == "ok":
        kept = len(read_feature_list(retained))
        print(f"retained {kept} feature(s)")
        print(_field("feature list", str(retained)))
        print(_field("final table", str(_table_path(config, f"final_{stem}"))))
        if not kept:
            # The chain finished and produced nothing. Say where the reason is written -
            # but only where there is something written: no selection, no evidence table.
            where = (
                "select_evidence.tsv carries the linear slope and the winning model for "
                "every feature and stratum"
                if config.select.enabled
                else f"final_{stem}_features.tsv says which rule removed each feature"
            )
            print(f"  nothing survived - {where}")
    return code


def _run_header(config: Config, planned: list[str]) -> list[str]:
    """What this run is about to do, and the two lines of `check` that can stop it.

    Not the whole of `check`: the `[columns]` split and the exposure axis are the two
    things nothing downstream can catch, and they are worth the second they cost. For
    the rest, `cdr-fs check`.
    """
    from cdr_fs.schema import read_header

    lines = [
        f"cdr_FS {__version__} - run",
        _field("config", str(config.path)),
        _field("output", str(config.output.dir)),
        _field("stages", ", ".join(planned)),
    ]
    absent = []
    if "select" not in planned:
        absent.append(
            "emd, fit, select  ([select] enabled is false, so nothing is fitted and "
            "every feature carries forward)"
        )
    if "correlation" not in planned:
        absent.append("correlation  ([correlation] enabled is false)")
    if "plot" not in planned:
        absent.append("plot  (this run writes no table any figure is drawn from)")
    # Two different facts wear the same label. With trimming on, `trim` is absent because
    # its output is an inspection copy no stage reads; with it off, the stage would refuse,
    # and pointing at a command that exits 3 would be worse than saying nothing.
    absent.append(
        "trim  (every stage trims from [input] table itself; `cdr-fs trim` "
        "writes the trimmed copy out for inspection)"
        if config.trim.enabled
        else "trim  ([trim] enabled is false, so nothing is trimmed anywhere)"
    )
    # One entry per line, the rest continuing under the same column.
    lines.append(_field("not run", absent[0]))
    lines += [_field("", note) for note in absent[1:]]
    lines.append(_field("exposure axis", _exposure_axis(config)))
    lines.append("")
    lines += _describe_schema(config, read_header(config.input.table, config.input.sep))
    return lines


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
