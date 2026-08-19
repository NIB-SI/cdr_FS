"""Command line entry point.

One subcommand per stage, plus `check`, which validates a configuration and reports the
metadata/feature split without running anything. Every subcommand takes `-c/--config`,
which is what removes both the `/PATH/TO/...` editing ritual and the
run-from-the-script's-own-directory constraint of the original scripts.

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

#: Stages still to be built, and where they land. Kept explicit so `--help` cannot
#: imply a stage exists before it does.
PENDING = {
    "run": "every stage in sequence (available once the stages are)",
}

_EXTENSIONS = {"tab": ".tsv", "comma": ".csv", "semicolon": ".csv"}


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
        description=(
            "Feature selection for high-content screening by fitting "
            "concentration/dose-response models to earth mover's distance scores."
        ),
        epilog=(
            "Start from examples/published.ini, which reproduces the run described in "
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

    prune = _add(
        subparsers, "prune", "collapse near-redundant features by correlation (optional)"
    )
    prune.set_defaults(handler=_prune)

    subset = _add(
        subparsers, "subset", "write the input table restricted to the retained features"
    )
    subset.add_argument(
        "--features",
        metavar="FILE",
        help=(
            "feature list to apply, one name per line; defaults to the pruned list when "
            "[prune] is enabled and the selected list otherwise"
        ),
    )
    subset.set_defaults(handler=_subset)

    for command, description in PENDING.items():
        pending = _add(subparsers, command, f"{description} [not yet implemented]")
        pending.set_defaults(handler=_pending, command_name=command)

    return parser


def _add(subparsers, name: str, help_text: str) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(name, help=help_text, description=help_text)
    parser.add_argument(
        "-c",
        "--config",
        required=True,
        metavar="FILE",
        help="configuration file (.ini)",
    )
    return parser


def _pending(args: argparse.Namespace) -> int:
    print(
        f"error: `cdr-fs {args.command_name}` is not implemented yet - "
        f"{PENDING[args.command_name]}",
        file=sys.stderr,
    )
    return 3


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
        _field("levels", f"{_abbreviate(design.levels)}   (low to high)"),
    ]
    if design.dose is not None:
        lines.append(_field("dose", _abbreviate(f"{value:g}" for value in design.dose)))
    else:
        lines.append(_field("dose", "(none - levels are used by rank)"))
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
        _field(
            "slope_positive",
            f"{config.select.slope_positive}  (a positive linear slope on "
            f"{'at least one' if config.select.slope_positive == 'any' else 'every'} stratum)",
        ),
        _field(
            "nonconstant",
            f"{config.select.nonconstant}  (the constant model is not the best fit on "
            f"{'at least one' if config.select.nonconstant == 'any' else 'every'} stratum)",
        ),
        _field(
            "strata",
            ", ".join(config.select.strata)
            if config.select.strata
            else "(every stratum present in the data)",
        ),
        "[prune]",
        _field("enabled", str(config.prune.enabled).lower()),
    ]
    if config.prune.enabled:
        lines += [
            _field("threshold", f"|r| >= {config.prune.threshold:g}"),
            _field("linkage", config.prune.linkage),
            _field("representative", config.prune.representative),
            _field(
                "aggregate_by",
                ", ".join(config.prune.aggregate_by)
                if config.prune.aggregate_by
                else "(none - the rows are correlated as they are)",
            ),
            _field("fill_missing", config.prune.fill_missing),
        ]
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
    # Name the members of every small group: the two organelle counts among 469
    # rp_norm_* features are exactly the case where a lazy metadata pattern goes wrong.
    for head, count in breakdown:
        if count <= 5:
            lines.append(
                _field(f"  {head}*", ", ".join(resolved.features_with_prefix(head)))
            )
    return lines


# -------------------------------------------------------------------------- handlers


def _check(args: argparse.Namespace) -> int:
    from cdr_fs.schema import read_header

    config = load_config(args.config)
    print(f"cdr_FS {__version__} - configuration check")
    print(_field("config", str(config.path)))
    for line in _describe(config):
        print(line)

    columns = read_header(config.input.table, config.input.sep)
    warnings = config.validate_columns(columns)
    print()
    for line in _describe_schema(config, columns):
        print(line)

    if args.scan:
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
    if config.schema.group_by:
        usecols.append(config.schema.group_by)
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
    print()
    print(f"[data]  {len(frame):,} row(s)")
    print(_field(config.schema.condition, _abbreviate(sorted(levels), limit=12)))
    if strata is not None:
        print(_field(config.schema.group_by, _abbreviate(sorted(strata), limit=12)))
    return config.validate_observed(levels=levels, strata=strata)


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


def _write(config: Config, table, name: str):
    destination = Path(config.output.dir)
    destination.mkdir(parents=True, exist_ok=True)
    target = destination / f"{name}{_EXTENSIONS[config.input.sep_keyword]}"
    table.to_csv(target, sep=config.input.sep, index=False)
    print(f"wrote {target}  ({_human_bytes(target.stat().st_size)})")
    return target


def _list_path(config: Config, name: str) -> Path:
    """Where a feature list lives. One convention, so the stages can find each other."""
    return Path(config.output.dir) / f"{name}.txt"


def _write_list(config: Config, features: list[str], name: str) -> Path:
    destination = Path(config.output.dir)
    destination.mkdir(parents=True, exist_ok=True)
    target = _list_path(config, name)
    target.write_text("\n".join(features) + "\n", encoding="utf-8")
    print(f"wrote {target}  ({len(features)} feature(s))")
    return target


def _read_stage(config: Config, name: str, produced_by: str):
    """Read a previous stage's output, or say which command produces it."""
    import pandas as pd

    target = (
        Path(config.output.dir) / f"{name}{_EXTENSIONS[config.input.sep_keyword]}"
    )
    if not target.exists():
        raise ConfigError(
            f"{target} is missing - run `cdr-fs {produced_by} -c {config.path}` first"
        )
    return pd.read_csv(target, sep=config.input.sep, dtype={"feature": str, "stratum": str})


def _trim(args: argparse.Namespace) -> int:
    config = load_config(args.config)
    if not config.trim.enabled:
        print(
            "error: [trim] enabled is false, so there is nothing to do",
            file=sys.stderr,
        )
        return 3

    frame, _ = _prepare(config)
    if args.dry_run:
        print("dry run - nothing written")
        return 0
    _write(config, frame, "trimmed")
    return 0


def _emd(args: argparse.Namespace) -> int:
    from cdr_fs.emd import compute_baseline, compute_contrasts

    config = load_config(args.config)
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
    from cdr_fs.fit import fit_series

    config = load_config(args.config)
    distances = _read_stage(config, "emd", "emd")
    table, report = fit_series(config, distances)
    print(report.summary())
    _write(config, table, "fit")
    return 0


def _select(args: argparse.Namespace) -> int:
    from cdr_fs.select import select_features

    config = load_config(args.config)
    fits = _read_stage(config, "fit", "fit")
    retained, evidence, report = select_features(config, fits)
    print(report.summary())

    _write(config, evidence, "select_evidence")
    _write_list(config, retained, "selected")
    return 0


def _prune(args: argparse.Namespace) -> int:
    from cdr_fs.prune import prune_features
    from cdr_fs.subset import read_feature_list

    config = load_config(args.config)
    if not config.prune.enabled:
        print(
            "error: [prune] enabled is false, so there is nothing to do",
            file=sys.stderr,
        )
        return 3

    selected = _list_path(config, "selected")
    if not selected.exists():
        raise ConfigError(
            f"{selected} is missing - run `cdr-fs select -c {config.path}` first"
        )
    features = read_feature_list(selected)
    frame, _ = _prepare(config)

    kept, clusters, tree, report = prune_features(config, frame, features)
    print(report.summary())
    _write(config, clusters, "prune_clusters")
    _write(config, tree, "prune_linkage")
    _write_list(config, kept, "pruned")
    return 0


def _subset(args: argparse.Namespace) -> int:
    from cdr_fs.subset import read_feature_list, subset_table

    config = load_config(args.config)
    if args.features:
        source = Path(args.features)
        if not source.is_file():
            raise ConfigError(f"--features file not found: {source}")
    else:
        stage = "pruned" if config.prune.enabled else "selected"
        source = _list_path(config, stage)
        if not source.exists():
            raise ConfigError(
                f"{source} is missing - run `cdr-fs {'prune' if config.prune.enabled else 'select'}"
                f" -c {config.path}` first, or name a list with --features"
            )
    print(f"applying {source.name}")

    features = read_feature_list(source)
    frame, resolved = _prepare(config)
    subset, quality, report = subset_table(frame, resolved.metadata, features)
    print(report.summary())

    name = f"subset_{source.stem}"
    _write(config, subset, name)
    _write(config, quality, f"{name}_features")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
