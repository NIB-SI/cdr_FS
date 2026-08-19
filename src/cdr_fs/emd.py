"""Earth mover's distance between populations, driven by a declared contrast set.

For each feature, each stratum and each contrast, this computes the 1-D Wasserstein
distance between the reference population's values and the compared population's values.
That is the measurement the whole selection rests on: a distance per exposure level, which
the fitter then treats as a concentration-response series.

Two contrast sets come out of one engine, because they differ only in how populations are
paired:

* **treatment contrasts** - control against each exposure level. This is the series that
  gets fitted.
* **the baseline set** - control against control, between biological replicates. This is
  the reproducibility floor: the distance two populations show when nothing was done to
  either of them, and therefore the yardstick for reading the treatment distances.

Both are returned as tidy long tables with the same columns, so downstream code and plots
treat them alike.

## Pooling

`pool_over` names the replicate column. With `per_replicate=False`, the published default,
all replicates at a level merge into one distribution and each contrast yields one
distance. With `per_replicate=True` the distance is computed within each replicate
separately, which puts batch variation between the points rather than inside the
distributions - and changes the numbers, so it no longer reproduces the article.

## Missing values

A population contributes only its non-missing values, so N varies per feature; both counts
are reported alongside every distance. When either side is empty the distance is undefined
and no row is emitted, matching the original scripts. Surviving infinities propagate into
the distance as NaN rather than being silently dropped - see `trim` on why they exist.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Iterable, Sequence

if TYPE_CHECKING:  # pragma: no cover
    import numpy as np
    import pandas as pd

    from cdr_fs.config import Config

__all__ = ["COLUMNS", "EmdReport", "compute_baseline", "compute_contrasts", "compute_emd"]

#: Column order of every table this module produces.
COLUMNS = (
    "feature",
    "stratum",
    "contrast",
    "group_a",
    "group_b",
    "replicate",
    "n_a",
    "n_b",
    "emd",
)


# eq=False: these hold numpy arrays, and a generated __eq__ or __hash__ over an
# array raises rather than answering. Identity is all that is ever needed.
@dataclass(frozen=True, eq=False)
class Population:
    """One side of a comparison: a label plus the rows it selects."""

    label: str
    rows: np.ndarray


# eq=False: these hold numpy arrays, and a generated __eq__ or __hash__ over an
# array raises rather than answering. Identity is all that is ever needed.
@dataclass(frozen=True, eq=False)
class Comparison:
    """A pair of populations to measure, tagged for the output table."""

    stratum: str
    contrast: str
    replicate: str
    a: Population
    b: Population


@dataclass(frozen=True)
class EmdReport:
    """What the engine measured, and what it could not."""

    comparisons: int
    features: int
    rows_emitted: int
    #: (feature, comparison) cells skipped because a population had no values at all.
    skipped_empty: int
    features_skipped: tuple[str, ...]
    #: Distances that came out non-finite, which happens when an infinity survived
    #: trimming. Such a feature cannot yield a complete series and will not be fitted.
    nonfinite: int
    features_nonfinite: tuple[str, ...]

    def summary(self) -> str:
        lines = [
            f"{self.rows_emitted:,} distance(s) over {self.comparisons} comparison(s) "
            f"x {self.features} feature(s)",
        ]
        if self.skipped_empty:
            lines.append(
                f"  {self.skipped_empty} (feature, comparison) cell(s) skipped - a "
                f"population held no values: {_listing(self.features_skipped)}"
            )
        if self.nonfinite:
            lines.append(
                f"  {self.nonfinite} non-finite distance(s) in "
                f"{len(self.features_nonfinite)} feature(s): "
                f"{_listing(self.features_nonfinite)}"
            )
            lines.append(
                "  an infinity survived trimming; these features cannot be fitted"
            )
        return "\n".join(lines)


def _listing(names: Sequence[str], limit: int = 3) -> str:
    shown = ", ".join(names[:limit])
    return shown if len(names) <= limit else f"{shown}, +{len(names) - limit} more"


# ------------------------------------------------------------------ population building


def _row_index(frame: pd.DataFrame, columns: Sequence[str]) -> dict[tuple[str, ...], np.ndarray]:
    """Positional row indices for every observed combination of `columns`."""
    import numpy as np

    if not columns:
        return {(): np.arange(len(frame), dtype=np.intp)}
    keys = list(zip(*(frame[column].to_numpy() for column in columns)))
    groups: dict[tuple[str, ...], list[int]] = {}
    for position, key in enumerate(keys):
        groups.setdefault(key, []).append(position)
    return {key: np.asarray(rows, dtype=np.intp) for key, rows in groups.items()}


def _strata(config: Config, frame: pd.DataFrame) -> list[str]:
    """The strata to work over, in the configured order where one is given."""
    present = list(dict.fromkeys(frame[config.schema.group_by])) if config.schema.group_by else [""]
    if config.select.strata is None:
        return present
    return [stratum for stratum in config.select.strata if stratum in present]


def _populations(
    config: Config, frame: pd.DataFrame
) -> dict[tuple[str, str, str], np.ndarray]:
    """(stratum, level, replicate) -> rows. `replicate` is "" for the pooled population."""
    import numpy as np

    columns = [config.schema.condition]
    if config.schema.group_by:
        columns.insert(0, config.schema.group_by)
    if config.schema.pool_over:
        columns.append(config.schema.pool_over)

    index = _row_index(frame, columns)
    populations: dict[tuple[str, str, str], list[np.ndarray]] = {}
    for key, rows in index.items():
        parts = dict(zip(columns, key))
        stratum = parts[config.schema.group_by] if config.schema.group_by else ""
        level = parts[config.schema.condition]
        replicate = parts[config.schema.pool_over] if config.schema.pool_over else ""
        # The pooled population is keyed by an empty replicate. Only add the per-replicate
        # entry when there is a replicate column, or the two keys coincide and every row
        # would be counted twice.
        if replicate:
            populations.setdefault((stratum, level, replicate), []).append(rows)
        populations.setdefault((stratum, level, ""), []).append(rows)

    return {
        key: np.sort(np.concatenate(parts)) for key, parts in populations.items()
    }


def contrast_comparisons(config: Config, frame: pd.DataFrame) -> list[Comparison]:
    """Control against each declared level, per stratum.

    With `per_replicate` the pairing happens inside each replicate; otherwise the
    replicates are pooled into one distribution per level.
    """
    populations = _populations(config, frame)
    replicates = (
        list(dict.fromkeys(frame[config.schema.pool_over]))
        if config.emd.per_replicate and config.schema.pool_over
        else [""]
    )

    comparisons = []
    for stratum in _strata(config, frame):
        for reference, level in config.emd.contrasts:
            for replicate in replicates:
                a = populations.get((stratum, reference, replicate))
                b = populations.get((stratum, level, replicate))
                if a is None or b is None:
                    continue
                comparisons.append(
                    Comparison(
                        stratum=stratum,
                        contrast=f"{reference}v{level}",
                        replicate=replicate,
                        a=Population(reference, a),
                        b=Population(level, b),
                    )
                )
    return comparisons


def baseline_comparisons(config: Config, frame: pd.DataFrame) -> list[Comparison]:
    """Every pair of replicates at the control level, per stratum.

    This is the reproducibility floor - the distance between two populations that were
    treated identically - so it is what a treatment distance has to exceed to mean
    anything.
    """
    if config.emd.baseline == "none" or not config.schema.pool_over:
        return []

    populations = _populations(config, frame)
    replicates = list(dict.fromkeys(frame[config.schema.pool_over]))
    control = config.design.control

    comparisons = []
    for stratum in _strata(config, frame):
        available = [
            replicate
            for replicate in replicates
            if (stratum, control, replicate) in populations
        ]
        for position, first in enumerate(available):
            for second in available[position + 1 :]:
                comparisons.append(
                    Comparison(
                        stratum=stratum,
                        contrast=f"{first}v{second}",
                        replicate="",
                        a=Population(f"{first}_{stratum}_{control}", populations[(stratum, control, first)]),
                        b=Population(f"{second}_{stratum}_{control}", populations[(stratum, control, second)]),
                    )
                )
    return comparisons


# -------------------------------------------------------------------------- the engine


def compute_emd(
    frame: pd.DataFrame,
    features: Sequence[str],
    comparisons: Iterable[Comparison],
) -> tuple[pd.DataFrame, EmdReport]:
    """Wasserstein distance for every (feature, comparison) pair.

    Looping features outermost keeps only one column in hand at a time - the full feature
    matrix of the published table is nearly 2 GB - and lets each population's values be
    extracted once per feature and reused, which matters because the control appears in
    every contrast of its stratum.
    """
    import numpy as np
    import pandas as pd
    from scipy.stats import wasserstein_distance

    features = list(features)
    comparisons = list(comparisons)

    records: list[tuple] = []
    skipped = 0
    skipped_features: dict[str, None] = {}
    nonfinite = 0
    nonfinite_features: dict[str, None] = {}

    for feature in features:
        column = frame[feature].to_numpy(dtype=np.float64, copy=False)
        extracted: dict[tuple[str, str, str], np.ndarray] = {}

        def values_of(stratum: str, population: Population, replicate: str) -> np.ndarray:
            key = (stratum, population.label, replicate)
            if key not in extracted:
                taken = column[population.rows]
                extracted[key] = taken[~np.isnan(taken)]
            return extracted[key]

        for comparison in comparisons:
            values_a = values_of(comparison.stratum, comparison.a, comparison.replicate)
            values_b = values_of(comparison.stratum, comparison.b, comparison.replicate)
            if values_a.size == 0 or values_b.size == 0:
                skipped += 1
                skipped_features.setdefault(feature)
                continue
            with np.errstate(invalid="ignore"):
                distance = wasserstein_distance(values_a, values_b)
            if not np.isfinite(distance):
                nonfinite += 1
                nonfinite_features.setdefault(feature)
            records.append(
                (
                    feature,
                    comparison.stratum,
                    comparison.contrast,
                    comparison.a.label,
                    comparison.b.label,
                    comparison.replicate,
                    values_a.size,
                    values_b.size,
                    distance,
                )
            )

    table = pd.DataFrame.from_records(records, columns=COLUMNS)
    report = EmdReport(
        comparisons=len(comparisons),
        features=len(features),
        rows_emitted=len(table),
        skipped_empty=skipped,
        features_skipped=tuple(skipped_features),
        nonfinite=nonfinite,
        features_nonfinite=tuple(nonfinite_features),
    )
    return table, report


def compute_contrasts(
    config: Config, frame: pd.DataFrame, features: Sequence[str]
) -> tuple[pd.DataFrame, EmdReport]:
    """The treatment contrast set: control against each level."""
    return compute_emd(frame, features, contrast_comparisons(config, frame))


def compute_baseline(
    config: Config, frame: pd.DataFrame, features: Sequence[str]
) -> tuple[pd.DataFrame, EmdReport]:
    """The baseline set: control against control, across replicates."""
    return compute_emd(frame, features, baseline_comparisons(config, frame))
