"""The reference data, reshaped into a structurally different experiment.

`examples/published.ini` is one experimental design, and every other test either uses it
or builds a table by hand. Neither shows the thing worth showing: that the design is read
from the configuration and nowhere else. A hardcoded `Metadata_Day`, an assumption that
levels sort numerically, a stratum label the output format cannot represent - all of those
survive a suite that only ever sees the RTgill-W1 layout.

So this module takes the committed slice of the published table and changes every
structural choice in it at once, leaving the measured values alone:

* **No time axis.** `Metadata_Day` is dropped rather than renamed, so four days pool into
  one stratum and `[schema] group_by` is empty. A single unnamed stratum is the case the
  output tables have to be able to round-trip.
* **Five exposure levels, not nine.** Consecutive concentrations are binned in pairs. Five
  fitted points cannot support BC5's five parameters, so `[fit] models` has to say so -
  which is itself the config carrying a structural fact.
* **Labels that do not sort into the response order.** The levels are named
  `trace, low, mid, high, max` and the control `vehicle`. Alphabetically that is
  `high, low, max, mid, trace`: anything that sorts the labels instead of reading
  `[design] levels` produces a different, wrong axis.
* **Different column names, metadata and features alike.** Nothing is called `Metadata_*`,
  `Concentration` or `rp_norm_*`. The features become two families, `alpha_*` and
  `beta_*`, which also gives `cdr-fs check`'s prefix breakdown something to report.
* **Different file shape.** Comma-separated rather than tab, columns interleaved so
  metadata does not come first, rows shuffled so no stage can rely on arrival order.

The values themselves are the published ones, so this is real data in a different shape -
which is the claim being tested. It is emphatically *not* a second biological dataset, and
the README says so.
"""

from __future__ import annotations

from pathlib import Path

from cases import SUBSET, write_config

__all__ = [
    "CONTROL",
    "LEVELS",
    "LEVEL_MAP",
    "METADATA",
    "RESHAPED",
    "build",
    "reshape",
]

#: Published label -> reshaped label. The five bins are consecutive pairs of the nine
#: exposure levels, keeping the low->high order; the top level stands alone.
LEVEL_MAP = {
    "11": "vehicle",
    "10": "trace",
    "9": "trace",
    "8": "low",
    "7": "low",
    "6": "mid",
    "5": "mid",
    "4": "high",
    "3": "high",
    "2": "max",
}
CONTROL = "vehicle"
#: Ordered low -> high. Sorted alphabetically this is high, low, max, mid, trace.
LEVELS = ("trace", "low", "mid", "high", "max")

#: The published table's metadata columns - the ones read as text rather than as floats.
_PUBLISHED_METADATA = (
    "Concentration",
    "counts_Cells",
    "counts_Cytoplasm",
    "counts_FilteredNuclei",
    "Metadata_Well",
    "Metadata_Day",
    "Metadata_Biorep",
    "Tech_replica",
    "Day_Well_BR",
    "cell_ID",
)


#: Published metadata column -> its name here. Columns absent from this map are dropped:
#: `Metadata_Day` because the time axis goes, `Tech_replica` and `Day_Well_BR` because they
#: are plate bookkeeping that a different experiment would not have.
METADATA = {
    "Concentration": "exposure",
    "Metadata_Biorep": "batch",
    "Metadata_Well": "sample",
    "cell_ID": "object_uid",
    "counts_Cells": "qc_cell_objects",
    "counts_Cytoplasm": "qc_cytoplasm_objects",
    "counts_FilteredNuclei": "qc_nucleus_objects",
}

#: A valid configuration for the reshaped table, in the shape `cases.config_text` renders.
#: Every value here differs from the published one, including the two percentiles and the
#: separator, so a default leaking through from `cases.BASE` would show up as a failure.
RESHAPED: dict[str, dict[str, str]] = {
    "input": {"table": "<table>", "sep": "comma"},
    "schema": {
        "metadata_patterns": "\n".join(
            ["^exposure$", "^batch$", "^sample$", "^object_uid$", "^qc_"]
        ),
        "condition": "exposure",
        "group_by": "",
        "pool_over": "batch",
    },
    "design": {
        "control": CONTROL,
        "levels": ",".join(LEVELS),
        "dose": "",
        "exclude_from_fit": "",
    },
    "trim": {
        "enabled": "true",
        "lower_percentile": "5",
        "upper_percentile": "95",
        "scope": "batch,sample",
    },
    "emd": {},
    # Five points cannot identify BC5's five parameters, and `load_config` refuses the
    # combination rather than fitting it - so a shorter series means a shorter model list.
    "fit": {"models": "BC4,LL4,WB1.4,Lin,Con"},
    "select": {},
    "prune": {"enabled": "true"},
    "subset": {"drop_missing": "true", "max_missing": "30"},
    "output": {"dir": "<output>"},
}


def reshape(frame):
    """Return `frame` with the published design replaced by the reshaped one.

    Deterministic: the row and column permutations come from a fixed seed, so a failure is
    reproducible and a diff of the written table is meaningful.
    """
    import numpy as np

    features = [column for column in frame.columns if column not in _PUBLISHED_METADATA]
    names = dict(METADATA)
    # Two families rather than one running index: `prefix_breakdown` groups by leading
    # token, and a single family would make that report vacuous.
    half = (len(features) + 1) // 2
    for position, feature in enumerate(features):
        family, index = ("alpha", position) if position < half else ("beta", position - half)
        names[feature] = f"{family}_{index + 1:02d}"

    reshaped = frame[[column for column in frame.columns if column in names]].rename(
        columns=names
    )
    reshaped["exposure"] = reshaped["exposure"].map(LEVEL_MAP)

    generator = np.random.default_rng(20260819)
    reshaped = reshaped.iloc[generator.permutation(len(reshaped))].reset_index(drop=True)
    columns = list(reshaped.columns)
    return reshaped[[columns[index] for index in generator.permutation(len(columns))]]


def build(directory: Path, overrides: dict[str, str | None] | None = None):
    """Write the reshaped table and a configuration for it into `directory`.

    Returns the configuration's path; the table sits beside it as `reshaped.csv`.
    """
    import pandas as pd

    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    table = directory / "reshaped.csv"
    reshape(pd.read_csv(SUBSET, sep="\t", dtype={name: str for name in _PUBLISHED_METADATA}))\
        .to_csv(table, index=False)
    return write_config(directory, overrides, table=table, base=RESHAPED)

