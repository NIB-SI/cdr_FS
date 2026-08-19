"""The contrast-driven EMD engine.

Two things are being pinned here. First, the arithmetic: for shifted uniform samples the
Wasserstein distance is the shift, which is checkable by hand. Second, and more
importantly, the *pairing* - which populations get compared, what happens when one of them
does not exist, and what pooling does. The pairing is what the original scripts hardcoded,
so it is what has to be read from the configuration now.

These configurations deliberately do not describe the RTgill-W1 design: three levels, no
time axis in some cases, renamed columns. If the engine only worked on the published shape,
these would fail.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from cases import COLUMNS_PUBLISHED, METADATA_PATTERNS, write_config

from cdr_fs.config import load_config
from cdr_fs.emd import (
    baseline_comparisons,
    compute_contrasts,
    contrast_comparisons,
)

SIMPLE = {
    "schema.metadata_patterns": "^Concentration$\n^Metadata_",
    "schema.condition": "Concentration",
    "schema.group_by": "Metadata_Day",
    "schema.pool_over": "Metadata_Biorep",
    "design.control": "ctrl",
    "design.levels": "low,mid,high",
    "design.exclude_from_fit": "",
    "design.dose": "",
    "trim.enabled": "false",
    "fit.models": "Lin,Con",
    "select.strata": "",
    "prune.enabled": "false",
}


def configured(tmp_path, frame, **extra):
    """Write `frame` to a table and return a Config describing it."""
    table = tmp_path / "table.tsv"
    frame.to_csv(table, sep="\t", index=False)
    return load_config(write_config(tmp_path, {**SIMPLE, **extra}, table=table))


def build(levels=("ctrl", "low", "mid", "high"), days=("T1",), reps=("R1", "R2"), n=40):
    """A tidy table with one feature that shifts by 10 per level."""
    offsets = {level: index * 10.0 for index, level in enumerate(levels)}
    rows = []
    for day in days:
        for level in levels:
            for rep in reps:
                rows.append(
                    pd.DataFrame(
                        {
                            "Concentration": level,
                            "Metadata_Day": day,
                            "Metadata_Biorep": rep,
                            "Metadata_Well": f"{level}_{rep}",
                            "shift": np.arange(n, dtype=float) + offsets[level],
                            "flat": np.arange(n, dtype=float),
                        }
                    )
                )
    return pd.concat(rows, ignore_index=True)


FEATURES = ["shift", "flat"]


def test_distance_of_a_pure_shift_is_the_shift(tmp_path):
    frame = build()
    config = configured(tmp_path, frame)
    table, _ = compute_contrasts(config, frame, FEATURES)
    shift = table[table["feature"] == "shift"].set_index("contrast")["emd"]
    assert shift["ctrlvlow"] == pytest.approx(10.0)
    assert shift["ctrlvmid"] == pytest.approx(20.0)
    assert shift["ctrlvhigh"] == pytest.approx(30.0)
    flat = table[table["feature"] == "flat"]["emd"]
    assert flat.to_numpy() == pytest.approx(0.0)


def test_baseline_pairs_every_replicate_at_the_control(tmp_path):
    frame = build(reps=("R1", "R2", "R3", "R4"))
    config = configured(tmp_path, frame)
    comparisons = baseline_comparisons(config, frame)
    assert len(comparisons) == 6  # 4 choose 2
    assert {c.contrast for c in comparisons} == {
        "R1vR2", "R1vR3", "R1vR4", "R2vR3", "R2vR4", "R3vR4",
    }
    # Labels carry replicate, stratum and level, as the published table does.
    assert comparisons[0].a.label == "R1_T1_ctrl"


def test_only_non_missing_values_contribute_and_counts_say_so(tmp_path):
    frame = build()
    frame.loc[frame["Concentration"] == "low", "shift"] = np.where(
        np.arange((frame["Concentration"] == "low").sum()) < 10,
        np.nan,
        frame.loc[frame["Concentration"] == "low", "shift"],
    )
    config = configured(tmp_path, frame)
    table, _ = compute_contrasts(config, frame, ["shift"])
    row = table[table["contrast"] == "ctrlvlow"].iloc[0]
    assert row["n_a"] == 80  # control untouched: 2 replicates x 40
    assert row["n_b"] == 70  # ten values dropped
    assert row["n_b"] < row["n_a"]


def test_the_published_grid_is_the_shape_the_configuration_asks_for(tmp_path):
    """The size of both published EMD tables, derived without reading a byte of the data.

    16,946 and 11,292 are otherwise asserted only by `test_golden.py`, which needs 3.9 GB
    and skips by default - so a change to the contrast set or the baseline pairing would sit
    unnoticed until someone downloaded the dataset. Both numbers are grid sizes: features
    times strata times comparisons, less the cells where a population was empty. Everything
    but the feature count comes from the design, and the feature count comes from the
    committed header.

    A skeleton frame is enough. The comparisons depend on which combinations of day, level
    and replicate exist, not on how many objects sit in each.
    """
    import itertools

    from cdr_fs.schema import resolve_schema

    columns = COLUMNS_PUBLISHED.read_text(encoding="utf-8").split()
    features = len(resolve_schema(columns, METADATA_PATTERNS).features)
    assert features == 471  # the published split, before the object indices became metadata

    days = ["D1", "D5", "D7", "D9"]
    replicates = ["BR1", "BR2", "BR3", "BR4"]
    levels = ["11", "10", "9", "8", "7", "6", "5", "4", "3", "2"]
    skeleton = pd.DataFrame(
        [
            {"Metadata_Day": day, "Metadata_Biorep": rep, "Concentration": level}
            for day, rep, level in itertools.product(days, replicates, levels)
        ]
    )
    config = load_config(write_config(tmp_path))

    # Control against each of the nine levels, on each of the four days. The withheld top
    # dose is still measured - only the fit skips it - so it is nine, not eight.
    treatment = len(contrast_comparisons(config, skeleton))
    assert treatment == 4 * 9
    # Every unordered pair of the four replicates, at the control, on each day.
    baseline = len(baseline_comparisons(config, skeleton))
    assert baseline == 4 * 6

    # The published tables sit just under the full grids: the shortfall is the
    # (feature, comparison) cells where one side had no values left after trimming.
    published_treatment, published_baseline = 16_946, 11_292
    assert features * treatment - published_treatment == 10
    assert features * baseline - published_baseline == 12
