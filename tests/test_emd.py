"""The contrast-driven EMD engine.

Two things are being pinned here. First, the arithmetic: for shifted uniform samples the
Wasserstein distance is the shift, which is checkable by hand. Second, and more
importantly, the *pairing* - which populations get compared, what happens when one of them
does not exist, and what pooling does. The pairing is where a generalised engine earns its
keep, because the original scripts hardcoded it.

These configurations deliberately do not describe the RTgill-W1 design: three levels, no
time axis in some cases, renamed columns. If the engine only worked on the published shape,
these would fail.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from cases import write_config

from cdr_fs.config import load_config
from cdr_fs.emd import (
    baseline_comparisons,
    compute_contrasts,
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
