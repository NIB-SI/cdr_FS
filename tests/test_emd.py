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
    COLUMNS,
    baseline_comparisons,
    compute_baseline,
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


def test_table_shape_and_columns(tmp_path):
    frame = build(days=("T1", "T2", "T3"))
    config = configured(tmp_path, frame)
    table, report = compute_contrasts(config, frame, FEATURES)
    assert list(table.columns) == list(COLUMNS)
    # 3 strata x 3 contrasts x 2 features
    assert len(table) == 3 * 3 * 2
    assert report.skipped_empty == 0
    assert set(table["stratum"]) == {"T1", "T2", "T3"}


def test_every_declared_level_is_contrasted_including_withheld_ones(tmp_path):
    # EMD is measured for the whole series; only the fit skips withheld levels, so a
    # withheld level must still appear here. Four levels, because withholding one still
    # has to leave more points than the widest model has parameters.
    frame = build(levels=("ctrl", "low", "mid", "high", "top"))
    config = configured(
        tmp_path,
        frame,
        **{"design.levels": "low,mid,high,top", "design.exclude_from_fit": "top"},
    )
    table, _ = compute_contrasts(config, frame, ["shift"])
    assert set(table["contrast"]) == {"ctrlvlow", "ctrlvmid", "ctrlvhigh", "ctrlvtop"}
    assert config.design.fitted_levels == ("low", "mid", "high")


def test_works_without_a_stratification_column(tmp_path):
    frame = build(days=("T1",)).drop(columns=["Metadata_Day"])
    config = configured(tmp_path, frame, **{"schema.group_by": "", "select.strata": ""})
    table, _ = compute_contrasts(config, frame, ["shift"])
    assert len(table) == 3
    assert set(table["stratum"]) == {""}


def test_works_without_a_replicate_column(tmp_path):
    # No pool_over at all: there is nothing to pool over, so each population is just the
    # rows at that level. The counts have to say so - the pooled and per-replicate keys
    # coincide here, and a careless build counts every row twice.
    frame = build(reps=("R1",)).drop(columns=["Metadata_Biorep"])
    config = configured(
        tmp_path, frame, **{"schema.pool_over": "", "emd.baseline": "none"}
    )
    table, _ = compute_contrasts(config, frame, ["shift"])
    assert len(table) == 3
    assert set(table["n_a"]) == {40}
    assert set(table["n_b"]) == {40}
    assert table.set_index("contrast")["emd"]["ctrlvlow"] == pytest.approx(10.0)


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


def test_baseline_of_identical_replicates_is_zero(tmp_path):
    frame = build(reps=("R1", "R2", "R3"))
    config = configured(tmp_path, frame)
    table, _ = compute_baseline(config, frame, FEATURES)
    assert table["emd"].to_numpy() == pytest.approx(0.0)


def test_baseline_can_be_switched_off(tmp_path):
    frame = build()
    config = configured(tmp_path, frame, **{"emd.baseline": "none"})
    assert baseline_comparisons(config, frame) == []


def test_per_replicate_splits_each_contrast(tmp_path):
    frame = build(reps=("R1", "R2", "R3", "R4"))
    pooled = configured(tmp_path, frame)
    split = configured(tmp_path, frame, **{"emd.per_replicate": "true"})
    assert len(contrast_comparisons(pooled, frame)) == 3
    assert len(contrast_comparisons(split, frame)) == 12  # 3 contrasts x 4 replicates
    table, _ = compute_contrasts(split, frame, ["shift"])
    assert set(table["replicate"]) == {"R1", "R2", "R3", "R4"}


def test_pooling_can_cancel_an_effect_both_replicates_show(tmp_path):
    """Why `per_replicate` is a methodological choice and not a refactor.

    Here each replicate shows a clean shift of 10, but in opposite directions. Pooled, the
    two mixtures very nearly coincide and the distance collapses to 2.5; per replicate it
    is 10 in both. Pooling is not averaging - it can hide a response that every replicate
    individually shows.

    Note this cuts both ways, which is why the published default stays pooled: the same
    cancellation is what makes pooling robust to a plate that is merely shifted, and the
    upstream row/plate standardization exists to remove exactly that.
    """
    base = np.arange(40, dtype=float)
    offsets = {"ctrl": 0.0, "low": 1.0, "mid": 2.0, "high": 3.0}
    rows = []
    for replicate, direction in (("R1", +1.0), ("R2", -1.0)):
        for level, step in offsets.items():
            rows.append(
                pd.DataFrame(
                    {
                        "Concentration": level,
                        "Metadata_Day": "T1",
                        "Metadata_Biorep": replicate,
                        "Metadata_Well": f"{level}_{replicate}",
                        "shift": base + direction * step * 10.0,
                    }
                )
            )
    frame = pd.concat(rows, ignore_index=True)

    pooled_table, _ = compute_contrasts(configured(tmp_path, frame), frame, ["shift"])
    split_table, _ = compute_contrasts(
        configured(tmp_path, frame, **{"emd.per_replicate": "true"}), frame, ["shift"]
    )
    pooled = pooled_table.set_index("contrast")["emd"]["ctrlvlow"]
    per_replicate = split_table[split_table["contrast"] == "ctrlvlow"]["emd"]

    assert per_replicate.to_numpy() == pytest.approx([10.0, 10.0])
    assert pooled == pytest.approx(2.5)
    assert pooled < per_replicate.mean() / 3


def test_missing_level_replicate_combination_is_skipped_not_fatal(tmp_path):
    # The published data has exactly this hole: the top dose has no cells at all in one
    # replicate on one day, because it killed them.
    frame = build(reps=("R1", "R2"))
    frame = frame[~((frame["Concentration"] == "high") & (frame["Metadata_Biorep"] == "R2"))]
    config = configured(tmp_path, frame, **{"emd.per_replicate": "true"})
    comparisons = contrast_comparisons(config, frame)
    assert ("ctrlvhigh", "R2") not in {(c.contrast, c.replicate) for c in comparisons}
    assert ("ctrlvhigh", "R1") in {(c.contrast, c.replicate) for c in comparisons}
    # Pooled, the level still exists: it just has half the cells.
    pooled = contrast_comparisons(configured(tmp_path, frame), frame)
    assert "ctrlvhigh" in {c.contrast for c in pooled}


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


def test_empty_population_is_skipped_and_reported(tmp_path):
    frame = build()
    frame.loc[frame["Concentration"] == "mid", "shift"] = np.nan
    config = configured(tmp_path, frame)
    table, report = compute_contrasts(config, frame, ["shift", "flat"])
    assert "ctrlvmid" not in set(table[table["feature"] == "shift"]["contrast"])
    assert "ctrlvmid" in set(table[table["feature"] == "flat"]["contrast"])
    assert report.skipped_empty == 1
    assert report.features_skipped == ("shift",)


def test_surviving_infinity_gives_a_non_finite_distance(tmp_path):
    # Trimming is off here, so the infinity reaches the distance. It has to be visible
    # rather than quietly turning into a number.
    frame = build()
    frame.loc[frame.index[0], "shift"] = np.inf
    config = configured(tmp_path, frame)
    table, report = compute_contrasts(config, frame, ["shift"])
    assert report.nonfinite == 3
    assert report.features_nonfinite == ("shift",)
    assert not np.isfinite(table["emd"]).all()
    assert "cannot be fitted" in report.summary()


def test_distance_is_symmetric(tmp_path):
    frame = build()
    config = configured(tmp_path, frame)
    forward, _ = compute_contrasts(config, frame, ["shift"])
    reversed_config = configured(
        tmp_path, frame, **{"design.control": "high", "design.levels": "low,mid,ctrl"}
    )
    backward, _ = compute_contrasts(reversed_config, frame, ["shift"])
    assert forward.set_index("contrast")["emd"]["ctrlvhigh"] == pytest.approx(
        backward.set_index("contrast")["emd"]["highvctrl"]
    )
