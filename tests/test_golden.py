"""Golden regression against the published run.

This is the test that decides whether the extraction is trustworthy: the same input, the
same configuration, and the numbers the article's pipeline actually produced.

It needs the 3.7 GB per-object table and the two published EMD tables, none of which can
live in git, so it skips with instructions when they are absent. It also skips when they
*are* present unless `CDR_FS_GOLDEN=1` is set, because the pipeline takes minutes and the
rest of the suite runs in seconds.

    CDR_FS_GOLDEN=1 pytest tests/test_golden.py -v

Required in `data/`, all from https://doi.org/10.5281/zenodo.17951792:

* `cell_ID_pooled_median_row_plate_standardization_cid.txt` - the untrimmed input
* `EMD_conc_2.5_97.5_well.txt` - published treatment contrasts
* `EMD_c11_2.5_97.5_well.txt` - published baseline set

## On the tolerance

Distances are compared to 1e-9 relative, not bit-for-bit. Population sizes are compared
exactly, and that is the sharper test: a count is an integer, so all 28,238 of them
matching means the trim step reproduced the published one value for value, with no
tolerance to hide behind.

The distances themselves agree to 8.5e-13 relative, with about a quarter bit-identical.
Each one is a sum of tens of thousands of non-negative terms, so summation order alone
admits a relative error up to N*eps, around 9e-12 here - that bound is the reason a
tolerance is correct rather than lax. Requiring bit-equality would tie the suite to one
scipy version and one history of how the intermediate table was written.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from cases import METADATA_PATTERNS, write_config
from tests_support import as_fit_table

from cdr_fs.config import load_config
from cdr_fs.emd import compute_baseline, compute_contrasts
from cdr_fs.schema import read_header, read_table, resolve_schema
from cdr_fs.trim import trim_extremes

DATA = Path(__file__).resolve().parents[1] / "data"
INPUT = DATA / "cell_ID_pooled_median_row_plate_standardization_cid.txt"
PUBLISHED_CONTRASTS = DATA / "EMD_conc_2.5_97.5_well.txt"
PUBLISHED_BASELINE = DATA / "EMD_c11_2.5_97.5_well.txt"
#: The published fit table, from the run behind the article - so BC4 is still the LL4
#: duplicate it was before the correction. Optional: only the selection tests use it.
PUBLISHED_FITS = DATA / "model_fit_results.txt"

TOLERANCE = 1e-9

# Measured from the published outputs, not quoted from the article.
EXPECTED = {
    "rows": 503_920,
    "columns": 481,
    "metadata": 10,
    "features": 471,
    "trim_groups": 531,
    "contrast_rows": 16_946,
    "baseline_rows": 11_292,
    # The two AreaShape_FormFactor organelle features carry infinities, so some of their
    # cells are emptied by trimming and some of their distances are non-finite.
    "contrast_nonfinite": 13,
    "baseline_nonfinite": 3,
}

missing = [path.name for path in (INPUT, PUBLISHED_CONTRASTS, PUBLISHED_BASELINE) if not path.exists()]
pytestmark = [
    pytest.mark.skipif(
        bool(missing),
        reason=(
            f"golden data absent from {DATA}: {', '.join(missing)} - "
            "fetch from https://doi.org/10.5281/zenodo.17951792"
        ),
    ),
    pytest.mark.skipif(
        os.environ.get("CDR_FS_GOLDEN") != "1",
        reason="golden run takes minutes; set CDR_FS_GOLDEN=1 to enable",
    ),
]


@pytest.fixture(scope="module")
def published_run(tmp_path_factory):
    """Read, trim and measure the published input once for the whole module."""
    directory = tmp_path_factory.mktemp("golden")
    config = load_config(
        write_config(
            directory,
            {"schema.metadata_patterns": "\n".join(METADATA_PATTERNS)},
            table=INPUT,
        )
    )
    columns = read_header(INPUT, config.input.sep)
    config.validate_columns(columns)
    schema = resolve_schema(columns, config.schema.compiled)

    frame = read_table(
        INPUT, config.input.sep, metadata=schema.metadata, features=schema.features
    )
    frame, trim_report = trim_extremes(
        frame,
        schema.features,
        config.trim.scope,
        config.trim.lower_percentile,
        config.trim.upper_percentile,
        inplace=True,
    )
    contrasts, contrast_report = compute_contrasts(config, frame, schema.features)
    baseline, baseline_report = compute_baseline(config, frame, schema.features)
    return {
        "columns": columns,
        "schema": schema,
        "trim": trim_report,
        "contrasts": contrasts,
        "contrast_report": contrast_report,
        "baseline": baseline,
        "baseline_report": baseline_report,
    }


#: which of the two contrast sets -> (published file, expected row count, expected NaNs)
SETS = {
    "contrasts": (PUBLISHED_CONTRASTS, "contrast_rows", "contrast_nonfinite"),
    "baseline": (PUBLISHED_BASELINE, "baseline_rows", "baseline_nonfinite"),
}


def published_for(which):
    """The published table for one contrast set, restricted to the pairs we declare."""
    path = SETS[which][0]
    frame = load_published(path)
    if which == "contrasts":
        # The published table holds all 45 level pairs per day; our contrast set declares
        # control against each level, so keep only those.
        wanted = {frozenset(("11", str(level))) for level in range(2, 11)}
        frame = frame[
            [
                frozenset((a, b)) in wanted
                for a, b in zip(frame["Population1"], frame["Population2"])
            ]
        ]
    return frame


def load_published(path):
    frame = pd.read_csv(
        path,
        sep="\t",
        dtype={"Feature": str, "Population1": str, "Population2": str, "Metadata_Day": str},
    )
    frame["key"] = [
        (feature, day, frozenset((a, b)))
        for feature, day, a, b in zip(
            frame["Feature"], frame["Metadata_Day"], frame["Population1"], frame["Population2"]
        )
    ]
    return frame


def joined(mine, published):
    mine = mine.copy()
    mine["key"] = [
        (feature, stratum, frozenset((a, b)))
        for feature, stratum, a, b in zip(
            mine["feature"], mine["stratum"], mine["group_a"], mine["group_b"]
        )
    ]
    merged = mine.merge(published, on="key", how="outer", indicator=True)
    assert (merged["_merge"] == "both").all(), (
        f"{(merged['_merge'] != 'both').sum()} row(s) present on only one side, e.g. "
        f"{merged[merged['_merge'] != 'both']['key'].head(3).tolist()}"
    )
    return merged


# ------------------------------------------------------------------------- the schema


def test_schema_matches_the_published_table(published_run):
    schema = published_run["schema"]
    assert len(published_run["columns"]) == EXPECTED["columns"]
    assert len(schema.metadata) == EXPECTED["metadata"]
    assert len(schema.features) == EXPECTED["features"]
    assert schema.prefix_breakdown() == [("rp_", 469), ("counts_", 2)]
    assert "counts_RelateLysoCell" in schema.feature_set
    assert "counts_RelateMitoCell" in schema.feature_set


def test_trim_covers_the_published_design(published_run):
    report = published_run["trim"]
    assert report.rows == EXPECTED["rows"]
    assert report.features == EXPECTED["features"]
    # One well of one replicate on one day; each well holds a single concentration.
    assert report.groups == EXPECTED["trim_groups"]
    assert report.rows_with_missing_scope == 0
    # Only the two AreaShape_FormFactor organelle features carry infinities.
    assert set(report.features_with_nonfinite) == {
        "rp_norm_AreaShape_FormFactor_RelateLysoCell",
        "rp_norm_AreaShape_FormFactor_RelateMitoCell",
    }


# --------------------------------------------------------------------- the distances


@pytest.mark.parametrize("which", list(SETS))
def test_distances_reproduce_the_published_tables(published_run, which):
    _, rows_key, nonfinite_key = SETS[which]
    mine = published_run[which]
    assert len(mine) == EXPECTED[rows_key]
    merged = joined(mine, published_for(which))

    theirs = merged["EMD_score"].to_numpy(float)
    ours = merged["emd"].to_numpy(float)

    # Non-finite distances must land in exactly the same cells.
    assert np.array_equal(np.isnan(ours), np.isnan(theirs))
    assert int(np.isnan(ours).sum()) == EXPECTED[nonfinite_key]

    finite = ~np.isnan(ours)
    scale = np.maximum(np.abs(ours[finite]), np.abs(theirs[finite]))
    relative = np.abs(ours[finite] - theirs[finite]) / np.where(scale > 0, scale, 1.0)
    assert relative.max() < TOLERANCE, f"worst relative difference {relative.max():.3e}"


@pytest.mark.parametrize("which", list(SETS))
def test_population_sizes_reproduce_exactly(published_run, which):
    """The sharper of the two comparisons: integers, so no tolerance can hide a difference.

    A population size is the number of values left in that (feature, population) cell after
    trimming. Every one of them matching means the trim step reproduced the published one
    exactly - including the cells that infinities emptied.
    """
    merged = joined(published_run[which], published_for(which))

    swapped = merged["group_a"] != merged["Population1"]
    theirs_a = np.where(swapped, merged["Count2"], merged["Count1"])
    theirs_b = np.where(swapped, merged["Count1"], merged["Count2"])
    assert np.array_equal(merged["n_a"].to_numpy(), theirs_a)
    assert np.array_equal(merged["n_b"].to_numpy(), theirs_b)


# ------------------------------------------ the fit, against the published fit table
#
# The selection half of this comparison lives in test_golden_selection.py, which needs
# only the 1 MB fit table and so runs without opting in.


def published_fit_table():
    """The published fit table, in this package's column layout."""
    return as_fit_table(PUBLISHED_FITS)


@pytest.mark.skipif(
    not PUBLISHED_FITS.exists(),
    reason=f"{PUBLISHED_FITS.name} absent from {DATA}",
)
def test_the_deterministic_models_match_the_published_fit(published_run, tmp_path):
    """Con and Lin agree to floating-point noise; the sigmoids cannot be expected to.

    Least squares for a constant and for a straight line has a closed form, so those two
    columns are a clean check on the whole chain that feeds them - the trim, the distances,
    the eight-point series, the AIC/BIC expressions. They agree to 1e-12.

    The four sigmoid models are fitted iteratively from the same starting values with the
    same evaluation budget, and land in slightly different places on a flat surface under a
    different scipy version. On the reference data that difference is worth exactly one
    feature out of 471, and BC4 differs outright because the published run's BC4 was the
    LL4 duplicate.
    """
    from cdr_fs.fit import fit_series
    from cdr_fs.config import load_config

    config = load_config(write_config(tmp_path, table=INPUT))
    mine, _ = fit_series(config, published_run["contrasts"])
    merged = mine.merge(
        published_fit_table(), on=["feature", "stratum", "model"], suffixes=("_mine", "_pub")
    )

    for model in ("Con", "Lin"):
        rows = merged[merged["model"] == model]
        assert len(rows) == 1880, model
        ours = rows["aic_mine"].to_numpy(float)
        theirs = rows["aic_pub"].to_numpy(float)
        finite = np.isfinite(ours) & np.isfinite(theirs)
        scale = np.maximum(np.abs(ours[finite]), np.abs(theirs[finite]))
        relative = np.abs(ours[finite] - theirs[finite]) / np.where(scale > 0, scale, 1.0)
        assert relative.max() < 1e-12, f"{model}: worst {relative.max():.3e}"
