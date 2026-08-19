"""The same code, a structurally different experiment.

Every other test in this suite is shaped by the RTgill-W1 design - four days, nine
concentrations, `Metadata_*` columns, tab-separated. That leaves one property of the code
unchecked, and it is the one the configuration file exists for: that the design comes from
the configuration and from nothing else.

`reshape.py` builds a table that shares its *values* with the published one and none of
its *structure*: no time axis, five exposure levels instead of nine, level labels that do
not sort into the response order, different column names on both sides of the split,
commas instead of tabs, interleaved columns and shuffled rows. This module runs the whole
chain on it and checks the answers are the ones the configuration asks for.

Two bugs came out of writing it, and both are regression-tested below:

* An experiment with no `[schema] group_by` has a single stratum labelled `""`. Written to
  a delimited file that is an empty field, and pandas reads an empty field back as NaN -
  which `groupby` then drops, silently, taking every row with it. `fit` produced an empty
  table, `select` reported "0 of 0", and nothing anywhere said why.
* `prune` on an empty feature list raised `IndexError` out of the clustering, because zero
  features had been treated as one empty cluster rather than as no clusters at all.

This says nothing about whether the method carries to another experiment. It is one dataset
in two shapes, not two datasets.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
import reshape
from cases import SUBSET
from reshape import CONTROL, LEVELS

from cdr_fs.cli import main
from cdr_fs.config import load_config
from cdr_fs.schema import read_stage_table

#: 20 features x 5 exposure levels, one stratum.
FEATURES = 20
#: The reshaped metadata: exposure, batch, sample, object_uid and three qc_ columns.
METADATA = 7


@pytest.fixture(scope="module")
def chain(tmp_path_factory):
    """Run every stage once on the reshaped table, and hold on to what it produced."""
    directory = tmp_path_factory.mktemp("reshaped")
    config = reshape.build(directory)
    for stage in ("check", "emd", "fit", "select", "prune", "subset", "plot"):
        assert main([stage, "-c", str(config)]) == 0, stage
    return SimpleNamespace(
        config=config,
        loaded=load_config(config),
        results=directory / "results",
    )


def read(chain, name: str) -> pd.DataFrame:
    return read_stage_table(chain.results / name, ",")


# --------------------------------------------------------------- the schema is read, not assumed


def test_the_reshaped_table_shares_no_column_name_with_the_published_one():
    """If it did, a leftover hardcoded name could still be doing the work."""
    published = set(pd.read_csv(reshape.SUBSET, sep="\t", nrows=0).columns)
    reshaped = set(reshape.reshape(pd.read_csv(reshape.SUBSET, sep="\t", nrows=8)).columns)
    assert published & reshaped == set()


def test_check_resolves_the_reshaped_schema(tmp_path, capsys):
    config = reshape.build(tmp_path)
    assert main(["check", "-c", str(config), "--scan"]) == 0
    printed = capsys.readouterr().out
    assert f"{METADATA + FEATURES} = {METADATA} metadata + {FEATURES} feature(s)" in printed
    # Two families, so the breakdown says something rather than naming one prefix.
    assert "alpha_* 10, beta_* 10" in printed
    assert "group_by          (none - a single stratum)" in printed
    assert "configuration is valid" in printed


def test_the_output_files_follow_the_configured_separator(chain):
    """`sep = comma` has to reach the written names, not only the reader."""
    produced = sorted(path.name for path in chain.results.iterdir() if path.suffix != ".png")
    assert produced == [
        "emd.csv",
        "emd_baseline.csv",
        "fit.csv",
        "prune_clusters.csv",
        "prune_linkage.csv",
        "pruned.txt",
        "select_evidence.csv",
        "selected.txt",
        "subset_pruned.csv",
        "subset_pruned_features.csv",
        "subset_pruned_retained.txt",
    ]


# ------------------------------------------------------------------------- the design is followed


def test_the_exposure_axis_follows_the_configured_order_not_the_alphabet(chain):
    """The sharp one. Sorted, the labels are high, low, max, mid, trace.

    `[design] levels` says trace, low, mid, high, max, and the fitted axis has to be that.
    A tool that sorted its levels would still produce a full series here, and a plausible
    one - it would simply be answering a different question.
    """
    from cdr_fs.fit import axis_positions, series_from_emd

    assert axis_positions(chain.loaded) == {level: float(i) for i, level in enumerate(LEVELS)}
    assert sorted(LEVELS) != list(LEVELS)  # otherwise this test proves nothing

    distances = read(chain, "emd.csv")
    by_level = {
        (row.feature, row.contrast.removeprefix(f"{CONTROL}v")): row.emd
        for row in distances.itertuples()
    }
    series = {feature: y for feature, _, _, y, _ in series_from_emd(chain.loaded, distances)}
    assert series
    for feature, y in series.items():
        assert list(y) == [by_level[(feature, level)] for level in LEVELS]


# ------------------------------------------------------------------------------- the two bugs


def test_a_single_unnamed_stratum_survives_the_round_trip(tmp_path):
    """An empty field is a label here, not a missing value.

    Read with a bare `pandas.read_csv` it comes back as NaN, and `groupby` drops NaN keys -
    so the table empties itself and only for the unstratified design.
    """
    path = tmp_path / "emd.csv"
    path.write_text("feature,stratum,emd\nalpha_01,,1.5\nbeta_01,,2.0\n", encoding="utf-8")

    assert pd.read_csv(path, dtype={"stratum": str})["stratum"].isna().all()

    table = read_stage_table(path, ",")
    assert list(table["stratum"]) == ["", ""]
    assert len(table.groupby(["feature", "stratum"])) == 2


def test_no_features_means_no_clusters():
    """Zero features used to come back as one empty cluster, which has no member to keep."""
    from cdr_fs.prune import cluster_features

    groups, tree = cluster_features(np.empty((0, 0)), [], cut=0.1, method="average")
    assert groups == []
    assert tree.shape == (0, 4)

    groups, tree = cluster_features(np.zeros((1, 1)), ["only"], cut=0.1, method="average")
    assert groups == [(0,)]
    assert tree.shape == (0, 4)


# ------------------------------------------------------------- the reshaped config is a real one


#: The configuration keys that carry the structure of the experiment. Every one of them
#: must differ from the published design, or this module is testing that design again
#: under a different name. `[trim] enabled` and `[prune] enabled` are deliberately absent:
#: both runs trim and prune, and that is a choice about method rather than about shape.
STRUCTURAL_KEYS = (
    "input.sep",
    "schema.metadata_patterns",
    "schema.condition",
    "schema.group_by",
    "schema.pool_over",
    "design.control",
    "design.levels",
    "design.exclude_from_fit",
    "trim.lower_percentile",
    "trim.upper_percentile",
    "trim.scope",
    "fit.models",
    "subset.drop_missing",
)


def test_the_fit_table_is_self_consistent(chain):
    """Every stored parameter set must reproduce the AIC and BIC stored beside it.

    This is what `plots` relies on - a curve is drawn from the fit table, not refitted - so
    a parameter that does not read back exactly means a figure that contradicts its own
    legend. Six percent of rows failed this when the parameters were written rounded to six
    digits: the models evaluate `log(x + 1e-10)`, a fitted `e` can settle a hair above
    -1e-10, and rounding collapses the argument of that logarithm onto zero.
    """
    import numpy as np

    from cdr_fs.fit import series_from_emd
    from cdr_fs.models import MODEL_FUNCTIONS, information_criteria
    from cdr_fs.plots import parse_parameters

    series = {
        (feature, stratum): (x, y)
        for feature, stratum, x, y, complete in series_from_emd(
            chain.loaded, read(chain, "emd.csv")
        )
        if complete
    }
    fits = read(chain, "fit.csv")
    assert len(fits) > 50  # the test is worthless if the chain produced nothing

    worst = 0.0
    for row in fits.itertuples():
        x, y = series[(row.feature, row.stratum)]
        parameters = parse_parameters(row.parameters)
        with np.errstate(over="ignore", invalid="ignore"):
            aic, bic = information_criteria(
                y - MODEL_FUNCTIONS[row.model](x, *parameters), len(parameters)
            )
        worst = max(worst, abs(aic - row.aic), abs(bic - row.bic))
    assert worst < 1e-9


def test_a_figure_that_refuses_to_draw_is_skipped_not_raised(chain, tmp_path, capsys):
    """`plot` draws several figures, so one unusable table must not lose the others.

    A figure function raises rather than draw something misleading - a distribution with no
    points reads as a result. That refusal has to reach the user as a skipped figure with a
    reason: `main` catches only `ConfigError`, so it surfaced as a traceback and exit 1.
    """
    import shutil

    config = reshape.build(tmp_path)
    results = tmp_path / "results"
    results.mkdir(exist_ok=True)
    for name in ("emd.csv", "fit.csv", "emd_baseline.csv"):
        shutil.copy(chain.results / name, results / name)
    # Header only: the shape `emd` used to write when a design had one replicate.
    header = (results / "emd_baseline.csv").read_text(encoding="utf-8").splitlines()[0]
    (results / "emd_baseline.csv").write_text(header + "\n", encoding="utf-8")

    assert main(["plot", "-c", str(config), "--only", "emd,baseline"]) == 0
    printed = capsys.readouterr()
    assert "skipped baseline - the distance table has no rows" in printed.err
    assert (results / "emd.png").exists()          # the other figure still drew
    assert not (results / "emd_baseline.png").exists()
