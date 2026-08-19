"""The stage-to-stage contract.

Every other module is tested directly, which says nothing about whether the stages can find
each other's output. That contract is entirely convention - a file name under `[output] dir` -
so it needs a test that runs the commands in order and looks at what appears on disk.
"""

from __future__ import annotations

from types import SimpleNamespace

import pandas as pd
import pytest
from cases import SUBSET, write_config
from tests_support import METADATA

from cdr_fs.cli import main


@pytest.fixture(scope="module")
def chain(tmp_path_factory):
    """Run the whole sequence once on the fixture table, and return the output directory."""
    directory = tmp_path_factory.mktemp("chain")
    config = write_config(directory, {"prune.aggregate_by": "Metadata_Well"})
    for stage in ("emd", "fit", "select", "prune", "subset"):
        assert main([stage, "-c", str(config)]) == 0, stage
    return SimpleNamespace(config=config, results=directory / "results")


def test_every_stage_leaves_its_output_where_the_next_one_looks(chain):
    produced = sorted(path.name for path in chain.results.iterdir())
    assert produced == [
        "emd.tsv",
        "emd_baseline.tsv",
        "fit.tsv",
        "prune_clusters.tsv",
        "prune_linkage.tsv",
        "pruned.txt",
        "select_evidence.tsv",
        "selected.txt",
        "subset_pruned.tsv",
        "subset_pruned_features.tsv",
        "subset_pruned_retained.txt",
    ]


def test_pruning_narrows_the_selected_list(chain):
    selected = (chain.results / "selected.txt").read_text(encoding="utf-8").split()
    pruned = (chain.results / "pruned.txt").read_text(encoding="utf-8").split()
    assert 0 < len(pruned) <= len(selected)
    assert set(pruned) <= set(selected)


def test_the_subset_is_named_after_the_list_it_applied(chain):
    """Two lists over one table must not overwrite each other's subset."""
    pruned = (chain.results / "pruned.txt").read_text(encoding="utf-8").split()
    columns = list(pd.read_csv(chain.results / "subset_pruned.tsv", sep="\t", nrows=1).columns)
    source = list(pd.read_csv(SUBSET, sep="\t", nrows=1).columns)
    # Every metadata column, the kept features, and nothing else.
    assert set(columns) == METADATA | set(pruned)
    # Metadata first, and both groups in the input table's own order.
    assert columns == [column for column in source if column in set(columns)]
    assert [column for column in columns if column in set(pruned)] == pruned


def test_the_retained_list_is_the_subset_header(chain):
    """Three lists, three counts: selected, pruned, and what survived the missing-data filter."""
    retained = (chain.results / "subset_pruned_retained.txt").read_text(encoding="utf-8").split()
    columns = list(pd.read_csv(chain.results / "subset_pruned.tsv", sep="\t", nrows=1).columns)
    assert [column for column in columns if column not in METADATA] == retained


def test_subset_can_be_pointed_at_any_list(chain, tmp_path):
    config = write_config(tmp_path, {"prune.aggregate_by": "Metadata_Well"})
    listing = tmp_path / "mine.txt"
    listing.write_text("counts_RelateLysoCell\n", encoding="utf-8")
    assert main(["subset", "-c", str(config), "--features", str(listing)]) == 0
    written = pd.read_csv(tmp_path / "results" / "subset_mine.tsv", sep="\t", nrows=1)
    assert "counts_RelateLysoCell" in written.columns
    assert "counts_RelateMitoCell" not in written.columns


def test_the_fit_table_is_self_consistent(chain):
    """Every stored parameter set must reproduce the AIC stored beside it.

    This is what `plots` relies on - a curve is drawn from `fit.tsv`, not refitted - so a
    parameter that does not read back exactly means a figure that contradicts its own legend.
    Six percent of these rows failed this when the parameters were written rounded.
    """
    import numpy as np

    from cdr_fs.config import load_config
    from cdr_fs.fit import series_from_emd
    from cdr_fs.models import MODEL_FUNCTIONS, information_criteria
    from cdr_fs.plots import parse_parameters

    config = load_config(chain.config)
    read = lambda name: pd.read_csv(  # noqa: E731
        chain.results / name, sep="\t", dtype={"feature": str, "stratum": str}
    )
    series = {
        (feature, stratum): (x, y)
        for feature, stratum, x, y, complete in series_from_emd(config, read("emd.tsv"))
        if complete
    }
    fits = read("fit.tsv")
    assert len(fits) > 100  # the test is worthless if the chain produced nothing

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


def test_the_plot_stage_draws_from_the_tables_on_disk(chain):
    assert main(["plot", "-c", str(chain.config), "--grid", "2"]) == 0
    drawn = sorted(path.name for path in chain.results.glob("*.png"))
    assert "dendrogram.png" in drawn
    assert "emd.png" in drawn and "emd_baseline.png" in drawn
    assert sum(name.startswith("fit_") for name in drawn) > 0
    assert all((chain.results / name).stat().st_size > 0 for name in drawn)


def test_plotting_says_what_it_could_not_draw(tmp_path, capsys):
    config = write_config(tmp_path, {"prune.aggregate_by": "Metadata_Well"})
    assert main(["plot", "-c", str(config)]) == 3
    printed = capsys.readouterr()
    assert "nothing to draw" in printed.err
    assert "run `cdr-fs emd`" in printed.err


def test_an_unknown_figure_name_is_refused(tmp_path, capsys):
    config = write_config(tmp_path, {"prune.aggregate_by": "Metadata_Well"})
    assert main(["plot", "-c", str(config), "--only", "fits,heatmap"]) == 2
    assert "unknown figure(s): heatmap" in capsys.readouterr().err


def test_a_stage_names_the_command_that_produces_its_input(tmp_path, capsys):
    config = write_config(tmp_path, {"prune.aggregate_by": "Metadata_Well"})
    assert main(["prune", "-c", str(config)]) == 2
    assert "run `cdr-fs select" in capsys.readouterr().err


def test_a_disabled_stage_says_so_rather_than_writing_nothing(tmp_path, capsys):
    config = write_config(tmp_path, {"prune.enabled": "false"})
    assert main(["prune", "-c", str(config)]) == 3
    assert "[prune] enabled is false" in capsys.readouterr().err


def test_an_unimplemented_stage_is_refused_not_ignored(tmp_path, capsys):
    config = write_config(tmp_path)
    assert main(["run", "-c", str(config)]) == 3
    assert "not implemented yet" in capsys.readouterr().err


def test_check_validates_without_touching_the_data(tmp_path, capsys):
    config = write_config(tmp_path)
    assert main(["check", "-c", str(config)]) == 0
    printed = capsys.readouterr().out
    assert "configuration is valid" in printed
    assert "aggregate_by" in printed
    assert not (tmp_path / "results").exists()
