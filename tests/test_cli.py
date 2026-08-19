"""The stage-to-stage contract.

Every other module is tested directly, which says nothing about whether the stages can find
each other's output. That contract is entirely convention - a file name under `[output] dir` -
so it needs a test that runs the commands in order and looks at what appears on disk.
"""

from __future__ import annotations

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
    return directory / "results"


def test_every_stage_leaves_its_output_where_the_next_one_looks(chain):
    produced = sorted(path.name for path in chain.iterdir())
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
    selected = (chain / "selected.txt").read_text(encoding="utf-8").split()
    pruned = (chain / "pruned.txt").read_text(encoding="utf-8").split()
    assert 0 < len(pruned) <= len(selected)
    assert set(pruned) <= set(selected)


def test_the_subset_is_named_after_the_list_it_applied(chain):
    """Two lists over one table must not overwrite each other's subset."""
    pruned = (chain / "pruned.txt").read_text(encoding="utf-8").split()
    columns = list(pd.read_csv(chain / "subset_pruned.tsv", sep="\t", nrows=1).columns)
    source = list(pd.read_csv(SUBSET, sep="\t", nrows=1).columns)
    # Every metadata column, the kept features, and nothing else.
    assert set(columns) == METADATA | set(pruned)
    # Metadata first, and both groups in the input table's own order.
    assert columns == [column for column in source if column in set(columns)]
    assert [column for column in columns if column in set(pruned)] == pruned


def test_the_retained_list_is_the_subset_header(chain):
    """Three lists, three counts: selected, pruned, and what survived the missing-data filter."""
    retained = (chain / "subset_pruned_retained.txt").read_text(encoding="utf-8").split()
    columns = list(pd.read_csv(chain / "subset_pruned.tsv", sep="\t", nrows=1).columns)
    assert [column for column in columns if column not in METADATA] == retained


def test_subset_can_be_pointed_at_any_list(chain, tmp_path):
    config = write_config(tmp_path, {"prune.aggregate_by": "Metadata_Well"})
    listing = tmp_path / "mine.txt"
    listing.write_text("counts_RelateLysoCell\n", encoding="utf-8")
    assert main(["subset", "-c", str(config), "--features", str(listing)]) == 0
    written = pd.read_csv(tmp_path / "results" / "subset_mine.tsv", sep="\t", nrows=1)
    assert "counts_RelateLysoCell" in written.columns
    assert "counts_RelateMitoCell" not in written.columns


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
