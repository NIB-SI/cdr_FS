"""Golden regression for correlation pruning, against the published pruned list.

The published run's own output for this stage was never a small file, so the check is made
against the *composition* of the list instead of against the list itself. HCS-proc's feature
categorization step is run on exactly this stage's output - the pruned list, before the
missing-data filter its dimension-reduction step applies - and every bar segment of the
published figure is labelled with its count. That makes the figure a 4 x 7 table of counts,
transcribed here as `fixtures/published_categories_all_days.tsv`.

So the assertion is not "99 features", which one wrong merge could satisfy by accident, but
"99 features distributed across four organelles and seven measurement families exactly as
published". Reproducing the count and the composition takes the same list.

Needs `data/all_days_trimmed_features.txt` - the published run's own selected-and-trimmed
subset, 1.4 GB, from https://doi.org/10.5281/zenodo.17951792. Starting from that file rather
than from the raw table isolates this stage: its 182 feature columns are the published
selection, and its values are already trimmed, so anything that fails here is pruning.

    CDR_FS_GOLDEN=1 pytest tests/test_golden_prune.py -v
"""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd
import pytest
from cases import METADATA_PATTERNS, write_config
from tests_support import feature_list

from cdr_fs.config import load_config
from cdr_fs.prune import prune_features
from cdr_fs.schema import read_header, read_table, resolve_schema

FIXTURES = Path(__file__).resolve().parent / "fixtures"
DATA = Path(__file__).resolve().parents[1] / "data"
#: The published all-days selection applied back to the data: 9 metadata columns + 182
#: features x 503,920 objects, already trimmed to [p2.5, p97.5] per well.
INPUT = DATA / "all_days_trimmed_features.txt"
SELECTED = FIXTURES / "published_selected_all_days.txt"
CATEGORIES = FIXTURES / "hcs_proc_feature_categories.tsv"
PUBLISHED_COMPOSITION = FIXTURES / "published_categories_all_days.tsv"

EXPECTED = {
    "features_in": 182,
    "units": 531,  # (day, replicate, well)
    "kept": 99,
    "singletons": 65,
    "largest": 11,
    "after_missing_filter": 97,
}

pytestmark = [
    pytest.mark.skipif(
        not INPUT.exists(),
        reason=(
            f"{INPUT.name} absent from {DATA} - "
            "fetch from https://doi.org/10.5281/zenodo.17951792"
        ),
    ),
    pytest.mark.skipif(
        os.environ.get("CDR_FS_GOLDEN") != "1",
        reason="reads a 1.4 GB table; set CDR_FS_GOLDEN=1 to enable",
    ),
]

#: Display names, from HCS-proc's scripts/categorization/config.ini.
CATEGORY_OF = {
    "Counts": "Counts",
    "Number": "Counts",
    "AreaShape": "AreaShape",
    "Granularity": "Granularity",
    "Intensity": "Intensity",
    "RadialDistribution": "Radial Distribution",
    "Texture": "Texture",
    "Correlation": "Other/Correlation",
}
ORGANELLE_OF = {
    "Nuclei": "Nuclei (Hoecst)",
    "Lysosomes": "Lysosomes (Quinacrine)",
    "Mitochondria": "Mitochondria (TMRM)",
    "Other": "Other",
}


@pytest.fixture(scope="module")
def pruned(tmp_path_factory):
    """Run the pruning stage once, on the published trimmed subset."""
    directory = tmp_path_factory.mktemp("golden-prune")
    config = load_config(
        write_config(
            directory,
            {
                "schema.metadata_patterns": "\n".join(METADATA_PATTERNS),
                # The input is already trimmed; trimming it again would trim a trim.
                "trim.enabled": "false",
                "trim.scope": None,
                "prune.aggregate_by": "Metadata_Day,Metadata_Biorep,Metadata_Well",
            },
            table=INPUT,
        )
    )
    columns = read_header(INPUT, config.input.sep)
    schema = resolve_schema(columns, config.schema.compiled)
    frame = read_table(
        INPUT, config.input.sep, metadata=schema.metadata, features=schema.features
    )
    kept, clusters, tree, report = prune_features(config, frame, schema.features)
    return {
        "features": schema.features,
        "kept": kept,
        "clusters": clusters,
        "linkage": tree,
        "report": report,
        "frame": frame,
        "metadata": schema.metadata,
        "config": config,
    }


def test_the_input_is_the_published_selection(pruned):
    """Sanity check on the file: its feature columns are the 182 the published run selected."""
    assert set(pruned["features"]) == set(feature_list(SELECTED))
    assert len(pruned["features"]) == EXPECTED["features_in"]


def test_the_published_number_of_clusters_is_reproduced(pruned):
    report = pruned["report"]
    assert report.features == EXPECTED["features_in"]
    assert report.units == EXPECTED["units"]
    assert report.kept == EXPECTED["kept"]
    assert report.singletons == EXPECTED["singletons"]
    assert report.largest == EXPECTED["largest"]
    assert len(pruned["kept"]) == EXPECTED["kept"]


def test_the_published_composition_is_reproduced(pruned):
    """The real gate: the same features, category by category and organelle by organelle."""
    lookup = pd.read_csv(CATEGORIES, sep="\t", comment="#", dtype=str).fillna("")
    lookup = {row.feature: (row.category, row.organelle) for row in lookup.itertuples()}
    published = pd.read_csv(
        PUBLISHED_COMPOSITION, sep="\t", comment="#", index_col="organelle"
    ).fillna(0)

    counts = pd.DataFrame(0, index=published.index, columns=published.columns)
    unplaced = []
    for feature in pruned["kept"]:
        category, organelle = lookup.get(feature, ("", ""))
        column, row = CATEGORY_OF.get(category), ORGANELLE_OF.get(organelle)
        if column is None or row is None:
            unplaced.append(feature)
            continue
        counts.at[row, column] += 1

    assert unplaced == []
    assert counts.to_numpy().sum() == EXPECTED["kept"]
    pd.testing.assert_frame_equal(
        counts, published.astype(int), check_dtype=False, check_names=False
    )


def test_one_organelle_count_survives_pruning(pruned):
    """`counts_RelateLysoCell` is retained and `counts_RelateMitoCell` is not.

    The far end of the `counts_` trap: a `^counts_` metadata pattern would have dropped a
    feature the published run kept all the way through.
    """
    assert "counts_RelateLysoCell" in pruned["kept"]
    assert "counts_RelateMitoCell" not in pruned["kept"]


def test_every_cluster_has_exactly_one_representative(pruned):
    clusters = pruned["clusters"]
    per_cluster = clusters.groupby("cluster")["representative"].sum()
    assert (per_cluster == 1).all()
    assert clusters["cluster"].nunique() == EXPECTED["kept"]
    assert len(clusters) == EXPECTED["features_in"]
    assert set(clusters.loc[clusters["representative"], "feature"]) == set(pruned["kept"])


def test_how_loose_the_collapsing_actually_is(pruned):
    """Average linkage cuts on cluster means, so a member can sit past the pair threshold.

    Worth measuring rather than assuming: the 83 dropped features sit a median 0.052 from the
    representative that stands for them - |r| = 0.95 - but the loosest sits at 0.215, which
    is |r| = 0.785, comfortably outside the 0.9 the threshold names. That is what chaining
    does, and it is the number to look at before trusting one feature to speak for a cluster.
    """
    within = pruned["clusters"].query("not representative")["distance_to_representative"]
    assert len(within) == EXPECTED["features_in"] - EXPECTED["kept"]
    assert within.median() == pytest.approx(0.0524, abs=5e-4)
    assert within.max() == pytest.approx(0.2150, abs=5e-4)


def test_the_missing_data_filter_takes_two_more(pruned):
    """The default 30% rule, over the whole table: 99 features become 97.

    The two it removes are missing 47.0% and 95.2%, and the next thinnest feature is missing
    13.0% - so the threshold sits on a wide plateau rather than on a knife edge. Anything from
    14% to 46% gives the same 97. That is the number worth knowing about the choice: it is not
    delicately tuned.
    """
    from cdr_fs.subset import subset_table

    _, quality, report = subset_table(
        pruned["frame"],
        pruned["metadata"],
        pruned["kept"],
        drop_missing=True,
        max_missing=30,
    )
    assert report.matched == EXPECTED["kept"]
    assert report.kept == EXPECTED["after_missing_filter"]
    assert [name for name, _ in report.dropped] == [
        "rp_norm_AreaShape_FormFactor_RelateMitoCell",
        "rp_norm_Mean_PunctaLyso_Distance_Minimum_Cytoplasm_FilteredNuclei",
    ]
    surviving = quality[~quality["dropped"]]["nonmissing_fraction"]
    assert surviving.min() == pytest.approx(0.870, abs=5e-4)


def test_treating_the_object_indices_as_metadata_lands_on_the_published_list(pruned):
    """The route `examples/published.ini` takes: 175 -> 97 -> 95.

    `Number_Object_Number` is CellProfiler's within-image object label. Seven of its eight
    columns passed the concentration-response gate, and the published run carried them through
    pruning and then removed them by name downstream. Declaring them metadata instead removes
    them here - and they turn out to form exactly two pure clusters of the 99, one of four
    members and one of three, so nothing else about the clustering changes.

    What comes out is 95 features: 94 `rp_norm_*` plus `counts_RelateLysoCell`, which is the
    length and the composition of the published retained list. The two features the published
    run also struck out by hand are in it, and survive on their own merits.
    """
    from cdr_fs.subset import subset_table

    labels = [name for name in pruned["features"] if "Number_Object_Number" in name]
    assert len(labels) == 7
    # All seven sit in clusters made only of labels, so removing them removes whole clusters.
    clusters = pruned["clusters"]
    label_clusters = set(clusters.loc[clusters["feature"].isin(labels), "cluster"])
    members = clusters[clusters["cluster"].isin(label_clusters)]
    assert sorted(members.drop_duplicates("cluster")["size"]) == [3, 4]
    assert set(members["feature"]) == set(labels)

    measurements = [name for name in pruned["features"] if name not in set(labels)]
    kept, _, _, report = prune_features(pruned["config"], pruned["frame"], measurements)
    assert (report.features, report.kept) == (175, 97)

    _, _, final = subset_table(
        pruned["frame"], pruned["metadata"], kept, drop_missing=True, max_missing=30
    )
    assert final.kept == 95
    surviving = set(kept) - {name for name, _ in final.dropped}
    assert sum(name.startswith("rp_norm_") for name in surviving) == 94
    assert surviving - {name for name in surviving if name.startswith("rp_norm_")} == {
        "counts_RelateLysoCell"
    }


def test_the_linkage_table_describes_the_same_tree(pruned):
    table = pruned["linkage"]
    leaves = table[table["label"].notna()]
    merges = table[table["label"].isna()]
    assert leaves["label"].tolist() == list(pruned["features"])
    assert len(merges) == EXPECTED["features_in"] - 1
    assert merges["height"].is_monotonic_increasing
    # The cut that produced 99 clusters: 181 merges minus the 82 taken below it.
    assert int((merges["height"] < 0.1).sum()) == EXPECTED["features_in"] - EXPECTED["kept"]
