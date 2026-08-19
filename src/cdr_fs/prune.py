"""Optional correlation pruning: one representative per cluster of near-redundant features.

Concentration-response selection asks of each feature whether it responds. It cannot ask
whether two features are saying the same thing, and in morphological profiling many are:
area, perimeter and bounding-box area of the same object move together almost exactly. This
step collapses those groups, so that what follows - a UMAP, a distance measure, a figure -
is not dominated by whichever measurement family happens to be the most redundant.

The method:

1. **Aggregate** the object-level table to one row per experimental unit, by median. In the
   published run the unit is one well of one replicate on one day, which is
   `[prune] aggregate_by`. Correlating raw objects would measure within-well shape noise;
   correlating unit medians measures whether two features respond alike across the
   experiment.
2. **Correlate** every pair of features (Pearson) across those units.
3. **Cluster** on the distance `1 - |r|`, so that a strong *negative* correlation counts as
   the redundancy it is. Agglomerative, average linkage, cut at `1 - [prune] threshold`.
4. **Keep one member** of each cluster - by default the alphabetically first, which is what
   the published run did.

The clustering is kept separate from the dendrogram, as fitting is kept separate from
plotting: this module returns the linkage matrix as a table and `plots.py` draws it. The
original script computed both, but with *different* linkage methods - `average` for the
clusters and `median` for the picture - so the tree that was drawn was not the tree that was
cut. One linkage serves both here.

## Missing values, and why the default fills them

Trimming leaves NaN behind (see `trim.py`), so a unit median is taken over however many
objects survived, and a feature whose well was emptied has no median at all.
`[prune] fill_missing = column_mean` substitutes the feature's overall mean for every
missing *object-level* value before aggregating, which is what the published run did. It is
worth knowing what that does: for a feature with many trimmed values it pulls unit medians
toward the global mean, shrinking the between-unit spread that the correlation is computed
from. On the reference dataset the choice is not cosmetic - it moves the cluster count by
eight - so it is configurable and reported, with `none` computing each median from the
values that are actually there.

## A correlation that cannot be computed is treated as no correlation

Two features whose correlation comes out undefined - a constant series, a column left
entirely empty, an infinity that survived trimming - get distance 1.0, the maximum. They
therefore never merge with anything and are kept as singletons. That is the published
behaviour, and it errs in the safe direction (keep, rather than silently discard), but it
does mean an all-missing feature *survives* pruning. `PruneReport` names those features for
exactly that reason.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Sequence

if TYPE_CHECKING:  # pragma: no cover
    import numpy as np
    import pandas as pd

    from cdr_fs.config import Config

__all__ = [
    "CLUSTER_COLUMNS",
    "LINKAGE_COLUMNS",
    "Cluster",
    "PruneReport",
    "aggregate_units",
    "cluster_features",
    "correlation_distance",
    "linkage_table",
    "prune_features",
]

#: Column order of the cluster membership table.
CLUSTER_COLUMNS = (
    "cluster",
    "size",
    "feature",
    "representative",
    "distance_to_representative",
)
#: Column order of the linkage table. Leaf rows carry a label and no merge; merge rows carry
#: `left`, `right` and `height` and no label, so the tree and its leaf order travel together
#: in one file.
LINKAGE_COLUMNS = ("node", "label", "left", "right", "height", "size")


@dataclass(frozen=True)
class Cluster:
    """One group of near-redundant features, and the member kept for it."""

    features: tuple[str, ...]
    representative: str

    @property
    def size(self) -> int:
        return len(self.features)


@dataclass(frozen=True)
class PruneReport:
    """What pruning collapsed, in enough detail to put in a run manifest."""

    aggregate_by: tuple[str, ...]
    units: int
    fill_missing: str
    threshold: float
    linkage: str
    representative: str
    features: int
    kept: int
    singletons: int
    largest: int
    #: Feature pairs whose correlation was undefined; each was given distance 1.0.
    undefined_pairs: int
    #: Features with no data at all after aggregation. They survive as singletons.
    features_all_missing: tuple[str, ...]
    #: Features aggregating to a non-finite value, which makes their correlations undefined
    #: and so keeps them out of every cluster.
    features_nonfinite: tuple[str, ...]

    @property
    def removed(self) -> int:
        return self.features - self.kept

    def summary(self) -> str:
        unit = (
            f"{self.units} unit(s) of {', '.join(self.aggregate_by)}"
            if self.aggregate_by
            else f"{self.units} row(s), not aggregated"
        )
        filled = "" if self.fill_missing == "none" else ", missing values filled"
        lines = [
            f"pruned {self.features} to {self.kept} feature(s), one per cluster",
            f"  |r| >= {self.threshold:g} on median profiles over {unit}{filled}",
            f"  {self.linkage} linkage cut at distance {1 - self.threshold:g}, keeping the "
            f"{self.representative} member of each cluster",
            f"  {self.singletons} feature(s) clustered with nothing, largest cluster "
            f"{self.largest}, {self.removed} removed as redundant",
        ]
        if self.undefined_pairs:
            lines.append(
                f"  {self.undefined_pairs:,} pair(s) had no computable correlation and were "
                f"given the maximum distance, so they clustered with nothing"
            )
        if self.features_all_missing:
            lines.append(
                f"  {len(self.features_all_missing)} feature(s) are entirely missing after "
                f"aggregation and were kept as singletons: "
                f"{_listing(self.features_all_missing)}"
            )
        if self.features_nonfinite:
            lines.append(
                f"  {len(self.features_nonfinite)} feature(s) aggregate to a non-finite "
                f"value: {_listing(self.features_nonfinite)}"
            )
        return "\n".join(lines)


def _listing(names: Sequence[str], limit: int = 3) -> str:
    shown = ", ".join(names[:limit])
    return shown if len(names) <= limit else f"{shown}, +{len(names) - limit} more"


# ------------------------------------------------------------------------- aggregation


def aggregate_units(
    frame: pd.DataFrame,
    features: Sequence[str],
    by: Sequence[str],
    *,
    fill_missing: str = "column_mean",
) -> pd.DataFrame:
    """Median feature profile per experimental unit.

    `by` names the columns whose combinations define a unit; empty means no aggregation, and
    the rows are correlated as they are. `fill_missing = "column_mean"` replaces missing
    object-level values with the feature's overall mean *before* aggregating, which is the
    published behaviour - see the module docstring for what it costs.
    """
    features = list(features)
    by = list(by)
    missing = [column for column in (*features, *by) if column not in frame.columns]
    if missing:
        raise KeyError(f"columns absent from the table: {', '.join(missing)}")
    if fill_missing not in ("column_mean", "none"):
        raise ValueError(
            f"fill_missing must be 'column_mean' or 'none', got {fill_missing!r}"
        )

    values = frame[features]
    if fill_missing == "column_mean":
        values = values.fillna(values.mean())
    if not by:
        return values

    # dropna=False so a unit whose key is itself missing is aggregated rather than dropped.
    units = frame.groupby(by, sort=True, observed=True, dropna=False).ngroup()
    aggregated = values.groupby(units, sort=True).median()
    aggregated.index.name = "unit"
    return aggregated


# ------------------------------------------------------------------------- correlation


def correlation_distance(matrix: pd.DataFrame) -> tuple[np.ndarray, int]:
    """`1 - |Pearson r|` between the columns of `matrix`, as a dense distance matrix.

    Returns the matrix and the number of pairs whose correlation was undefined. Those are
    given the maximum distance 1.0, which is what keeps a constant or empty feature out of
    every cluster instead of merging it with an arbitrary partner.
    """
    import numpy as np

    correlation = matrix.corr(method="pearson").to_numpy(dtype=np.float64)
    undefined = int(np.isnan(correlation[np.triu_indices_from(correlation, k=1)]).sum())

    distance = 1.0 - np.abs(correlation)
    # NaN for an incomputable correlation. An infinity cannot come out of a correlation
    # coefficient, but a degenerate input has produced one on some pandas versions, so both
    # are pinned to the maximum distance rather than trusted.
    distance = np.nan_to_num(distance, nan=1.0, posinf=1.0, neginf=1.0)
    # Rounding can leave 1 - |r| a few ulp below zero, and a negative distance is not one.
    np.clip(distance, 0.0, None, out=distance)
    # A column with no computable self-correlation would otherwise sit at distance 1.0 from
    # itself. Nothing reads the diagonal, but leaving a lie there invites something that does.
    np.fill_diagonal(distance, 0.0)
    # Enforce exact symmetry: the condensed form keeps one triangle, and which one it keeps
    # must not be able to change the answer.
    return np.minimum(distance, distance.T), undefined


# -------------------------------------------------------------------------- clustering


def cluster_features(
    distance: np.ndarray,
    features: Sequence[str],
    cut: float,
    method: str = "average",
) -> tuple[list[tuple[int, ...]], np.ndarray]:
    """Agglomerative clustering of a precomputed distance matrix, cut at `cut`.

    Returns the clusters as tuples of column indices - ordered by their first member, so the
    numbering cannot depend on the order the merges happened to be recorded in - and the
    scipy linkage matrix, so that the same tree can be drawn.

    Two clusters merge only while their linkage distance is **strictly below** `cut`, so a
    pair sitting exactly at the threshold stays apart. That is the published rule
    (scikit-learn's `distance_threshold` is the distance at or above which clusters are not
    merged), and it is why the cut is expressed as a distance rather than as `|r|` throughout
    this module.
    """
    import numpy as np
    from scipy.cluster.hierarchy import linkage as scipy_linkage
    from scipy.spatial.distance import squareform

    count = len(features)
    if distance.shape != (count, count):
        raise ValueError(
            f"{count} feature name(s) for a distance matrix of shape {distance.shape}"
        )
    if count < 2:
        # squareform has no condensed form for a single point, and scipy's linkage rejects
        # an empty observation set.
        return [tuple(range(count))], np.empty((0, 4), dtype=np.float64)

    tree = scipy_linkage(squareform(distance, checks=False), method=method)

    # Union-find over the merges below the cut. Reading the tree directly, rather than
    # calling fcluster, is what makes the boundary rule ours: fcluster's 'distance'
    # criterion merges at exactly t.
    parent = list(range(count + len(tree)))

    def find(node: int) -> int:
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    for index, row in enumerate(tree):
        if row[2] >= cut:
            # Merge heights are non-decreasing, so nothing further down can still qualify.
            break
        root = find(int(row[0]))
        parent[find(int(row[1]))] = root
        parent[count + index] = root

    members: dict[int, list[int]] = {}
    for column in range(count):
        members.setdefault(find(column), []).append(column)
    return [tuple(group) for group in sorted(members.values())], tree


def linkage_table(tree: np.ndarray, features: Sequence[str]) -> pd.DataFrame:
    """The linkage matrix and its leaf labels as one table, for `plots.py` to draw.

    Nodes `0 .. n-1` are the leaves, carrying a `label` and no merge; nodes `n .. 2n-2` are
    the merges, carrying `left`, `right` and `height`. Keeping the leaf order in the same
    file is what lets the dendrogram be redrawn without recomputing a single correlation.
    """
    import pandas as pd

    count = len(features)
    leaves = pd.DataFrame(
        {
            "node": range(count),
            "label": list(features),
            "left": pd.NA,
            "right": pd.NA,
            "height": pd.NA,
            "size": 1,
        }
    )
    merges = pd.DataFrame(
        {
            "node": range(count, count + len(tree)),
            "label": pd.NA,
            "left": tree[:, 0].astype("int64"),
            "right": tree[:, 1].astype("int64"),
            "height": tree[:, 2],
            "size": tree[:, 3].astype("int64"),
        }
    )
    return pd.concat([leaves, merges], ignore_index=True)[list(LINKAGE_COLUMNS)]


# ------------------------------------------------------------------------------ stage


def prune_features(
    config: Config,
    frame: pd.DataFrame,
    features: Sequence[str],
) -> tuple[list[str], pd.DataFrame, pd.DataFrame, PruneReport]:
    """Cluster `features` by correlation and keep one representative per cluster.

    Returns the kept features in the input's own order, the cluster membership table, the
    linkage table, and a report.
    """
    import numpy as np
    import pandas as pd

    features = list(features)
    spec = config.prune
    units = aggregate_units(
        frame, features, spec.aggregate_by, fill_missing=spec.fill_missing
    )
    distance, undefined = correlation_distance(units)
    groups, tree = cluster_features(
        distance, features, cut=1.0 - spec.threshold, method=spec.linkage
    )

    clusters = []
    for group in groups:
        names = tuple(sorted(features[column] for column in group))
        # "first" is the input's own order, which is the table's column order.
        # "alphabetical" is what the published run used, by accident of a sorted
        # pandas Index.
        representative = (
            names[0] if spec.representative == "alphabetical" else features[group[0]]
        )
        clusters.append(Cluster(features=names, representative=representative))

    position = {name: index for index, name in enumerate(features)}
    rows = []
    for number, cluster in enumerate(clusters, start=1):
        centre = position[cluster.representative]
        for name in cluster.features:
            rows.append(
                {
                    "cluster": number,
                    "size": cluster.size,
                    "feature": name,
                    "representative": name == cluster.representative,
                    "distance_to_representative": distance[position[name], centre],
                }
            )
    membership = pd.DataFrame(rows, columns=list(CLUSTER_COLUMNS))

    kept = sorted(
        (cluster.representative for cluster in clusters), key=position.__getitem__
    )
    sizes = np.array([cluster.size for cluster in clusters])
    values = units.to_numpy(dtype=np.float64)
    report = PruneReport(
        aggregate_by=tuple(spec.aggregate_by),
        units=len(units),
        fill_missing=spec.fill_missing,
        threshold=spec.threshold,
        linkage=spec.linkage,
        representative=spec.representative,
        features=len(features),
        kept=len(kept),
        singletons=int((sizes == 1).sum()),
        largest=int(sizes.max()),
        undefined_pairs=undefined,
        features_all_missing=tuple(
            name for name, column in zip(features, np.isnan(values).T) if column.all()
        ),
        features_nonfinite=tuple(
            name for name, column in zip(features, np.isinf(values).T) if column.any()
        ),
    )
    return kept, membership, linkage_table(tree, features), report
