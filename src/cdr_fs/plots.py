"""Diagnostic figures, drawn from the tables the other stages wrote.

Nothing here computes a distance, fits a model or clusters a feature. Every figure is a view
of a file: `emd.tsv` and `fit.tsv` for the fit panels, `emd_baseline.tsv` for the
reproducibility floor, `correlation_linkage.tsv` and `correlation_clusters.tsv` for the
dendrogram. That is the whole point of the split - in the original pipeline one function
fitted, accumulated results and drew, so a figure could disagree with the table beside it and
nothing would notice. Here a figure that disagrees with its table is impossible, because the
table is where the figure comes from.

Three figures:

* **Fit panels** - one panel per feature and stratum: the distance points and the six fitted
  curves, with a legend ordered by information criterion. This is the figure the article's
  Figure 4 was composed from.
* **Distance distribution** - every distance of a table, per feature, features ordered by
  median. Run on `emd_baseline.tsv` it is the between-replicate reproducibility floor; run
  on `emd.tsv` it is the treatment distances to read against that floor.
* **Dendrogram** - the tree the correlation stage cut, with the cut drawn on it.

## Curves are drawn through the fitted points, not on a dense grid

A model curve here is the fitted function evaluated at the exposure levels and joined up,
which is what the published figures show. It is also the honest choice under `x_scale =
rank`, where the lowest level sits at x = 0 and the sigmoid models evaluate `log(x + 1e-10)`:
on a dense grid the fitted curve swings through several orders of magnitude between the first
two levels, which says more about the axis than about the feature. The parameters are in
`fit.tsv` for anyone who wants to evaluate them anywhere else.

## No pyplot

Figures are built through `matplotlib.figure.Figure` rather than `pyplot`, so there is no
global current-figure state to leak between calls and no backend to select. A caller that
draws four hundred pages in a loop does not have to remember to close them.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Iterable, Sequence

if TYPE_CHECKING:  # pragma: no cover
    import pandas as pd
    from matplotlib.axes import Axes

    from cdr_fs.config import Config

__all__ = [
    "MODEL_COLOURS",
    "SPLIT_PERCENTILE",
    "parse_parameters",
    "plot_dendrogram",
    "plot_distribution",
    "plot_fit_panels",
]

#: Model -> colour, as published: ColorBrewer RdYlBu, warm for the hormesis models through
#: to cool for the null ones.
MODEL_COLOURS = {
    "BC4": "#d73027",
    "BC5": "#f46d43",
    "LL4": "#fdae61",
    "WB1.4": "#abd9e9",
    "Lin": "#74add1",
    "Con": "#4575b4",
}

#: Where the distribution figure's linear panel gives way to a logarithmic one, as a
#: percentile of the distances plotted. The published figure split at a hand-picked 7, which
#: is the 83rd percentile of both published tables (7.6 and 6.6). Carrying the percentile
#: rather than the 7 at least rescales with the data; it is still one number read off one
#: experiment.
SPLIT_PERCENTILE = 83.0

#: Inches per fit panel, and how many leaves an inch of dendrogram carries. Both as published.
PANEL_SIZE = 8.0
LEAVES_PER_INCH = 7.3


def parse_parameters(text: str) -> list[float]:
    """`"b=1.2,d=3.4"` -> `[1.2, 3.4]`, the form `fit.tsv` stores fitted parameters in.

    Positional, because that is how the model functions take them; the names are there to be
    read by a person, and are checked against the model's own list by the caller.
    """
    if not text or not str(text).strip():
        return []
    values = []
    for item in str(text).split(","):
        _, _, value = item.partition("=")
        values.append(float(value))
    return values


# ------------------------------------------------------------------------- fit panels


def _panel(
    axes: Axes,
    feature: str,
    stratum: str,
    x,
    y,
    fits: pd.DataFrame,
    axis: Sequence[float],
    labels: Sequence[str],
    score: str = "aic_plus_bic",
) -> None:
    """One feature on one stratum: the distance points, the fitted curves, the legend.

    `score` is the column the legend is ordered by, which is `[fit] rank_by`'s column so
    that the panel lists the models in the order the retention rule compared them.
    """
    import numpy as np

    from cdr_fs.models import MODEL_FUNCTIONS, MODEL_PARAMETER_NAMES

    axes.plot(x, y, "o", markersize=10)

    # Best first by the configured ranking column, so the order a reader sees is the order
    # the retention rule read them in.
    entries = []
    for row in fits.sort_values(score).itertuples():
        function = MODEL_FUNCTIONS.get(row.model)
        parameters = parse_parameters(row.parameters)
        if function is None or len(parameters) != len(MODEL_PARAMETER_NAMES[row.model]):
            continue
        colour = MODEL_COLOURS.get(row.model, "gray")
        # A fitted sigmoid can overflow when evaluated: `log(x + 1e-10)` at x = 0 is -23, and
        # a steep b sends the exponential past the float range. The point simply does not
        # draw, which is the right outcome, so the warning is noise.
        with np.errstate(over="ignore", invalid="ignore"):
            curve = function(x, *parameters)
        axes.plot(x, curve, "-", color=colour, alpha=0.9, lw=3)
        label = f"{row.model} AIC/BIC: {row.aic:.2f}/{row.bic:.2f}"
        if row.model == "Lin":
            label += f", Slope: {row.slope:.2f}"
        entries.append((label, colour))

    title = f"{feature}, {stratum}" if stratum else feature
    axes.set_title(title, fontsize=16, wrap=True)
    axes.set_ylabel("EMD score", fontsize=16)
    # Ticks come from the configured axis, not from the points present, so that a series
    # with a hole shows the hole instead of closing up.
    axes.set_xticks(list(axis))
    axes.set_xticklabels(labels, rotation=45, fontsize=6)
    axes.tick_params(axis="y", labelsize=16)

    if entries:
        from matplotlib.lines import Line2D

        axes.legend(
            [Line2D([0], [0], color=colour, lw=3, alpha=0.9) for _, colour in entries],
            [label for label, _ in entries],
            fontsize=16,
        )
    else:
        # A series with a hole is not fitted at all, so the panel would otherwise be points
        # with no explanation. See `fit.py` on why an incomplete series is skipped.
        axes.text(
            0.5,
            0.95,
            "not fitted: incomplete exposure series",
            transform=axes.transAxes,
            ha="center",
            va="top",
            fontsize=14,
            color="#b2182b",
        )


def plot_fit_panels(
    config: Config,
    emd_table: pd.DataFrame,
    fit_table: pd.DataFrame,
    destination: str | Path,
    *,
    grid: int = 3,
    features: Iterable[str] | None = None,
    dpi: int = 100,
) -> list[Path]:
    """Draw every fitted series as a panel, `grid` x `grid` panels per page.

    One page series per stratum, named `fit_<stratum>_part_<n>.png`. `features` restricts and
    orders the panels; the default is every feature the distance table holds, in its order.
    Returns the paths written.
    """
    from matplotlib.figure import Figure

    from cdr_fs.fit import axis_positions, series_from_emd
    from cdr_fs.select import SCORE_COLUMNS

    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    if grid < 1:
        raise ValueError(f"grid must be at least 1, got {grid}")

    # An empty list is not the same as no list: `None` means "every feature in the
    # distance table", where `[]` means a list was supplied and held nothing. Drawing
    # zero panels for it is the quietest possible failure - the caller asked for the
    # figure, both tables were there, and nothing came back. Refuse instead, which the
    # CLI reports as a skipped figure with this reason.
    if features is not None and not list(features):
        raise ValueError("the feature list is empty, so there is nothing to draw")
    wanted = None if features is None else list(dict.fromkeys(features))
    order = (
        {name: index for index, name in enumerate(wanted)}
        if wanted is not None
        else {
            name: index
            for index, name in enumerate(dict.fromkeys(emd_table["feature"]))
        }
    )
    positions = axis_positions(config)
    axis = [positions[level] for level in config.design.fitted_levels]
    labels = [
        f"{config.design.control}v{level}" for level in config.design.fitted_levels
    ]
    indexed = fit_table.set_index(["feature", "stratum"]).sort_index()
    score = SCORE_COLUMNS[config.fit.rank_by]

    # Collected per stratum first, so the pages of one stratum are consecutive and complete.
    by_stratum: dict[str, list[tuple]] = {}
    for feature, stratum, x, y, _ in series_from_emd(config, emd_table):
        if feature in order:
            by_stratum.setdefault(stratum, []).append((order[feature], feature, x, y))

    per_page = grid * grid
    written = []
    for stratum, panels in by_stratum.items():
        panels.sort(key=lambda entry: entry[0])
        for page, start in enumerate(range(0, len(panels), per_page), start=1):
            figure = Figure(figsize=(PANEL_SIZE * grid, PANEL_SIZE * grid))
            axes_grid = figure.subplots(grid, grid, squeeze=False)
            for offset, (_, feature, x, y) in enumerate(panels[start : start + per_page]):
                axes = axes_grid[offset // grid][offset % grid]
                try:
                    fits = indexed.loc[[(feature, stratum)]].reset_index()
                except KeyError:
                    fits = fit_table.iloc[0:0]
                _panel(axes, feature, stratum, x, y, fits, axis, labels, score)
            # Leave the unused cells of a partly filled page empty rather than framed.
            for offset in range(len(panels[start : start + per_page]), per_page):
                axes_grid[offset // grid][offset % grid].set_axis_off()

            figure.tight_layout()
            name = f"fit_{stratum}_part_{page}.png" if stratum else f"fit_part_{page}.png"
            figure.savefig(destination / name, dpi=dpi)
            written.append(destination / name)
    return written


# ---------------------------------------------------------------- distance distribution


def plot_distribution(
    table: pd.DataFrame,
    path: str | Path,
    *,
    split: float | None = None,
    title: str | None = None,
    figsize: tuple[float, float] = (28.8, 15.0),
    dpi: int = 300,
    label_features: bool = True,
) -> Path:
    """Every distance in `table`, one column of points per feature, ordered by median.

    The y-axis is broken: a linear panel below `split` and a logarithmic one above it, because
    the distribution is heavy enough that one scale hides either the bulk or the tail. `split`
    defaults to the `SPLIT_PERCENTILE`th percentile of the finite distances.

    The same function draws either contrast set, because they are the same view of the same
    quantity: on `emd.tsv` the distances from the control to each exposure level, on
    `emd_baseline.tsv` the distances between control replicates - the reproducibility floor the
    first has to clear.
    """
    import numpy as np
    from matplotlib.figure import Figure
    from matplotlib.gridspec import GridSpec

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    absent = [column for column in ("feature", "emd") if column not in table.columns]
    if absent:
        raise ValueError(
            f"the distance table has no {' or '.join(absent)} column; it has "
            f"{', '.join(map(str, table.columns))}"
        )

    values = table["emd"].to_numpy(dtype=np.float64)
    finite = values[np.isfinite(values)]
    # An empty figure is worse than an error: it looks like a result. Refuse it, and say which
    # of the two ways the table was empty, since they have different causes - no rows at all
    # usually means the table was filtered down to nothing on the way in.
    if not len(table):
        raise ValueError("the distance table has no rows, so there is nothing to plot")
    if not finite.size:
        raise ValueError(
            f"none of the {len(table):,} distances in the table are finite, so there is "
            f"nothing to plot"
        )
    if split is None:
        split = float(np.percentile(finite, SPLIT_PERCENTILE))

    medians = table.groupby("feature")["emd"].median().sort_values()
    position = {feature: index for index, feature in enumerate(medians.index)}
    x = table["feature"].map(position).to_numpy(dtype=np.float64)

    figure = Figure(figsize=figsize)
    grid = GridSpec(2, 1, height_ratios=[1, 3], hspace=0.05, figure=figure)
    upper = figure.add_subplot(grid[0])
    lower = figure.add_subplot(grid[1])

    small = figsize[0] < 12
    size = 0.2 if small else 6
    alpha = 0.6 if small else 0.5
    scale = 0.8 if small else 1.5

    above = np.isfinite(values) & (values > split)
    below = np.isfinite(values) & (values <= split)
    upper.scatter(x[above], values[above], c="black", s=size, alpha=alpha)
    lower.scatter(x[below], values[below], c="black", s=size, alpha=alpha)

    upper.set_yscale("log")
    top = float(finite.max()) * 1.1 if finite.size else split * 2
    upper.set_ylim(split, max(top, split * 1.1))
    lower.set_ylim(0, split)
    for axes in (upper, lower):
        axes.set_xlim(-1, len(position))
        axes.spines["top"].set_visible(False)
        axes.spines["right"].set_visible(False)
    upper.set_xticks([])

    if label_features and not small:
        lower.set_xticks(range(len(position)))
        lower.set_xticklabels(list(medians.index), rotation=90, fontsize=4.2 * scale)
        lower.tick_params(axis="x", pad=1)
    else:
        lower.set_xticks([])
    for axes in (upper, lower):
        axes.tick_params(axis="y", labelsize=14 * scale)

    # supylabel rather than a hand-placed figure.text: the label belongs to both panels, and
    # placing it by hand puts it on top of the upper panel's tick labels at some figure sizes.
    figure.supylabel(
        f"EMD score  (linear below {split:.3g}, log above)", fontsize=13 * scale
    )
    figure.supxlabel("Features", fontsize=13 * scale)
    if title:
        figure.suptitle(title, fontsize=16 * scale)
    figure.subplots_adjust(left=0.08, top=0.94, bottom=0.15, right=0.99)
    figure.savefig(path, dpi=dpi, bbox_inches="tight")
    return path


# -------------------------------------------------------------------------- dendrogram


def _cluster_colours(count: int) -> list[str]:
    """`count` evenly spaced hues.

    The published figure used seaborn's `husl`. This is the same idea - walk the hue circle -
    without the dependency; the hues are not perceptually equalised, which changes nothing
    about which leaves share a colour.
    """
    from matplotlib import colormaps
    from matplotlib.colors import to_hex

    if count < 1:
        return []
    wheel = colormaps["hsv"]
    return [to_hex(wheel(index / count)) for index in range(count)]


def plot_dendrogram(
    linkage_table: pd.DataFrame,
    path: str | Path,
    *,
    cut: float,
    clusters: pd.DataFrame | None = None,
    min_coloured: int = 3,
    dpi: int = 300,
    default_colour: str = "gray",
) -> Path:
    """Draw the tree the correlation stage cut, from `correlation_linkage.tsv`.

    Pass `correlation_clusters.tsv` as `clusters` to colour each cluster of at least `min_coloured`
    members - links and leaf labels alike, and consistently with each other. Both come from
    the same run, so the colours are the clustering rather than a redrawing of it: a link is
    coloured only when every leaf below it belongs to one cluster.
    """
    import numpy as np
    from matplotlib.figure import Figure
    from scipy.cluster.hierarchy import dendrogram

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    leaves = linkage_table[linkage_table["label"].notna()]
    merges = linkage_table[linkage_table["label"].isna()]
    names = leaves["label"].tolist()
    count = len(names)
    tree = merges[["left", "right", "height", "size"]].to_numpy(dtype=np.float64)
    if len(tree) != max(count - 1, 0):
        raise ValueError(
            f"{count} leaf/leaves and {len(tree)} merge(s): a linkage table needs n - 1 merges"
        )

    of_cluster: dict[str, int] = {}
    sizes: dict[int, int] = {}
    if clusters is not None:
        for row in clusters.itertuples():
            of_cluster[row.feature] = row.cluster
            sizes[row.cluster] = row.size
    big = [number for number, size in sorted(sizes.items()) if size >= min_coloured]
    palette = dict(zip(big, _cluster_colours(len(big))))
    leaf_colour = {
        name: palette.get(of_cluster.get(name, -1), default_colour) for name in names
    }

    # Every node's cluster, or None where its leaves disagree. Bottom-up, so each merge sees
    # its children already resolved.
    node_cluster: list[int | None] = [of_cluster.get(name, -1) for name in names]
    for left, right, _, _ in tree:
        first, second = node_cluster[int(left)], node_cluster[int(right)]
        node_cluster.append(first if first == second else None)
    link_colour = {
        node: palette.get(cluster, default_colour)
        for node, cluster in enumerate(node_cluster)
        if cluster is not None
    }

    width = max(10.0, count / LEAVES_PER_INCH)
    figure = Figure(figsize=(width, max(6.0, 0.6 * width)))
    axes = figure.add_subplot(111)
    dendrogram(
        tree,
        labels=names,
        ax=axes,
        leaf_rotation=90,
        leaf_font_size=10,
        link_color_func=lambda node: link_colour.get(node, default_colour),
    )
    for label in axes.get_xmajorticklabels():
        label.set_color(leaf_colour.get(label.get_text(), default_colour))

    axes.set_title("Feature clustering by correlation distance")
    axes.set_xlabel("Feature")
    axes.set_ylabel("Distance  (1 - |r|)")
    axes.axhline(y=cut, color="r", linestyle="--")
    figure.tight_layout()
    figure.savefig(path, dpi=dpi)
    return path
