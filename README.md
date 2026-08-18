# cdr_FS — feature selection by concentration/dose–response model fitting

Selects morphological features in high-content screening by fitting
concentration/dose–response models to earth mover's distance scores between control and
treated cell populations.

> **Status: under construction.** The scaffold and the configuration, schema and trimming
> stages are in place; the EMD engine, the model fitting, the retention rule, the plots and
> the pruning are not yet. See [Status](#status) for what runs today. Until the golden
> regression passes, treat this repository as work in progress rather than a tool to
> depend on.

## The method

Given a tidy per-object table — one row per cell, columns split into metadata and measured
features — the selection asks of each feature a single question: **does its distributional
distance from the control grow with exposure in a way a concentration–response curve
describes better than a flat line?**

1. For every feature and every exposure level, compute the earth mover's (Wasserstein)
   distance between the control population and the exposed population. The distance
   between control populations of different biological replicates gives the
   reproducibility floor to read those numbers against.
2. Fit six models to the resulting distance-versus-exposure series: Brain–Cousens
   hormesis (BC4, BC5), four-parameter log-logistic (LL4), four-parameter Weibull
   (WB1.4), linear (Lin) and constant (Con).
3. Retain a feature when its linear slope is positive and the constant model is *not* the
   best fit by AIC + BIC. A feature whose distance from the control is flat, or shrinks,
   carries no concentration-dependent signal and is dropped.
4. Optionally collapse the survivors by correlation, keeping one representative per
   cluster of near-redundant features.

Using distributions rather than per-well means is the point: it detects a subpopulation
shifting while the population mean stays put, which is what subcytotoxic exposure tends to
look like.

### Why `cdr_FS`

**c**oncentration/**d**ose **r**esponse **F**eature **S**election.

The slash is deliberate. *Concentration–response* is the ecotoxicology term and
*dose–response* the pharmacological one for the same idea, and the method needs neither
field's assumptions — it applies to any ordered exposure series, including drug-discovery
screens. Machine-read fields (`pyproject.toml`, package metadata) use a plain hyphen and
prose uses the en dash; the package itself is plain lowercase `cdr_fs`.

## Status

| Stage | Module | State |
|---|---|---|
| Configuration loading and validation | `config.py` | done |
| Metadata/feature resolution, table reading | `schema.py` | done |
| Extreme-value trimming (optional) | `trim.py` | done |
| Contrast-driven EMD engine | `emd.py` | not started |
| The six models, AIC/BIC | `models.py` | not started |
| Fitting to a results table | `fit.py` | not started |
| Retention rule | `select.py` | not started |
| Diagnostic plots | `plots.py` | not started |
| Correlation pruning (optional) | `prune.py` | not started |
| Applying a selected list to the data | `subset.py` | not started |

Working today:

```bash
cdr-fs check -c config.ini          # validate the configuration, report the schema
cdr-fs check -c config.ini --scan   # also confirm the design occurs in the data
cdr-fs trim  -c config.ini          # write the trimmed table  (--dry-run to just report)
```

Every other subcommand exists in `--help` and refuses to run, naming the stage that is
missing. The nine original scripts are kept unmodified in [`legacy/`](legacy/README.md) as
the reference each rewritten module is checked against; that directory disappears when the
last of them has been reproduced.

## Installation

```bash
pip install -e .
```

Python 3.10 or newer. Requires numpy, pandas, scipy, matplotlib and scikit-learn; `pip
install -e ".[dev]"` adds pytest.

## Configuration

One `.ini` file describes a run, and `-c/--config` points at it from anywhere — no editing
paths inside scripts, and no requirement to run from a particular directory.
[`examples/published.ini`](examples/published.ini) is the reference configuration, with the
reasoning behind each value in its comments; start by copying it.

`.ini` carries no type information, so every value is parsed and cross-checked before any
data is read. A wrong configuration fails in the first second with a message naming the
section, the key and the fix, rather than producing plausible-looking wrong numbers. All
nine sections must be present even when every key in them is left at its default, so that
a mistyped section name is an error instead of a silent fall-back; unknown keys are
rejected for the same reason.

| Section | Key | Default | Meaning |
|---|---|---|---|
| `input` | `table` | *required* | path to the per-object table |
| | `sep` | `tab` | `tab`, `comma` or `semicolon` — a keyword, because `.ini` has no escapes, so `\t` would arrive as a literal backslash-t |
| `schema` | `metadata_patterns` | *required* | one regex per line; **everything not matched is a feature** |
| | `condition` | *required* | column holding the exposure level |
| | `group_by` | *(empty)* | stratification, e.g. time point; empty means one stratum |
| | `pool_over` | *(empty)* | replicate column merged into one distribution |
| `design` | `control` | *required* | the control level's label |
| | `levels` | *required* | exposure levels, ordered **low to high** — this is the response axis |
| | `dose` | *(empty)* | actual doses, index-matched to `levels`; needed only for `x_scale = dose`/`log_dose` |
| | `exclude_from_fit` | *(empty)* | levels withheld from fitting but still measured |
| `trim` | `enabled` | *required* | trimming changes the data, so it must be stated |
| | `lower_percentile`, `upper_percentile` | *required when enabled* | the interval kept is closed |
| | `scope` | *required when enabled* | columns whose groups the percentiles are computed within |
| `emd` | `contrasts` | `control_vs_each` | or an explicit list, `11v10,11v9,...` |
| | `baseline` | `control_across_replicates` | the between-replicate reproducibility floor, or `none` |
| | `per_replicate` | `false` | see [Pooling replicates](#pooling-replicates) |
| `fit` | `models` | all six | any subset of `BC4,BC5,LL4,WB1.4,Lin,Con` |
| | `x_scale` | `rank` | `rank`, `dose` or `log_dose` |
| | `rank_by` | `aic_plus_bic` | `aic_plus_bic`, `aic` or `bic` |
| `select` | `slope_positive` | `any` | positive linear slope on `any` or `all` strata |
| | `nonconstant` | `all` | constant model beaten on `any` or `all` strata |
| | `strata` | *(all present)* | which strata the rule is applied over |
| `prune` | `enabled` | *required* | |
| | `threshold` | `0.9` | on \|r\|; features are clustered by `1 - |r|` |
| | `linkage` | `average` | `average`, `complete` or `single` |
| | `representative` | `alphabetical` | which member of a cluster to keep |
| `output` | `dir` | *required* | where results are written |

Three things this configuration makes explicit that were previously implicit in the
scripts:

- **The metadata list.** It lived in three hardcoded copies; now one pattern list drives
  every stage. Because the split is by inversion, a new measured channel needs no
  configuration change — but see the warning below.
- **The two halves of the retention rule.** They use *different* quantifiers: a positive
  slope on **any** stratum, non-constant on **all** strata. Separate keys make that
  visible, and let you choose `all`/`all` for the stricter rule.
- **The fit's x-axis.** `rank` fits against the position of a level in the series.
  For a geometric dose series that is an exact affine transform of log-dose, so it is
  equivalent to `log_dose`; for an irregular series it is not, and `log_dose` is what you
  want.

### One warning worth reading before you write `metadata_patterns`

Patterns name the *metadata*, so anything you fail to match becomes a feature — and
anything you match too broadly stops being one. In the reference dataset five columns begin
with `counts_`, and only three of them are metadata:

```
counts_Cells, counts_Cytoplasm, counts_FilteredNuclei     segmentation QC  -> metadata
counts_RelateLysoCell, counts_RelateMitoCell              organelle counts -> FEATURES
```

The last two are the per-cell lysosomal and mitochondrial counts the article appends to
each cell's profile, and `counts_RelateLysoCell` survives the entire selection into the
published retained list. A tidy-looking `^counts_` pattern would classify all five as
metadata and silently discard a feature the published run kept.

`cdr-fs check` therefore prints the resolved split grouped by leading name token, so the
minority groups are visible:

```
[columns]  481 = 10 metadata + 471 feature(s)
  features          rp_* 469, counts_* 2
    counts_*        counts_RelateLysoCell, counts_RelateMitoCell
```

It also warns about any pattern that matched nothing, since that is how a metadata column
quietly turns into a feature.

### Pooling replicates

When computing the distance between control and exposure level *C*, what counts as "the
population"? `per_replicate = false`, the default, merges all replicates at that level into
one distribution: one distance per feature, stratum and contrast, so the curve is fitted to
as many points as there are exposure levels. This is the published behaviour, and it is
defensible here because the upstream two-step row/plate standardization exists precisely to
remove the batch structure that pooling would otherwise bake in.

`per_replicate = true` instead computes one distance per replicate. That puts batch
variation *between* the points rather than inside the distributions, gives each point a
spread, and multiplies the number of points the curve is fitted to — which matters, because
a five-parameter model fitted to eight points is at the edge of what AIC/BIC comparison
supports. It also changes the numbers, so it no longer reproduces the article. Both are
available; the default reproduces the paper.

## Scope

**In:** the schema declaration, optional trimming, contrast definition, the EMD
computation, multi-model fitting with information-criterion ranking, the retention rule
with explicit quantifiers, per-feature diagnostic panels, optional correlation pruning, and
a manifest recording the configuration and input hashes.

**Out:** image quality control, CellProfiler, segmentation, per-object pooling, and
row/plate standardization — all of which are tied to a particular plate design and belong
upstream. Also out: UMAP, MMD and Mahalanobis analyses, which consume this tool's output.
[HCS-proc](https://github.com/NIB-SI/HCS-proc) covers all of those.

## Origin

This is an extraction of `scripts/feature_selection/` from
[NIB-SI/HCS-proc](https://github.com/NIB-SI/HCS-proc), the pipeline published with:

> Tome, M.; Jozef, B.; Mosimann, S. L.; Kosnik, M.; Schirmer, K.; Županič, A.
> *A High-Content Imaging Pipeline to Investigate Subcytotoxic Effects in RTgill-W1 Cells.*
> **Environmental Science & Technology** **2026**, *60* (31), 21402–21416.
> <https://doi.org/10.1021/acs.est.5c18316>

The git history of those scripts is preserved here. HCS-proc remains the citable record of
the published pipeline; this repository generalises one stage of it into a tool that is
configuration-driven and installable, and that does not assume the RTgill-W1 experimental
design.

The design is experiment-agnostic by construction — nothing about days, replicates, wells
or nine concentrations is hardcoded — but it has so far been exercised on **one**
experimental design. Structural generality is tested by re-shaping that dataset (dropping
the time axis, changing the number of levels, renaming columns); a second real dataset has
not been run.

## Data

The reference dataset is 121 GB in total and lives on Zenodo:
<https://doi.org/10.5281/zenodo.17951792> (CC-BY-4.0). The single file this tool starts
from is `cell_ID_pooled_median_row_plate_standardization_cid.txt`, 3.7 GB — the untrimmed,
row/plate standardized per-cell table, 503,920 cells x 481 columns. It is deliberately
*untrimmed*: trimming lives in this tool, so a run reproduces that step rather than
inheriting it.

Nothing that large belongs in git. Put it in `data/`, which is ignored; the committed
fixtures under [`tests/fixtures/`](tests/fixtures/README.md) are small subsets that let the
suite run without it.

## Tests

```bash
pytest
```

The suite needs no large data. `tests/fixtures/columns_published.txt` carries the real 481
column names, so the metadata/feature arithmetic is asserted against the published schema
in CI; trimming is checked against the original per-group `np.nanpercentile` loop as its
oracle; and `tests/cases.py` holds the table of misconfigurations that must be rejected,
one entry per validation rule.

## Licence and citation

MIT — see [LICENSE](LICENSE). Copyright National Institute of Biology.

If you use this tool, please cite both it and the article in which the method was first
published; [CITATION.cff](CITATION.cff) has the machine-readable metadata for both.
