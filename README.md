# cdr_FS — feature selection by concentration/dose–response model fitting

Selects morphological features in high-content screening by fitting
concentration/dose–response models to earth mover's distance scores between control and
treated cell populations.

> **Status: under construction.** Every stage runs except the diagnostic plots. See
> [Status](#status) for what works today and [Reproducing the published
> run](#reproducing-the-published-run) for which published numbers come back.

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
   cluster of near-redundant features, and drop any that are left too sparse to use.

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
| Contrast-driven EMD engine | `emd.py` | done |
| The six models, AIC/BIC | `models.py` | done (lifted verbatim) |
| Fitting to a results table | `fit.py` | done |
| Retention rule | `select.py` | done |
| Correlation pruning (optional) | `prune.py` | done |
| Applying a selected list to the data | `subset.py` | done |
| Diagnostic plots | `plots.py` | not started |

Working today:

```bash
cdr-fs check -c config.ini          # validate the configuration, report the schema
cdr-fs check -c config.ini --scan   # also confirm the design occurs in the data
cdr-fs trim  -c config.ini          # write the trimmed table  (--dry-run to just report)
cdr-fs emd   -c config.ini          # distances per feature, stratum and contrast
cdr-fs fit   -c config.ini          # fit the six models to each distance series
cdr-fs select -c config.ini         # apply the retention rule
cdr-fs prune -c config.ini          # collapse near-redundant features   (optional)
cdr-fs subset -c config.ini         # write the table restricted to what survived
```

Every stage that needs cell-level data reads `[input] table` and applies the configured
trim itself, so there is no requirement to materialise a trimmed copy — on the reference
dataset that copy would be another 3.7 GB. `cdr-fs trim` exists to write it out for
inspection, not because later stages need it.

Each stage reads the previous one's output from `[output] dir` under a stable name:

| Stage | Writes |
|---|---|
| `emd` | `emd.tsv` — control against each level; `emd_baseline.tsv` — control against control, between replicates |
| `fit` | `fit.tsv` |
| `select` | `selected.txt`; `select_evidence.tsv` — per feature and stratum, which model won, by how much, and what the linear slope was |
| `prune` | `pruned.txt`; `prune_clusters.tsv` — every feature with its cluster and how far it sits from the representative; `prune_linkage.tsv` — the tree and its leaf order, so the dendrogram can be drawn without recomputing a correlation |
| `subset` | `subset_<list>.tsv` — the table restricted to that list; `subset_<list>_features.tsv` — how much data each column holds and whether the missing-data filter dropped it; `subset_<list>_retained.txt` — the features that survived, one per line |

`cdr-fs run` exists in `--help` and refuses to run, naming what is missing. The nine original
scripts are kept unmodified in [`legacy/`](legacy/README.md) as the reference each rewritten
module is checked against; that directory disappears when the last of them has been
reproduced.

## Installation

```bash
pip install -e .
```

Python 3.10 or newer. Requires numpy, pandas, scipy and matplotlib; `pip install -e ".[dev]"`
adds pytest.

## Configuration

One `.ini` file describes a run, and `-c/--config` points at it from anywhere — no editing
paths inside scripts, and no requirement to run from a particular directory.
[`examples/published.ini`](examples/published.ini) is the reference configuration, with the
reasoning behind each value in its comments; start by copying it.

`.ini` carries no type information, so every value is parsed and cross-checked before any
data is read. A wrong configuration fails in the first second with a message naming the
section, the key and the fix, rather than producing plausible-looking wrong numbers. All
ten sections must be present even when every key in them is left at its default, so that
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
| | `dose` | *(empty)* | actual doses, index-matched to `levels`; needed only for `x_scale = dose` |
| | `exclude_from_fit` | *(empty)* | levels withheld from fitting but still measured |
| `trim` | `enabled` | *required* | trimming changes the data, so it must be stated |
| | `lower_percentile`, `upper_percentile` | *required when enabled* | the interval kept is closed |
| | `scope` | *required when enabled* | columns whose groups the percentiles are computed within |
| `emd` | `contrasts` | `control_vs_each` | or an explicit list, `11v10,11v9,...` |
| | `baseline` | `control_across_replicates` | the between-replicate reproducibility floor, or `none` |
| | `per_replicate` | `false` | see [Pooling replicates](#pooling-replicates) |
| `fit` | `models` | all six | any subset of `BC4,BC5,LL4,WB1.4,Lin,Con` |
| | `x_scale` | `rank` | `rank` (position in the series) or `dose` |
| | `rank_by` | `aic_plus_bic` | `aic_plus_bic`, `aic` or `bic` |
| `select` | `slope_positive` | `any` | positive linear slope on `any` or `all` strata |
| | `nonconstant` | `all` | constant model beaten on `any` or `all` strata |
| | `strata` | *(all present)* | which strata the rule is applied over |
| `prune` | `enabled` | *required* | |
| | `threshold` | `0.9` | on \|r\|; features are clustered by `1 - |r|` |
| | `linkage` | `average` | `average`, `complete` or `single` |
| | `representative` | `alphabetical` | which member of a cluster to keep — `alphabetical` or `first` |
| | `aggregate_by` | *(`[trim] scope`)* | columns identifying one experimental unit; correlations are computed between unit medians. Explicitly empty correlates the object rows as they are |
| | `fill_missing` | `column_mean` | `column_mean` or `none` — see [Collapsing redundant features](#collapsing-redundant-features) |
| `subset` | `drop_missing` | *required* | drop features too empty to carry; a filter on the output must be stated |
| | `max_missing` | `30` | percent; a feature missing this much of the table **or more** is dropped |
| `output` | `dir` | *required* | where results are written |

Three things this configuration makes explicit that were previously implicit in the
scripts:

- **The metadata list.** It lived in three hardcoded copies; now one pattern list drives
  every stage. Because the split is by inversion, a new measured channel needs no
  configuration change — but see the warning below.
- **The two halves of the retention rule.** They use *different* quantifiers: a positive
  slope on **any** stratum, non-constant on **all** strata. Separate keys make that
  visible, and let you choose `all`/`all` for the stricter rule.
- **The fit's x-axis.** See [Spacing the exposure axis](#spacing-the-exposure-axis). It
  deserves its own section, because the obvious reading of it is wrong.

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

### Spacing the exposure axis

`[fit] x_scale` decides where along the x-axis the exposure levels sit. `rank` puts them at
0, 1, 2, ... - evenly spread, which is what the published run did. `dose` puts them at the
values in `[design] dose`.

The catch is that the four sigmoid models evaluate `log(x)` themselves, so the axis they see
is not the axis you supplied. For the reference dilution series:

```
x_scale = rank  ->  log(x) = [-23.03, 0, 0.69, 1.10, 1.39, 1.61, 1.79, 1.95]
x_scale = dose  ->  log(x) = [  2.99, 3.55, 4.11, 4.67, 5.23, 5.79, 6.35, 6.91]
                             spacing 0.559616 throughout - exactly even
```

So **`dose` is what spaces the levels evenly in log-concentration**, and `rank` is not: it
flings the lowest exposure out to `log(1e-10)` while the other seven span two units. Neither
choice dominates, because the two halves of the retention rule want different things:

| the series really is... | `rank` | `dose` |
|---|---|---|
| a logistic in log-dose (textbook DRC) | LL4 AIC -41 | LL4 AIC **-298** |
| a straight line in log-dose | Lin AIC **-573** | Lin AIC -19 |

`rank` suits the **linear** model, because for a geometric series the rank *is* log-dose up to
an affine shift - and the linear slope is what the retention rule tests the sign of. `dose`
suits the **sigmoids**, which become proper log-logistic curves in concentration. The sign of
the linear slope is positive either way, so the slope test is unaffected by the choice; what
changes is how easily a sigmoid beats the constant model.

`rank` stays the default because it is what was published. Switch to `dose` when the exposure
series is not geometrically spaced - rank is then not a transform of log-dose at all, and the
sigmoid fits are being handed a meaningless axis. `examples/published.ini` already carries the
dose vector, so `x_scale = dose` works there with no further edits.

There is deliberately no `log_dose`: it would hand `log(dose)` to models that take a logarithm
themselves, and it is NaN for any dose below 1.

## Pooling replicates

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

## Collapsing redundant features

`cdr-fs prune` aggregates the objects to one median profile per experimental unit -
`[prune] aggregate_by`, one well of one replicate on one day in the published run - correlates
the features across those units, and clusters them on `1 - |r|` so that a strong *negative*
correlation counts as the redundancy it is. Average linkage, cut at `1 - [prune] threshold`.
One member of each cluster survives.

Two things about it are worth knowing before trusting the output.

**Missing values are filled by default, and it is not a small effect.** Trimming leaves gaps,
so `fill_missing = column_mean` substitutes the feature's overall mean for each missing object
value before aggregating - which is what the published run did. For a feature with many
trimmed values that pulls its unit medians toward the global mean and shrinks the spread the
correlation sees. On the reference dataset the choice moves the cluster count by eight, from
99 to 107. `none` medians whatever values are actually present.

**Average linkage cuts on cluster means, so a member can sit further from its representative
than the threshold suggests.** On the reference dataset the 83 dropped features sit a median
0.052 from the feature that stands for them — `|r|` = 0.95, as advertised — but the loosest
sits at 0.215, which is `|r|` = 0.785. That is what chaining does, and it is why
`prune_clusters.tsv` records the distance for every member rather than only the cluster
number: it is the column to look at before trusting one feature to speak for a group.

## Dropping features with too little data

Trimming removes values rather than rows, so a feature can be missing for most objects and
still be present as a column. `[subset] drop_missing` removes the ones too empty to be worth
carrying: a feature missing `[subset] max_missing` percent of the table **or more** is left out,
30% by default.

It is applied over the **whole** table, before anything downstream subsamples. That is one
decision for the whole experiment, so a feature is either in the analysis or out of it —
whereas filtering each subsample separately lets the same feature be present in one day's file
and absent from another, which makes the resulting sets incomparable.

On the reference dataset the default takes 99 features to 97, dropping two that are missing 47%
and 95% of their values. The next thinnest surviving feature is missing 13%, so anything from
14% to 46% gives the same 97 — the threshold sits on a wide plateau rather than on a knife
edge, which is the useful thing to know about choosing it.

Nothing else is filtered. A constant feature, or one whose surviving values all come from a
single object, is named in the report and flagged in `subset_<list>_features.tsv`; what to do
about that depends on the analysis, and guessing is worse than saying so.

## Reproducing the published run

`examples/published.ini` is the configuration of the run described in the article. Three
stages have a published output to check against, and all three return it:

| Stage | Published output | Result |
|---|---|---|
| `emd` | the two EMD tables — 16,946 treatment and 11,292 baseline distances | every one of the 28,238 population sizes exact, distances within 8.5e-13 relative |
| `select` | the two retained feature lists — **182** across all days, **374** on D5 alone | both reproduced as identical sets |
| `prune` | the all-days list after pruning — **99** features | reproduced, and its composition matches the published categorization across all 4 organelles x 7 measurement families |

The selection gate runs off the published fit table, so it isolates the retention rule from
the curve fitting; `tests/test_golden_selection.py` needs only a 1 MB committed fixture and
runs in a second. The other two need the large inputs and are opt-in — see [Tests](#tests).

Two notes on what those numbers are and are not.

**The pruning check is on composition, not just on the count.** A single wrong merge would
still give a plausible number, so the assertion is against the published feature
categorization, every cell of which is labelled in
[HCS-proc's figure](https://github.com/NIB-SI/HCS-proc). Reproducing the count *and* the
composition takes the same 99 features.

**The published run's final list is shorter than 99, and the last step of getting there moved.**
The pipeline dropped features with 30% or more missing values — plus two excluded by hand — but
it did so while building the subsets for its dimension-reduction analyses: *per day-subset*,
downstream of the selection. `cdr-fs subset` applies the same 30% rule over the whole table
instead, which takes 99 to **97**. Applying it once, before anything subsamples, is the change
worth making: a feature is then either in the analysis or out of it, where a per-subsample
decision lets the same feature be present in one day's file and absent from another. See
[Dropping features with too little data](#dropping-features-with-too-little-data). The two
hand-excluded features are not reproduced; a hardcoded exclusion list is not a rule.

### An aside that is really a check

The published run produced two selections, one per gate: across all days, and on D5 alone.
Neither list contains the other - `counts_RelateLysoCell` is retained across all days but not
on D5, and `counts_RelateMitoCell` the reverse. Only the hybrid retention rule can do that. A
positive slope is required on **any** stratum, so the all-days gate has four chances where D5
has one; if both tests used `all`, the D5 set would necessarily contain the all-days set. The
asymmetry is therefore visible in the published output itself, not merely in a reading of the
code.

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
inheriting it. One more file from the same record, `all_days_trimmed_features.txt` (1.4 GB),
is the published run's own intermediate and is what the pruning gate checks against.

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
oracle, on data with and without infinities; the EMD engine is exercised on designs that
deliberately are *not* the RTgill-W1 one; the retention rule and the pruning are checked on
hand-built cases where the right answer is known by construction; and `tests/cases.py` holds
the table of misconfigurations that must be rejected, one entry per validation rule.

Two golden regressions compare against the published run itself, so they need the large inputs
in `data/`. Each skips when its input is absent, and skips when it is present unless you opt
in, because reading gigabytes takes minutes:

```bash
CDR_FS_GOLDEN=1 pytest tests/test_golden.py tests/test_golden_prune.py -v
```

`test_golden.py` reproduces the published EMD tables from the untrimmed 3.7 GB input: 16,946
treatment distances and 11,292 baseline distances, every population size matching exactly and
every distance within 1e-9 relative. Population sizes are the sharper comparison — they are
integers, so no tolerance can hide a difference, and their matching is what shows the trim step
was reproduced value for value.

`test_golden_prune.py` starts instead from the published run's own selected-and-trimmed subset,
`all_days_trimmed_features.txt`, which isolates the pruning stage: its 182 feature columns are
the published selection and its values are already trimmed, so anything that fails there is
pruning and nothing else.

## Licence and citation

MIT — see [LICENSE](LICENSE). Copyright National Institute of Biology.

If you use this tool, please cite both it and the article in which the method was first
published; [CITATION.cff](CITATION.cff) has the machine-readable metadata for both.
