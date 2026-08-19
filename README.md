# cdr_FS — feature selection by concentration/dose–response model fitting

A method, and a working implementation of it. For each morphological feature, measure how
far its distribution moves away from the control as exposure rises; then keep the features
where a concentration–response curve describes that movement better than a flat line does.
The distance is the earth mover's (Wasserstein) distance between cell populations, and the
curves are six standard concentration–response models ranked by AIC and BIC.

> **Status: every stage runs.** The method was developed for, and has been applied to, one
> dataset — the RTgill-W1 screen of [Tome et al. 2026](#origin). The experimental design is
> read from a configuration file rather than written into the code, and that is tested; it
> is a statement about the code, not evidence that the method carries to another organism,
> assay or chemistry. [Reproducing the published run](#reproducing-the-published-run) says
> which published numbers come back and how. Still to do: a tagged release.

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

Using distributions rather than per-well means is the point: a subpopulation can shift
while the population mean stays put, and a distance between distributions registers that
where a difference of means does not. That is the effect the article set out to find at
subcytotoxic exposures, and it is why the method is built on distances rather than on
summaries.

**The whole thing rests on one assumption**, and it is worth being blunt about: that a
higher exposure pushes a feature's distribution further from the control than a lower one
does. That is what makes the distances a series worth fitting a curve to, and it is what the
positive-slope gate tests. A response that is not monotone in exposure — one that saturates,
reverses, or only appears in a middle band — is not what this looks for.

### Why `cdr_FS`

**c**oncentration/**d**ose **r**esponse **F**eature **S**election.

The slash is deliberate. *Concentration–response* is the ecotoxicology term and
*dose–response* the pharmacological one for the same idea, and the name commits to neither.
What the code asks of an experiment is an ordered exposure series and a control; whether
that series is a dilution in a fish-gill assay or something else entirely is not a thing it
inspects, and not a thing it has been tried on. Machine-read fields (`pyproject.toml`,
package metadata) use a plain hyphen and prose uses the en dash; the package itself is plain
lowercase `cdr_fs`.

## Status

| Stage | Module | State |
|---|---|---|
| Configuration loading and validation | `config.py` | done |
| Metadata/feature resolution, table reading | `schema.py` | done |
| Extreme-value trimming (optional) | `trim.py` | done |
| Contrast-driven EMD engine | `emd.py` | done |
| The six models, AIC/BIC | `models.py` | done (taken unchanged) |
| Fitting to a results table | `fit.py` | done |
| Retention rule | `select.py` | done |
| Correlation pruning (optional) | `prune.py` | done |
| Applying a selected list to the data | `subset.py` | done |
| Diagnostic figures | `plots.py` | done |

The stages, in order:

```bash
cdr-fs check -c config.ini          # validate the configuration, report the schema
cdr-fs check -c config.ini --scan   # also confirm the design occurs in the data
cdr-fs trim  -c config.ini          # write the trimmed table  (--dry-run to just report)
cdr-fs emd   -c config.ini          # distances per feature, stratum and contrast
cdr-fs fit   -c config.ini          # fit the six models to each distance series
cdr-fs select -c config.ini         # apply the retention rule
cdr-fs prune -c config.ini          # collapse near-redundant features   (optional)
cdr-fs subset -c config.ini         # write the table restricted to what survived
cdr-fs plot  -c config.ini          # draw the figures from whatever tables exist
```

Every stage that needs cell-level data reads `[input] table` and applies the configured
trim itself, so there is no requirement to materialise a trimmed copy — on the reference
dataset that copy would be another 3.9 GB. `cdr-fs trim` exists to write it out for
inspection, not because later stages need it.

Each stage reads the previous one's output from `[output] dir` under a stable name:

| Stage | Writes |
|---|---|
| `emd` | `emd.tsv` — control against each level; `emd_baseline.tsv` — control against control, between replicates |
| `fit` | `fit.tsv` |
| `select` | `selected.txt`; `select_evidence.tsv` — per feature and stratum, which model won, by how much, and what the linear slope was |
| `prune` | `pruned.txt`; `prune_clusters.tsv` — every feature with its cluster and how far it sits from the representative; `prune_linkage.tsv` — the tree and its leaf order, so the dendrogram can be drawn without recomputing a correlation |
| `plot` | `fit_<stratum>_part_<n>.png` — the fit panels; `emd.png` and `emd_baseline.png` — the distance distributions; `dendrogram.png` — the tree pruning cut |
| `subset` | `subset_<list>.tsv` — the table restricted to that list; `subset_<list>_features.tsv` — how much data each column holds and which rule, if any, removed it; `subset_<list>_retained.txt` — the features that survived, one per line |

A stage refuses rather than write an artefact that would read as a result: `fit` will not
write a table when no series was complete, and `prune` and `subset` will not run on an empty
feature list. Each exits 3 and says what came up empty.

The file names above use `.tsv`; with `[input] sep = comma` they are `.csv`, since the stage
tables are written with the same separator the input uses. Read them back with
`cdr_fs.schema.read_stage_table` rather than with a bare `pandas.read_csv`: an experiment
with no `[schema] group_by` has one stratum whose label is the empty string, and pandas reads
an empty field as `NaN` — which `groupby` then drops, taking the whole table with it.

There is no single command that chains the stages. Run them in the order above; each
reports what it did and refuses when its input is not there yet.

The nine scripts this package was extracted from were carried along in a `legacy/` directory
while each one was reproduced and checked. All nine now are, so that directory is gone; the
scripts remain in this repository's history, and in
[HCS-proc](https://github.com/NIB-SI/HCS-proc)'s `scripts/feature_selection/` as the published
record. `models.py` was extracted from the script holding the model functions in the commit that
removed it; git does not score that as a rename, so the original author's history is reached
through the old path — `git log --oneline -- legacy/plots_emd_model_drc.py`.

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
section, the key and — where there is one to name — the fix, rather than producing
plausible-looking wrong numbers. All
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
| | `exclude` | *(empty)* | exact feature names to leave out of the final table, whatever else says |
| `output` | `dir` | *required* | where results are written |

Three things this configuration makes explicit that were previously implicit in the
scripts:

- **The metadata list.** It lived as a literal in five scripts, in two versions; now one drives
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
[columns]  481 = 18 metadata + 463 feature(s)
  features          rp_* 461, counts_* 2
    counts_*        counts_RelateLysoCell, counts_RelateMitoCell
```

It also warns about any pattern that matched nothing, since that is how a metadata column
quietly turns into a feature.

The mistake runs the other way too, and it is harder to see. A column can be named like a
measurement and still be a label:

```
rp_norm_Number_Object_Number_*      CellProfiler's object index -> metadata
```

`Number_Object_Number` is the running number of an object within its image. It has no
biological content, but it *does* respond to exposure — the number of objects per image changes
with dose, which shifts the index distribution, which moves the earth mover's distance. Seven of
its eight columns passed the concentration-response gate in the published run and two survived
correlation pruning; the pipeline removed them by name at the very end. Declaring them metadata
is the same decision made where it belongs, and `examples/published.ini` does. Nothing in the
tool guesses this for you: only a pattern you write keeps a label out of the analysis.

### Spacing the exposure axis

`[fit] x_scale` decides where along the x-axis the exposure levels sit. `rank` puts them at
0, 1, 2, ... - evenly spread, which is what the published run did. `dose` puts them at the
values in `[design] dose`.

The catch is that the four sigmoid models evaluate `log(x)` themselves, so the axis they see
is not the axis you supplied. For the reference dilution series:

```
x_scale = rank  ->  log(x) = [-23.03, 0, 0.69, 1.10, 1.39, 1.61, 1.79, 1.95]
x_scale = dose  ->  log(x) = [  2.43, 2.99, 3.55, 4.11, 4.67, 5.23, 5.79, 6.35]
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
suits the **sigmoids**, which become proper log-logistic curves in concentration. What changes
most is how easily a sigmoid beats the constant model.

The slope test is not strictly invariant either, though it is close enough to be. A
least-squares slope is `cov(x, y)/var(x)`, so re-spacing x changes it — and `Lin` takes no
logarithm, which means `x_scale = dose` hands it raw concentration, exponential in rank rather
than affine to it. On arbitrary eight-point series the two axes then disagree about the sign
one time in seven. For a series that rises monotonically with exposure, which is the case the
rule exists to find, the sign is positive on either axis.

`rank` stays the default because it is what was published. Switch to `dose` when the exposure
series is not geometrically spaced - rank is then not a transform of log-dose at all, and the
sigmoid fits are being handed a meaningless axis. `examples/published.ini` already carries the
dose vector, so `x_scale = dose` works there with no further edits.

There is deliberately no `log_dose`: it would hand `log(dose)` to models that take a logarithm
themselves, and it is NaN for any dose below 1.

## Describing your experiment

Of the thirty-odd keys across the ten sections, twelve describe the **experiment**. The
rest describe the **method** — thresholds, which models to fit, linkage, the two
quantifiers — and can be left alone while you get the first twelve right.

| Key | The fact about the experiment it states |
|---|---|
| `[schema] metadata_patterns` | where the boundary between bookkeeping and measurement lies |
| `[schema] condition` | which column carries each object's exposure level |
| `[schema] group_by` | the one stratification factor — day, cell line, plate. Empty means a single stratum |
| `[schema] pool_over` | which column is the replicate axis. Empty means no replicate structure |
| `[design] control` | which label is the unexposed arm. Exactly one |
| `[design] levels` | the ordered exposure axis, **low to high**. This order *is* the response axis |
| `[design] dose` | the actual exposure magnitudes, index-matched to `levels` — the spacing, not just the order |
| `[design] exclude_from_fit` | levels measured but outside the regime the curve should describe |
| `[trim] scope` | the columns identifying one physical unit of the assay |
| `[emd] contrasts` | which comparisons the design supports |
| `[select] strata` | which strata the run covers — see the warning below |
| `[prune] aggregate_by` | the columns identifying one experimental unit for correlation |

`[fit] models` sits between the two: it is a method choice, but the design constrains it,
because a curve needs more points than it has parameters.

### When a piece of the design is missing

None of it is required to exist. What changes is what you can ask.

| You have… | Write | What you lose |
|---|---|---|
| **no time course** — one batch, one readout | `group_by =` (empty) | Nothing. There is one stratum, and the `any`/`all` quantifiers collapse to the same rule. `[select] strata` must then be empty too, and the tool says so |
| **a stratification that is not time** — cell line, plate, passage | `group_by = <that column>` | Nothing; no part of the code is time-aware. But strata are never *compared*: one curve is fitted per stratum and the quantifiers combine the verdicts |
| **no replicates at all** | omit `pool_over`, **and** set `[emd] baseline = none` | The between-replicate reproducibility floor, which is the yardstick the treatment distances are read against. `[emd] per_replicate` also goes. Leaving `baseline` at its default here is an error, and the message names the key to change |
| **two replicates** | `pool_over = <column>` | Nothing structural — but the floor is then one number per feature and stratum rather than a spread |
| **five exposure levels** | `[fit] models = BC4,LL4,WB1.4,Lin,Con` | BC5. Five points cannot identify five parameters |
| **three or four levels** | `[fit] models = Lin,Con` | Every sigmoid. You are asking whether a line beats a flat line — still a question, but not concentration–response modelling |
| **two levels** | — | This design cannot be run: `[select]` requires both `Lin` and `Con`, and `Lin` needs three points |
| **no need to trim** | `[trim] enabled = false` | Nothing, except that `[prune] aggregate_by` no longer has `[trim] scope` to default to, so set it explicitly |
| **no wish to prune** | `[prune] enabled = false` | The redundancy collapse. `cdr-fs subset` then applies `selected.txt` instead of `pruned.txt`, and says which it used |

Two of these deserve spelling out.

**Level labels are opaque strings, and nothing sorts them.** The order in `[design] levels`
is the exposure axis, whatever the labels look like. The reference dataset's own labels run
`11, 10, … 2` with `11` the control and `2` the top dose, so its `levels` line reads
`10,9,8,…,2` — descending as text, ascending in exposure. That is also the quietest way to
get a wrong answer: a `levels` list written high to low runs to completion, produces an
identical distance table, negates every linear slope, and retains nothing. The only
automatic guard is `[design] dose` — when you supply it, the doses must rise along `levels`,
and a contradiction is refused. **Supply the doses even when you fit on `rank`.** It costs
one line and it is the only check that can catch a reversed axis.

**`[prune] aggregate_by` must name a unit that varies along the exposure axis.**
Correlations are computed between unit medians, so if each unit spans the whole dilution
series then every unit looks alike and pruning collapses nothing. In the published plate
layout one well holds one concentration, which is why `Metadata_Day,Metadata_Biorep,
Metadata_Well` works there; if your wells each hold a whole series, add the level column.
The report line `|r| >= 0.9 on median profiles over N unit(s) of …` is where you check that
N is what you expected.

### The smallest configuration that runs

Every section must be present. Nearly every key can be left out — except that `[emd]
baseline` defaults to comparing replicates with each other, so an experiment without them
has to say so.

```ini
[input]
table = /PATH/TO/objects.tsv
[schema]
metadata_patterns = ^level$
condition = level
[design]
control = C0
levels = L1,L2,L3,L4,L5,L6
[trim]
enabled = false
[emd]
baseline = none
[fit]
[select]
[prune]
enabled = false
[subset]
drop_missing = false
[output]
dir = /PATH/TO/results
```

Six exposure levels, because that is what all six models need; one metadata column, one
control, everything else defaulted. `cdr-fs check` lists the twenty-odd keys it filled in.

### What the configuration cannot say

Some of these fail loudly. The ones that do not are the ones worth knowing.

**Refused, with a message.** More than one control: `[design] control` is a single label
and every contrast is control-versus-level, so a vehicle-versus-untreated comparison is
rejected — as is any other contrast that is not control-versus-level. A two-arm design, as
above. A level label containing a comma, since `levels` is comma-separated.

**Accepted, and quietly not what you meant.** What the configuration describes is one
exposure axis and at most one stratification factor. Anything richer — a chemical A ×
chemical B factorial, day *and* cell line, plates nested inside batches inside days — has to
be flattened into a single column upstream, and if you flatten it wrongly the run still
completes. There is no pairing either: the distance is between two independent empirical
distributions, so a repeated-measures design runs with its pairing discarded.

One of these deserves stating on its own, because it is a property of the statistic and not
of the format. **The earth mover's distance is unsigned.** A feature whose values *fall* with
exposure is retained exactly like one whose values rise — "positive slope" describes the
distance from the control growing, not the measurement growing. Which direction a retained
feature moved in is a question for the data, and this tool does not answer it.

One more, because it surprises people: **`[select] strata` narrows the whole run, not only
the selection.** `emd` and `fit` read it too, so setting it to one stratum means the distance
table holds only that stratum, with nothing in the file to say it is partial.

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

## The figures

`cdr-fs plot` draws three things, from the tables the other stages wrote and never by
recomputing anything:

- **Fit panels** — one panel per feature and stratum, with the distance points, the six fitted
  curves and a legend ordered by information criterion. This is the figure the article's
  Figure 4 was composed from. `--grid N` sets panels per row and column, `--features FILE`
  restricts which features are drawn; without it a full run is several hundred pages.
- **Distance distributions** — every distance in a table, one column of points per feature,
  features ordered by median, on a broken axis (linear below a split, logarithmic above). On
  `emd_baseline.tsv` this is the between-replicate reproducibility floor; on `emd.tsv` it is
  the treatment distances to read against that floor.
- **Dendrogram** — the tree correlation pruning cut, with the cut drawn on it and each cluster
  of three or more members coloured.

A curve in a fit panel is the fitted function evaluated at the exposure levels and joined up,
as the published figures show it. The `fit.tsv` parameters are written at full precision so a
drawn curve reproduces the AIC printed beside it — which is asserted, and which fails for about
one fit in sixteen if the parameters are rounded to six digits, because a fitted inflection can
settle right against the `log(x + 1e-10)` guard.

The dendrogram colours links and leaf labels from the same cluster table, so a coloured subtree
is a cluster. The original drew the tree with one linkage method and cut it with another, so
its colours were only approximately the clustering.

## Dropping features from the final table

Trimming removes values rather than rows, so a feature can be missing for most objects and
still be present as a column. `[subset] drop_missing` removes the ones too empty to be worth
carrying: a feature missing `[subset] max_missing` percent of the table **or more** is left out,
30% by default.

It is applied over the **whole** table, before anything downstream subsamples. That is one
decision for the whole experiment, so a feature is either in the analysis or out of it —
whereas filtering each subsample separately lets the same feature be present in one day's file
and absent from another, which makes the resulting sets incomparable.

On the reference dataset the default takes 99 features to 97, dropping two that are missing 47%
and 95% of their values. The next thinnest surviving feature is missing 13%, and the rule
drops at the threshold or above, so every threshold in (13%, 47%] gives the same 97 — it sits
on a wide plateau rather than on a knife edge, which is the useful thing to know about
choosing one.

`[subset] exclude` is the other half: exact feature names to leave out whatever their quality —
the escape hatch for a judgement no rule expresses, such as a feature known to be an artifact of
one assay. Exact names rather than patterns, because a regex that quietly takes a second feature
with it is the wrong tool for a decision made one feature at a time; a name that matches nothing
is reported by `cdr-fs check`, since the failure mode is that the feature stays in and nobody
notices. `subset_<list>_features.tsv` records which rule removed each feature, and an explicit
exclusion is reported as an exclusion even when the missing-data rule would also have caught it.

Nothing else is filtered. A constant feature, or one whose surviving values all come from a
single object, is named in the report and flagged in `subset_<list>_features.tsv`; what to do
about that depends on the analysis, and guessing is worse than saying so.

## Reproducing the published run

`examples/published.ini` is the configuration for the dataset the method was published on.
Four published outputs are checked against, each stage given the published input to the stage
before it rather than this tool's own output, so that a difference is attributable to the
stage under test:

| Stage | Published output | Result |
|---|---|---|
| `emd` | the two EMD tables — 16,946 treatment and 11,292 baseline distances | both population sizes on every one of the 28,238 rows exact, distances within 8.5e-13 relative |
| `select` | the two retained feature lists — **182** across all days, **374** on D5 alone | both reproduced as identical sets |
| `prune` | the all-days list after pruning — **99** features | reproduced, and its composition matches the published categorization across all 4 organelles x 7 measurement families |
| `subset` | the final retained list — **95** features | reproduced: 94 `rp_norm_*` plus `counts_RelateLysoCell` |

The selection gate runs off the published fit table, so it isolates the retention rule from
the curve fitting; `tests/test_golden_selection.py` needs only a 1 MB committed fixture and
runs in a second. The others need the large inputs and are opt-in — see [Tests](#tests).

Three notes on what those numbers are and are not.

**The pruning check is on composition, not just on the count.** A single wrong merge would
still give a plausible number, so the assertion is against the published feature
categorization, every cell of which is labelled in
[HCS-proc's figure](https://github.com/NIB-SI/HCS-proc). Reproducing the count *and* the
composition takes the same 99 features.

**The 182 and the 99 belong to the published metadata split, and the 95 does not.** The
published run treated the eight `Number_Object_Number` columns as measured features; they are
CellProfiler's within-image object index, so `examples/published.ini` declares them metadata —
see the warning [above](#one-warning-worth-reading-before-you-write-metadata_patterns). The two
gates above are asserted against the published split, because reproducing a published
intermediate means reproducing the choices that produced it. With the object indices out, the
same rules give **463 features → 175 after selection → 97 after pruning → 95**.

**The last step of the published route moved, and it lands in the same place.** The pipeline
dropped features with 30% or more missing values while building the subsets for its
dimension-reduction analyses: *per day-subset*, downstream of the selection, along with two
features excluded by hand. `cdr-fs subset` applies the same 30% rule once, over the whole table,
before anything subsamples — so a feature is either in the analysis or out of it, where a
per-subsample decision lets the same feature be present in one day's file and absent from
another. That takes 97 to **95**, and needs neither hand exclusion: with the object indices
classified correctly, both of those features survive on their own merits. The mechanism for
excluding a feature by name is there, in `[subset] exclude`, and ships empty.

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
with explicit quantifiers, per-feature diagnostic panels, and optional correlation pruning.
Every stage prints a report of what it measured and what it could not.

**Out:** image quality control, CellProfiler, segmentation, per-object pooling, and
row/plate standardization — all of which are tied to a particular plate design and belong
upstream. Also out: UMAP, MMD and Mahalanobis analyses, which consume this tool's output.
[HCS-proc](https://github.com/NIB-SI/HCS-proc) covers all of those.

And out: any claim about experiments other than the one this was built on. The
configuration can describe a different design and the code will run it — that is what
[Describing your experiment](#describing-your-experiment) is about — but running is not
having been validated, and no second dataset has been through it.

## Origin

This is an extraction of `scripts/feature_selection/` from
[NIB-SI/HCS-proc](https://github.com/NIB-SI/HCS-proc), the pipeline published with:

> Tome, M.; Jozef, B.; Mosimann, S. L.; Kosnik, M.; Schirmer, K.; Županič, A.
> *A High-Content Imaging Pipeline to Investigate Subcytotoxic Effects in RTgill-W1 Cells.*
> **Environmental Science & Technology** **2026**, *60* (31), 21402–21416.
> <https://doi.org/10.1021/acs.est.5c18316>

The git history of those scripts is preserved here. HCS-proc remains the citable record of
the published pipeline. This repository takes one stage of it — the concentration-response
feature selection — and makes it installable and configuration-driven, so that the
experimental design is declared in a file instead of edited into five scripts.

Nothing about days, replicates, wells or nine concentrations is hardcoded, and that is
checked rather than asserted: `tests/test_generality.py` runs the whole chain on the same
data re-shaped into a different experiment, from configuration alone. It is a property of
the code. The method itself was developed for one dataset and has been applied to one
dataset — one organism, one assay, one chemistry — and nothing here is evidence about a
second.

## Data

The reference dataset is 121 GB in total and lives on Zenodo:
<https://doi.org/10.5281/zenodo.17951792> (CC-BY-4.0). The single file this tool starts
from is `cell_ID_pooled_median_row_plate_standardization_cid.txt`, 3.9 GB — the untrimmed,
row/plate standardized per-cell table, 503,920 cells x 481 columns. It is deliberately
*untrimmed*: trimming lives in this tool, so a run reproduces that step rather than
inheriting it. One more file from the same record, `all_days_trimmed_features.txt` (1.5 GB),
is the published run's own intermediate and is what the pruning gate checks against.

Nothing that large belongs in git. Put it in `data/`, which is ignored; the committed
fixtures under [`tests/fixtures/`](tests/fixtures/README.md) are small subsets that let the
suite run without it.

## Tests

```bash
pytest
```

Forty-three tests, about twenty seconds, no large data needed. The suite is deliberately
small: every test in it either backs one of the four numbers above, or pins a boundary that
a one-character edit would move, or is the regression for a bug that actually happened.
Tests that re-checked what another test already covered, or that asserted an implementation
detail rather than a promise, were removed rather than carried into the release.

What that leaves, by what it protects:

- **The reproduction gates**, four of them, described above.
- **The trim, against its own oracle.** The vectorised implementation is compared value for
  value with a transcription of the original per-group `np.nanpercentile` loop, on data with
  and without infinities. This is what lets the published population counts be claimed
  without the 3.9 GB file.
- **The boundaries.** The merge cut is exclusive, tested one ULP either side; a constant fit
  that merely *ties* drops its feature; the missing-data rule drops at the threshold and not
  above it. Each of these is one character away from a different published number.
- **The retention rule's asymmetry** — a positive slope on *any* stratum, the constant model
  beaten on *all* — on a table built so that a feature passes on one and fails on the other.
- **The `counts_` trap**, from the failing side: a tidy `^counts_` metadata pattern is shown
  to swallow two features the published run retained.
- **The configuration**, as one sweep over fifty misconfigurations, each asserting the
  message a user would have to act on — and, as its counterweight, that the shipped example
  still validates. Without the second, validation can only ever get stricter until the run
  this package exists to reproduce stops loading.
- **The shipped example itself**, resolved against the real 481-column header.
- **The two bugs that were real**: an unstratified design emptying itself on the way through
  its own output format, and pruning an empty list raising out of the clustering.

### The design comes from the configuration

Every other test is shaped by the RTgill-W1 design, which leaves one property of the code
unchecked — that the experiment comes from the configuration and from nowhere else.
`tests/reshape.py` rebuilds the committed slice of the published table as a different
experiment, keeping the measured values and changing every structural choice at once:

| | published | re-shaped |
|---|---|---|
| stratification | four days | **none** — one unnamed stratum |
| exposure levels | 9 (top one withheld) | **5**, all fitted |
| level labels | `11` control, `10 … 2` | `vehicle` control, `trace, low, mid, high, max` — which sort into a *different* order |
| models | all six | five: 5 points cannot identify BC5's 5 parameters, so the configuration has to say so |
| columns | `Metadata_*`, `Concentration`, `rp_norm_*` | no name in common, features included |
| file | tab-separated, metadata first | comma-separated, columns interleaved, rows shuffled |
| trim | `[p2.5, p97.5]` within day, replicate, well | `[p5, p95]` within batch, sample |

`tests/test_generality.py` runs `check`, `emd`, `fit`, `select`, `prune`, `subset` and
`plot` over it and asserts the answers are the configured ones — that the fitted axis
follows `[design] levels` rather than the labels' sort order, that the outputs take the
configured separator, and that a single unnamed stratum survives every table on disk.

It found two bugs the day it was written, both of which had gone unnoticed by a suite six
times this size. That is what it is for. It is one dataset in two shapes: evidence that the
code reads the design from the configuration, and not evidence about a second experiment.

### The golden regressions

Two of the four gates compare against the published run itself, so they need the large
inputs in `data/`. Each skips when its input is absent, and skips when it is present unless
you opt in, because reading gigabytes takes minutes:

```bash
CDR_FS_GOLDEN=1 pytest tests/test_golden.py tests/test_golden_prune.py -v
```

`test_golden.py` reproduces the published EMD tables from the untrimmed 3.9 GB input: 16,946
treatment distances and 11,292 baseline distances, both population sizes on every row
matching exactly and every distance within 1e-9 relative. The population sizes are the
sharper comparison — they are integers, so no tolerance can hide a difference, and their
matching is what shows the trim step was reproduced value for value.

`test_golden_prune.py` starts instead from the published run's own selected-and-trimmed
subset, `all_days_trimmed_features.txt`, which isolates the pruning stage: its 182 feature
columns are the published selection and its values are already trimmed, so anything that
fails there is pruning and nothing else.

## Licence and citation

MIT — see [LICENSE](LICENSE). Copyright National Institute of Biology.

If you use this tool, please cite both it and the article in which the method was first
published; [CITATION.cff](CITATION.cff) has the machine-readable metadata for both.
