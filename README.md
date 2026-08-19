# cdr_FS — feature selection by concentration/dose–response model fitting

For each morphological feature in a high-content screen, measure how far its distribution
moves away from the control as exposure rises; then keep the features where a
concentration–response curve describes that movement better than a flat line does. The
distance is the earth mover's (Wasserstein) distance between cell populations, and the
curves are six standard concentration–response models ranked by AIC and BIC. The name is
**c**oncentration/**d**ose–**r**esponse **F**eature **S**election.

> The method was developed for, and has been applied to, one dataset — the RTgill-W1 screen
> of [Tome et al. 2026](#origin). The experimental design is read from a configuration file
> rather than written into the code, and that is tested; it is a statement about the code,
> not evidence that the method carries to another organism, assay or chemistry.
> [Reproducing the published run](#reproducing-the-published-run) says which published
> numbers come back.

## The method

Given a per-object table — one row per cell, columns split into metadata and measured
features — the selection asks of each feature a single question: **does its distributional
distance from the control grow with exposure in a way a concentration–response curve
describes better than a flat line?**

1. For every feature and every exposure level, compute the earth mover's (Wasserstein)
   distance between the control population and the exposed population. The distance between
   control populations of different biological replicates gives the reproducibility floor to
   read those numbers against.
2. Fit six models to the resulting distance-versus-exposure series: Brain–Cousens hormesis
   (BC4, BC5), four-parameter log-logistic (LL4), four-parameter Weibull (WB1.4), linear
   (Lin) and constant (Con).
3. Retain a feature when its linear slope is positive **and** the constant model is not the
   best fit by AIC + BIC. A feature whose distance from the control is flat, or shrinks,
   carries no concentration-dependent signal.
4. Optionally collapse the survivors by correlation, keeping one representative per cluster
   of near-redundant features, and drop any left too sparse to use.

Using distributions rather than per-well means is the point: a subpopulation can shift while
the population mean stays put, and a distance between distributions registers that where a
difference of means does not.

**The whole thing rests on one assumption** — that a higher exposure pushes a feature's
distribution further from the control than a lower one does. That is what makes the
distances a series worth fitting a curve to, and what the positive-slope gate tests. A
response that is not monotone in exposure — one that saturates, reverses, or appears only in
a middle band — is not what this looks for.

**And the distance is unsigned.** A feature whose values *fall* with exposure is retained
exactly like one whose values rise: "positive slope" describes the distance from the control
growing, not the measurement growing. Which way a retained feature moved is a question for
the data, and this tool does not answer it.

## Installation

```bash
pip install -e .
```

Python 3.10 or newer. Requires numpy, pandas, scipy and matplotlib; `pip install -e ".[dev]"`
adds pytest.

## Quickstart

The repository ships a configuration you can run without editing anything. It points at
`tests/fixtures/subset.tsv` — 1,272 cells × 30 columns, a real slice of the published
dataset, committed so that there is something to run before you download 121 GB. From the
repository root:

```bash
cdr-fs check  -c examples/quickstart.ini --scan
cdr-fs emd    -c examples/quickstart.ini
cdr-fs fit    -c examples/quickstart.ini
cdr-fs select -c examples/quickstart.ini
cdr-fs prune  -c examples/quickstart.ini
cdr-fs subset -c examples/quickstart.ini
cdr-fs plot   -c examples/quickstart.ini --only fits --features results/pruned.txt
```

About fifty seconds in total, most of it in `fit`. Results land in `results/`.

**Start with `check`.** It reads the configuration and the table's header, and prints how the
columns resolved. This is the line to look at:

```
[columns]  30 = 10 metadata + 20 feature(s)
  metadata          Concentration, counts_Cells, counts_Cytoplasm, ...
  features          rp_* 18, counts_* 2
    counts_*        counts_RelateLysoCell, counts_RelateMitoCell
```

If that feature count is not what you expect, stop and fix `[schema] metadata_patterns`
before running anything else — see [Traps](#traps). `--scan` additionally reads the exposure
and stratum columns to confirm the levels you declared actually occur in the data.

**Then the chain.** Each stage prints what it measured and what it could not:

```
$ cdr-fs emd -c examples/quickstart.ini
reading tests/fixtures/subset.tsv ...
  1,272 row(s), 20 feature(s)
714 distance(s) over 36 comparison(s) x 20 feature(s)
  6 (feature, comparison) cell(s) skipped - a population held no values: rp_norm_Mean_PunctaLyso_...
wrote results/emd.tsv  (60.72 KiB)
458 distance(s) over 24 comparison(s) x 20 feature(s)
wrote results/emd_baseline.tsv  (45.50 KiB)

$ cdr-fs fit -c examples/quickstart.ini
fitted 76 series of 9 point(s) over 20 feature(s) x 4 stratum/strata
  4 series not fitted for want of a complete exposure series, across 1 feature(s): rp_norm_Mean_...
  20 fit(s) did not converge: BC4 3, BC5 3, LL4 6, WB1.4 8
wrote results/fit.tsv  (81.21 KiB)

$ cdr-fs select -c examples/quickstart.ini
retained 9 of 19 feature(s)
  rule: positive slope on any stratum/strata, constant model beaten on all, ranked by aic_plus_bic
  strata: D9, D1, D5, D7
  rejected: 0 for slope, 10 for being no better than constant
wrote results/select_evidence.tsv  (9.20 KiB)
wrote results/selected.txt  (9 feature(s))
```

Those counts are worth checking against the design rather than taking on trust. 36
comparisons is 4 strata × 9 exposure levels; 36 × 20 features, less 6 cells where a
population was empty, is 714 distances. The baseline's 24 is 4 strata × the 6 unordered pairs
of 4 replicates. 76 series is 20 × 4 less the 4 that were incomplete, and `select` says "of
19" rather than "of 20" because one feature was never fitted on any stratum — `fit` names it
in the line above.

**What came out.** The funnel for this run:

```
30 columns  →  10 metadata + 20 features   (check)
            →  19 reach selection          (fit: one feature never fitted)
            →   9 selected                 (select: 10 no better than a constant)
            →   9 after correlation pruning
            →   9 after the 30% missing-data filter
```

`results/subset_pruned_retained.txt` is the answer — one feature name per line. When a run
retains nothing, `select_evidence.tsv` is where to look first: it carries the linear slope
and the winning model for every feature and stratum.

Three things about this run that are normal and look alarming:

- **`prune` removed nothing.** With 20 features and a `|r| ≥ 0.9` threshold, no pair is
  redundant. On the published 471 it collapses 175 features to 97.
- **`fit` prints a scipy `OptimizeWarning` to stderr.** Fits that genuinely fail are counted
  on stdout — "20 fit(s) did not converge" — and are expected; some shapes do not fit some
  series.
- **`plot` without `--features` draws every feature in the distance table.** Here that is 15
  figures and 8 MB; on a full run it is hundreds of pages.

This configuration is a smoke test of the tool, not a reproduction of the method: it turns
trimming off, because the fixture holds one to eight cells per well and a within-well
percentile on two values discards both, and it fits all nine exposure levels rather than
withholding the top one. [`examples/published.ini`](examples/published.ini) is the method as
it was published, and it is the file to copy for real work.


## The pipeline

One `.ini` file describes a run, and `-c/--config` points at it from anywhere. Each stage
reads the previous one's output from `[output] dir` under a stable name.

```bash
cdr-fs check  -c config.ini          # validate the configuration, report the schema
cdr-fs check  -c config.ini --scan   # also confirm the design occurs in the data
cdr-fs trim   -c config.ini          # write the trimmed table  (optional — see below)
cdr-fs emd    -c config.ini          # distances per feature, stratum and contrast
cdr-fs fit    -c config.ini          # fit the models to each distance series
cdr-fs select -c config.ini          # apply the retention rule
cdr-fs prune  -c config.ini          # collapse near-redundant features   (optional)
cdr-fs subset -c config.ini          # write the table restricted to what survived
cdr-fs plot   -c config.ini          # draw the figures from whatever tables exist
```

| Stage | Writes |
|---|---|
| `emd` | `emd.tsv` — control against each level; `emd_baseline.tsv` — control against control, between replicates |
| `fit` | `fit.tsv` |
| `select` | `selected.txt`; `select_evidence.tsv` — per feature and stratum, which model won, by how much, and what the linear slope was |
| `prune` | `pruned.txt`; `prune_clusters.tsv` — every feature with its cluster and its distance to the representative; `prune_linkage.tsv` — the tree and its leaf order |
| `subset` | `subset_<list>.tsv` — the table restricted to that list; `subset_<list>_features.tsv` — how much data each column holds and which rule removed it; `subset_<list>_retained.txt` |
| `plot` | `fit_<stratum>_part_<n>.png`; `emd.png` and `emd_baseline.png`; `dendrogram.png` |

**`trim` is optional.** Every stage that needs cell-level data reads `[input] table` and
applies the configured trim itself, so there is no need to materialise a trimmed copy — on
the reference dataset that copy would be another 3.9 GB. `cdr-fs trim` exists to write it out
for inspection.

**Which file is the answer?** `selected.txt` if you stop after `select`; `pruned.txt` if you
prune; `subset_<list>_retained.txt` if you also apply the missing-data filter. Each is one
feature name per line. When a run retains nothing, `select_evidence.tsv` is where to look
first — it carries the slope and the winning model for every feature and stratum.

A stage refuses rather than write an artefact that would read as a result: `fit` will not
write a table when no series was complete, and `prune` and `subset` will not run on an empty
feature list. Exit codes are 0 for success, 2 for a configuration error, 3 for nothing to do.

One module per stage, named for it: `emd.py`, `fit.py`, `select.py`, `prune.py`, `subset.py`,
`trim.py`, `plots.py`, with `config.py` and `schema.py` underneath.

## Traps

Six things that are quiet when they go wrong. Each is explained fully on the page it links to.

**A `^counts_` metadata pattern silently discards features.** `[schema] metadata_patterns`
names the *metadata*, and everything unmatched is a feature — so an over-broad pattern
removes real measurements. In the reference dataset five columns begin with `counts_` and
only three are metadata:

```
counts_Cells, counts_Cytoplasm, counts_FilteredNuclei     segmentation QC  -> metadata
counts_RelateLysoCell, counts_RelateMitoCell              organelle counts -> FEATURES
```

`counts_RelateLysoCell` survives the entire selection into the published retained list. Run
`cdr-fs check` and read the `[columns]` line before trusting a pattern.
→ [Configuration](docs/configuration.md#one-warning-worth-reading-before-you-write-metadata_patterns)

**The mistake runs the other way too.** A column can be named like a measurement and be a
label: CellProfiler's `Number_Object_Number` is an object index, and it responds to exposure
because objects-per-image changes with dose. Only a pattern you write keeps it out.
→ [Configuration](docs/configuration.md)

**`[design] levels` must run low to high, and the labels do not tell the tool which end is
which.** A reversed list runs to completion, produces an identical distance table, negates
every linear slope, and retains nothing. Supplying `[design] dose` catches the case where the
two lists disagree — a dose vector that falls while the levels are declared to rise is
refused — but it cannot catch both being wrong the same way, because only you know which
label was your top dose. So the real check is one line of `cdr-fs check`:

```
exposure axis     10 is the LOWEST exposure, 2 the highest   (11.3683 -> 1000)
```

Read it against your own plate map before you trust a run.
→ [Experiment design](docs/experiment-design.md)

**Read the stage tables with `cdr_fs.schema.read_stage_table`, not a bare
`pandas.read_csv`.** An experiment with no `[schema] group_by` has one stratum whose label is
the empty string; pandas reads an empty field as `NaN`, and `groupby` then drops it — taking
the whole table with it. → [Testing](docs/testing.md)

**`[select] strata` narrows the whole run, not only the selection.** `emd` and `fit` read it
too, so setting it to one stratum leaves a distance table holding only that stratum, with
nothing in the file to say it is partial. → [Experiment design](docs/experiment-design.md)

**`[prune] aggregate_by` must name a unit that varies along the exposure axis.** Correlations
are computed between unit medians, so if each unit spans the whole dilution series then every
unit looks alike and pruning collapses nothing. Check the `N unit(s) of …` line in the
report. → [Method notes](docs/method-notes.md)

## Reproducing the published run

[`examples/published.ini`](examples/published.ini) is the configuration for the dataset the
method was published on. Four published outputs are checked against, each stage given the
published input to the stage before it rather than this tool's own output, so that a
difference is attributable to the stage under test:

| Stage | Published output | Result |
|---|---|---|
| `emd` | the two EMD tables — 16,946 treatment and 11,292 baseline distances | reproduced; both population sizes on every row exact, distances within 8.5e-13 relative |
| `select` | the two retained feature lists — **182** across all days, **374** on D5 alone | both reproduced as identical sets |
| `prune` | the all-days list after pruning — **99** features | reproduced, and its composition matches the published categorization across all 4 organelles × 7 measurement families |
| `subset` | the final retained list — **95** features | reproduced: 94 `rp_norm_*` plus `counts_RelateLysoCell` |

The `select` gate runs off the published fit table and needs only a 1 MB committed fixture,
so it runs in CI; the others need the large inputs and are opt-in.
→ [Reproducing the published run](docs/reproducing.md)

## Data

The reference dataset is 121 GB in total and lives on Zenodo:
<https://doi.org/10.5281/zenodo.17951792> (CC-BY-4.0). The single file this tool starts from
is `cell_ID_pooled_median_row_plate_standardization_cid.txt`, 3.9 GB — the untrimmed,
row/plate standardized per-cell table, 503,920 cells × 481 columns. It is deliberately
*untrimmed*: trimming lives in this tool, so a run reproduces that step rather than
inheriting it.

Nothing that large belongs in git. Put it in `data/`, which is ignored; the committed
fixtures under [`tests/fixtures/`](tests/fixtures/README.md) are small subsets that let the
suite — and the quickstart above — run without it.

## Tests

```bash
pytest
```

Forty-three tests, about twenty seconds, no large data needed. Two of the four reproduction
gates need the Zenodo files and are opt-in via `CDR_FS_GOLDEN=1`.
→ [Testing](docs/testing.md)

## Documentation

| Page | What is on it |
|---|---|
| [Configuration](docs/configuration.md) | every key, what it defaults to, and the two patterns worth getting right |
| [Experiment design](docs/experiment-design.md) | which keys state your design, what to write when a piece of it is missing, and what the format cannot express |
| [Method notes](docs/method-notes.md) | pooling replicates, correlation pruning, the missing-data filter, and the figures |
| [Reproducing the published run](docs/reproducing.md) | what each of the four gates checks, and what the numbers are and are not |
| [Testing](docs/testing.md) | what the suite covers and why it is small |

## Scope

**In:** the schema declaration, optional trimming, contrast definition, the EMD computation,
multi-model fitting with information-criterion ranking, the retention rule with explicit
quantifiers, per-feature diagnostic panels, and optional correlation pruning. Every stage
prints a report of what it measured and what it could not.

**Out:** image quality control, CellProfiler, segmentation, per-object pooling, and row/plate
standardization — all tied to a particular plate design and belonging upstream. Also out:
UMAP, MMD and Mahalanobis analyses, which consume this tool's output.
[HCS-proc](https://github.com/NIB-SI/HCS-proc) covers all of those.

## Origin

This is an extraction of `scripts/feature_selection/` from
[NIB-SI/HCS-proc](https://github.com/NIB-SI/HCS-proc), the pipeline published with:

> Tome, M.; Jozef, B.; Mosimann, S. L.; Kosnik, M.; Schirmer, K.; Županič, A.
> *A High-Content Imaging Pipeline to Investigate Subcytotoxic Effects in RTgill-W1 Cells.*
> **Environmental Science & Technology** **2026**, *60* (31), 21402–21416.
> <https://doi.org/10.1021/acs.est.5c18316>

HCS-proc remains the citable record of the published pipeline. This repository takes one
stage of it and makes it installable and configuration-driven, so that the experimental
design is declared in a file instead of edited into five scripts. The original author's
history is preserved here; because git does not score the extraction as a rename, reach it
through the old path — `git log --oneline -- legacy/plots_emd_model_drc.py`.

## Licence and citation

MIT — see [LICENSE](LICENSE). Copyright National Institute of Biology.

If you use this tool, please cite both it and the article in which the method was first
published; [CITATION.cff](CITATION.cff) has the machine-readable metadata for both.
