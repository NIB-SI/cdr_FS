# cdr_FS — feature selection by concentration/dose–response model fitting

For each morphological feature in a high-content screen, measure how far its distribution
moves away from the control as exposure rises; then keep the features where a
concentration–response curve describes that movement better than a flat line does. The
distance is the earth mover's (Wasserstein) distance between cell populations, and the
curves are six standard concentration–response models ranked by AIC and BIC. The name is
**c**oncentration/**d**ose–**r**esponse **F**eature **S**election.

> The method was developed for, and has been applied to, one dataset — the RTgill-W1 screen
> of [Tome et al. 2026](#origin-and-the-data). The experimental design is read from a configuration file
> rather than written into the code — a statement about the code, not evidence that the
> method carries to another organism, assay or chemistry.
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

### The limits of the question

**The whole thing rests on one assumption** — that a higher exposure pushes a feature's
distribution further from the control than a lower one does. That is what makes the
distances a series worth fitting a curve to, and what the positive-slope gate tests. A
response that is not monotone in exposure — one that saturates, reverses, or appears only in
a middle band — is not what this looks for.

**And the distance is unsigned.** A feature whose values *fall* with exposure is retained
exactly like one whose values rise: "positive slope" describes the distance from the control
growing, not the measurement growing. Which way a retained feature moved is a question for
the data, and this tool does not answer it.

### What it does not do

It starts from a per-object table and stops at a list of feature names. Image quality control,
CellProfiler, segmentation, per-object pooling and row/plate standardization sit upstream, all
tied to a particular plate design; the UMAP, MMD and Mahalanobis analyses that consume the
selected features sit downstream. [HCS-proc](https://github.com/NIB-SI/HCS-proc) covers those.

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
cdr-fs check -c examples/quickstart.ini --scan
cdr-fs run   -c examples/quickstart.ini
```

About fifty seconds in total, most of it in `fit`. Results land in `results/`.

`run` is the chain in one command, and every stage is also a command of its own. These six
are the same run:

```bash
cdr-fs emd          -c examples/quickstart.ini
cdr-fs fit          -c examples/quickstart.ini
cdr-fs select       -c examples/quickstart.ini
cdr-fs correlation  -c examples/quickstart.ini
cdr-fs drop_missing -c examples/quickstart.ini
cdr-fs plot         -c examples/quickstart.ini --features results/final_representatives_retained.txt
```

### Start with `check`

It reads the configuration and the table's header, and prints how the columns resolved. This
is the line to look at:

```
[columns]  30 = 10 metadata + 20 feature(s)
  metadata          Concentration, counts_Cells, counts_Cytoplasm, ...
  features          rp_* 18, counts_* 2
    counts_*        counts_RelateLysoCell, counts_RelateMitoCell
```

If that feature count is not what you expect, stop and fix `[schema] metadata_patterns` before
running anything else. The patterns name the *metadata*, and everything they fail to match is a
feature — so a pattern that matches too much quietly removes real measurements, and nothing
later says so. Five of the reference dataset's columns begin with `counts_` and only three are
metadata: a tidy `^counts_` discards `counts_RelateLysoCell`, a feature the published run kept.
→ [Configuration](docs/configuration.md#one-warning-worth-reading-before-you-write-metadata_patterns)

`check` prints a second line worth reading, for the same reason — nothing else can catch this
one either:

```
exposure axis     10 is the LOWEST exposure, 2 the highest   (11.3683 -> 1000)
```

`[design] levels` is the response axis and it runs low to high, but the labels themselves tell
the tool nothing: the reference dataset's own labels count *down* as exposure climbs. A list
written the wrong way round runs to completion, produces an identical distance table, negates
every slope and retains nothing. Read that line against your own plate map.
→ [Describing your experiment](docs/experiment-design.md)

`--scan` additionally reads the exposure and stratum columns to confirm the levels you
declared actually occur in the data.

### Then the chain

`run` executes the stages in order and prints each one's own report under a rule. Each stage
says what it measured and what it could not:

```
$ cdr-fs emd -c examples/quickstart.ini
reading tests/fixtures/subset.tsv ...
  1,272 row(s), 20 feature(s)
714 distance(s) over 36 comparison(s) x 20 feature(s)
  6 (feature, comparison) cell(s) skipped - a population held no values: rp_norm_Mean_PunctaLyso_...
wrote results/emd.tsv  (60.72 KiB)
458 distance(s) over 24 comparison(s) x 20 feature(s)
  22 (feature, comparison) cell(s) skipped - a population held no values: rp_norm_Mean_PunctaLyso_...
wrote results/emd_baseline.tsv  (45.50 KiB)

$ cdr-fs fit -c examples/quickstart.ini
fitted 76 series of 9 point(s) over 20 feature(s) x 4 stratum/strata
  4 series not fitted for want of a complete exposure series, across 1 feature(s): rp_norm_Mean_...
  1 feature(s) were never fitted on any stratum and so cannot be selected: rp_norm_Mean_PunctaLyso_...
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

Every count in those reports is a product of your design, so check them against it rather
than taking them on trust: 36 comparisons is 4 strata × 9 exposure levels, and the baseline's
24 is 4 strata × the 6 unordered pairs of 4 replicates. Where a number is short, the line
above says why — `select` says "of 19" rather than "of 20" because `fit` could not fit one
feature on any stratum.

A run closes on a ledger of what ran and the file that answers the question:

```
-- summary -------------------------------------------------------------------
  emd               ok
  fit               ok
  select            ok
  correlation       ok
  drop_missing      ok
  plot              ok
  trim              not run
retained 9 feature(s) -> results/final_representatives_retained.txt
```

### What came out

The funnel for this run:

```
30 columns  →  10 metadata + 20 features   (check)
            →  19 reach selection          (fit: one feature never fitted)
            →   9 selected                 (select: 10 no better than a constant)
            →   9 after correlation collapsing
            →   9 after the 30% missing-data filter
```

`results/final_representatives_retained.txt` is the answer — one feature name per line. When a run
retains nothing, `select_evidence.tsv` is where to look first: it carries the linear slope
and the winning model for every feature and stratum.

### Three things about this run that are normal and look alarming

- **`correlation` removed nothing.** With 20 features and a `|r| ≥ 0.9` threshold, no pair is
  redundant. On the reference dataset it collapses 175 selected features to 97.
- **`fit` prints a scipy `OptimizeWarning` to stderr.** Fits that genuinely fail are counted
  on stdout — "20 fit(s) did not converge" — and are expected; some shapes do not fit some
  series.
- **`plot` on its own draws every feature in the distance table.** Here that is 15 figures
  and 8 MB; on a full run it is hundreds of pages. Inside `run` it is given the retained
  list, which is 7 figures here.

**This configuration is a smoke test of the tool, not a reproduction of the method.** It
turns trimming off, because the fixture holds one to eight cells per well and a within-well
percentile on two values discards both, and it fits all nine exposure levels rather than
withholding the top one. The method as it was published is
[`examples/published.ini`](examples/published.ini). For a dataset that is not that one, copy
[`examples/template.ini`](examples/template.ini) instead: the same ten sections, every switch
written out beside its default, and a comment on each saying what it means for the
experiment.

## The pipeline

One `.ini` file describes a run — start from [`examples/template.ini`](examples/template.ini)
— and `-c/--config` points at it from anywhere. Each stage reads the previous one's output
from `[output] dir` under a stable name.

```bash
cdr-fs check        -c config.ini          # validate the configuration, report the schema
cdr-fs check        -c config.ini --scan   # also confirm the design occurs in the data
cdr-fs run          -c config.ini          # the six stages below, in order
cdr-fs trim         -c config.ini          # write the trimmed table  (optional — see below)
cdr-fs emd          -c config.ini          # distances per feature, stratum and contrast
cdr-fs fit          -c config.ini          # fit the models to each distance series
cdr-fs select       -c config.ini          # apply the retention rule
cdr-fs correlation  -c config.ini          # collapse near-redundant features   (optional)
cdr-fs drop_missing -c config.ini          # drop sparse features, write the final table
cdr-fs plot         -c config.ini          # draw the figures from whatever tables exist
```

| Stage | Writes |
|---|---|
| `emd` | `emd.tsv` — control against each level; `emd_baseline.tsv` — control against control, between replicates |
| `fit` | `fit.tsv` |
| `select` | `selected.txt`; `select_evidence.tsv` — per feature and stratum, which model won, by how much, and what the linear slope was |
| `correlation` | `representatives.txt`; `correlation_clusters.tsv` — every feature with its cluster and its distance to the representative; `correlation_linkage.tsv` — the tree and its leaf order |
| `drop_missing` | `final_<list>.tsv` — the table restricted to that list; `final_<list>_features.tsv` — how much data each column holds and which rule removed it; `final_<list>_retained.txt` |
| `plot` | `fit_<stratum>_part_<n>.png`; `emd.png` and `emd_baseline.png`; `dendrogram.png` |

### What `run` runs, and what it leaves out

`emd`, `fit`, `select`, `correlation`, `drop_missing`, `plot`, reading the configuration once
so that a bad one fails before the first stage rather than between two of them. Three
departures from "all of them", each stated in the run's own header:

- **`correlation` runs only when `[correlation] enabled` is true.** Turned off it is reported
  as skipped, and `drop_missing` applies `selected.txt` instead of `representatives.txt`.
- **`drop_missing` always runs**, whatever its switch says. The switch decides whether the
  missing-data filter drops anything; the stage writes the final table either way, so a run
  always ends in one.
- **`trim` is never part of a run**, for the reason below.

### `trim` is optional

Every stage that needs cell-level data reads `[input] table` and applies the configured trim
itself, so there is no need to materialise a trimmed copy — on the reference dataset that
copy would be another 3.9 GB. `cdr-fs trim` exists to write it out for inspection, which is
why `run` does not call it.

### Which file is the answer?

`selected.txt` if you stop after `select`; `representatives.txt` if you collapse correlated
features; `final_<list>_retained.txt` if you also apply the missing-data filter. Each is one
feature name per line.

### When a stage refuses

A stage refuses rather than write an artefact that would read as a result: `fit` will not
write a table when no series was complete, and `correlation` and `drop_missing` will not run
on an empty feature list. Exit codes are 0 for success, 2 for a configuration error, 3 for
nothing to do. Inside `run` a refusal stops the chain: the stage's message is printed, the
summary marks the rest `not reached`, and the run exits 3. A stage the configuration turned
off is not a refusal — it is a declared outcome, and the run still exits 0.

There is one module per stage, named for it — `emd.py`, `fit.py`, `select.py`,
`correlation.py`, `drop_missing.py`, `trim.py`, `plots.py` — with `config.py` and `schema.py`
underneath.

## Reproducing the published run

[`examples/published.ini`](examples/published.ini) is the configuration for the dataset the
method was published on. Four published outputs are checked against, each stage given the
published input to the stage before it rather than this tool's own output, so that a
difference is attributable to the stage under test:

| Stage | Published output | Result |
|---|---|---|
| `emd` | the two EMD tables — 16,946 treatment and 11,292 baseline distances | reproduced; both population sizes on every row exact, distances within 8.5e-13 relative |
| `select` | the two retained feature lists — **182** across all days, **374** on D5 alone | both reproduced as identical sets |
| `correlation` | the all-days list after the correlation stage — **99** features | reproduced, and its composition matches the published categorization across all 4 organelles × 7 measurement families |
| `drop_missing` | the final retained list — **95** features | reproduced: 94 `rp_norm_*` plus `counts_RelateLysoCell` |

The **182** and the **99** belong to the published metadata split, which carried CellProfiler's
object-index columns as features; with those declared metadata, as `examples/published.ini`
does, the same rules give 463 features → 175 → 97 → **95**. Both routes end at the same final
list. [Reproducing the published run](docs/reproducing.md) sets the two out side by side.

The `select` gate runs off the published fit table and needs only a 1 MB committed fixture,
so it needs no download; the others read the Zenodo files and are opt-in.
→ [Reproducing the published run](docs/reproducing.md)

### Checking it yourself

```bash
pytest                                    # the checks that need no data, about two seconds
CDR_FS_GOLDEN=1 pytest                    # and the ones that read the Zenodo files
```

The suite is deliberately small: it exists to hold these four numbers, and little else. The
gates that compare against the published run itself need the dataset in `data/`, so they skip
unless you opt in.

## Documentation

| Page | What is on it |
|---|---|
| [Configuration](docs/configuration.md) | every key, what it defaults to, and the two patterns worth getting right |
| [Describing your experiment](docs/experiment-design.md) | which keys state your design, what to write when a piece of it is missing, and what the format cannot express |
| [Method notes](docs/method-notes.md) | pooling replicates, correlation collapsing, the missing-data filter, and the figures |
| [Reproducing the published run](docs/reproducing.md) | what each of the four gates checks, and what the numbers are and are not |

## Origin, and the data

This is an extraction of `scripts/feature_selection/` from
[NIB-SI/HCS-proc](https://github.com/NIB-SI/HCS-proc), the pipeline published with:

> Tome, M.; Jozef, B.; Mosimann, S. L.; Kosnik, M.; Schirmer, K.; Županič, A.
> *A High-Content Imaging Pipeline to Investigate Subcytotoxic Effects in RTgill-W1 Cells.*
> **Environmental Science & Technology** **2026**, *60* (31), 21402–21416.
> <https://doi.org/10.1021/acs.est.5c18316>

HCS-proc remains the citable record of the published pipeline. This repository takes one stage
of it and makes it installable and configuration-driven, so that the experimental design is
declared in a file instead of edited into five scripts. The original author's history is
preserved here; because git does not score the extraction as a rename, reach it through the old
path — `git log --oneline --follow -- legacy/plots_emd_model_drc.py`.

**The data is that work's data**, published alongside the article on Zenodo:
<https://doi.org/10.5281/zenodo.17951792> (CC-BY-4.0), 121 GB in total. The single file this
tool starts from is `cell_ID_pooled_median_row_plate_standardization_cid.txt`, 3.9 GB — the
untrimmed, row/plate standardized per-cell table, 503,920 cells × 481 columns. It is
deliberately *untrimmed*: trimming lives in this tool, so a run reproduces that step rather
than inheriting it. Put it in `data/`, which is gitignored, and see
[Reproducing the published run](#reproducing-the-published-run).

The quickstart needs none of that. It runs on `tests/fixtures/subset.tsv`, a committed
1,272-cell slice of the same table — see [`tests/fixtures/`](tests/fixtures/README.md).

## Licence and citation

MIT — see [LICENSE](LICENSE). Copyright National Institute of Biology.

If you use this tool, please cite both it and the article in which the method was first
published; [CITATION.cff](CITATION.cff) has the machine-readable metadata for both.
