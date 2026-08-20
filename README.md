# cdr_FS: feature selection by concentration/dose–response model fitting

For each feature in a high-dimensional dataset where we expect a concentration/dose response,
we measure how far its distribution moves away from the control as exposure rises; then keep
the features where a concentration–response curve describes that movement better than a flat
line does. The distance is the earth mover's (Wasserstein) distance between populations, and
the curves are six standard concentration–response models ranked by AIC and BIC. The name is
**c**oncentration/**d**ose–**r**esponse **F**eature **S**election.

> The method was developed for, and has been applied to, one dataset: the RTgill-W1 screen of
> [Tome et al. 2026](#origin-and-the-data). The experimental design is read from a
> configuration file rather than written into the code.
> [Reproducing the published run](#reproducing-the-published-run) says which published
> numbers come back.

## The method

Given a per-object table (one row per object, columns split into metadata and measured
features), the selection asks of each feature a single question: **does its distributional
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

**The whole thing rests on one assumption**: that a higher exposure pushes a feature's
distribution further from the control than a lower one does. That is what makes the
distances a series worth fitting a curve to, and what the positive-slope gate tests. A
response that is not monotone in exposure (one that saturates, reverses, or appears only in
a middle band) is not what this looks for.

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
`tests/fixtures/subset.tsv`, 1,272 cells × 30 columns, a real slice of the published
dataset, committed so that there is something to run before you download 121 GB. From the
repository root:

```bash
cdr-fs check -c examples/quickstart.ini --scan
cdr-fs run   -c examples/quickstart.ini
```

About a minute in total, most of it in `fit`. Results land in `results/`.

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

### The two lines to read first

Two lines of `cdr-fs check` carry mistakes nothing downstream can catch. Read them before
running anything, and against your own plate map:

```
exposure axis     10 is the LOWEST exposure, 2 the highest   (11.3683 -> 1000)
[columns]  30 = 10 metadata + 20 feature(s)
```

→ [Quickstart](docs/quickstart.md) says what each one means and what goes wrong.

### What came out

```
30 columns  →  10 metadata + 20 features   (check)
            →  19 reach selection          (fit: one feature never fitted)
            →   9 selected                 (select: 10 no better than a constant)
            →   9 after correlation collapsing
            →   9 after the 30% missing-data filter
```

`results/final_representatives_retained.txt` is the answer, one feature name per line.

**This is a smoke test of the tool, not a reproduction of the method**: it turns trimming off
and fits all nine exposure levels rather than withholding the top one. For real work, copy
[`examples/template.ini`](examples/template.ini) if your dataset is not the published one, and
[`examples/published.ini`](examples/published.ini) if it is.
→ [Quickstart](docs/quickstart.md) annotates the run and the reports it prints.

## The pipeline

One `.ini` file describes a run, and `-c/--config` points at it from anywhere; start from
[`examples/template.ini`](examples/template.ini). Each stage reads the previous one's output
from `[output] dir` under a stable name.

```bash
cdr-fs check        -c config.ini          # validate the configuration, report the schema
cdr-fs check        -c config.ini --scan   # also confirm the design occurs in the data
cdr-fs run          -c config.ini          # the whole chain, in order (see below)
cdr-fs trim         -c config.ini          # write the trimmed table  (optional, see below)
cdr-fs emd          -c config.ini          # distances per feature, stratum and contrast
cdr-fs fit          -c config.ini          # fit the models to each distance series
cdr-fs select       -c config.ini          # apply the retention rule
cdr-fs correlation  -c config.ini          # collapse near-redundant features  (optional)
cdr-fs drop_missing -c config.ini          # drop sparse features, write the final table
cdr-fs plot         -c config.ini          # draw the figures from whatever tables exist
```

| Stage | Writes |
|---|---|
| `emd` | `emd.tsv`: control against each level; `emd_baseline.tsv`: control against control, between replicates |
| `fit` | `fit.tsv` |
| `select` | `selected.txt`; `select_evidence.tsv`: per feature and stratum, which model won, by how much, and what the linear slope was |
| `correlation` | `representatives.txt`; `correlation_clusters.tsv`: every feature with its cluster and its distance to the representative; `correlation_linkage.tsv`: the tree and its leaf order. All three take an `all_` prefix when there was no selection to collapse |
| `drop_missing` | `final_<list>.tsv`: the table restricted to that list; `final_<list>_features.tsv`: how much data each column holds and which rule removed it; `final_<list>_retained.txt` |
| `plot` | `fit_<stratum>_part_<n>.png`; `emd.png` and `emd_baseline.png`; `dendrogram.png`, or `all_dendrogram.png` beside its tree |

### What `run` runs, and what it leaves out

`emd`, `fit`, `select`, `correlation`, `drop_missing`, `plot`, reading the configuration once
so that a bad one fails before the first stage rather than between two of them. Four
departures from "all of them", each stated in the run's own header:

- **`emd`, `fit` and `select` run only when `[select] enabled` is true.** Turned off, no
  concentration–response selection happens and every feature carries forward into the
  filtering stages: the tool as a plain redundancy and sparsity filter. All three go
  together, because `select` needs `fit` and `fit` needs `emd`.
- **`correlation` runs only when `[correlation] enabled` is true.** Turned off, it is
  reported as `off` rather than as a failure, and `drop_missing` falls back to
  `selected.txt`, or to every feature if the selection is off as well.
- **`drop_missing` always runs**, whatever its switch says. The switch decides whether the
  missing-data filter drops anything; the stage writes the final table either way, so a run
  always ends in one.
- **`trim` is never part of a run.** Every stage that needs cell-level data reads
  `[input] table` and applies the configured trim itself, so a trimmed copy is an artefact
  no stage reads. On the reference dataset it would be another 3.9 GB. `cdr-fs trim` writes it
  out for inspection.

`plot` is given only the figures this run's own tables can support, so a report saying
nothing was fitted cannot contain a previous run's distances. With nothing left to draw, `plot`
is dropped from the plan and said to be, and a figure that refuses is reported as `no figures`
rather than failing the run: the run's product is the table, the figures are diagnostics.

### Which file is the answer?

`selected.txt` if you stop after `select`; `representatives.txt` if you collapse correlated
features; `final_<list>_retained.txt` if you also apply the missing-data filter. Each is one
feature name per line, and `<list>` is whichever list narrowed the features last.

A correlation of every feature is not a correlation of the selected ones, so with the selection
off those outputs take an `all_` prefix. The list becomes `all_representatives.txt`, and a run
ends in `final_all_representatives_retained.txt`. Two chains can therefore share one
`[output] dir` without either overwriting the other's answer:

| | `[select]` on | `[select]` off |
|---|---|---|
| `[correlation]` on | `representatives.txt` → `final_representatives*` | `all_representatives.txt` → `final_all_representatives*` |
| `[correlation]` off | `selected.txt` → `final_selected*` | (none) → `final_all*` |

## Reproducing the published run

[`examples/published.ini`](examples/published.ini) is the configuration for the dataset the
method was published on. Four of its published outputs come back:

| Stage | Published output | Result |
|---|---|---|
| `emd` | the two EMD tables, 16,946 treatment and 11,292 baseline distances | reproduced; both population sizes exact on every row, distances within 8.5e-13 relative |
| `select` | the two retained feature lists, **182** across all days and **374** on D5 alone | both reproduced as identical sets |
| `correlation` | the all-days list after correlation collapsing, **99** features | reproduced, and its composition matches the published categorization |
| `drop_missing` | the final retained list, **95** features | reproduced: 94 `rp_norm_*` plus `counts_RelateLysoCell` |

The 182 and the 99 belong to the published metadata split, which carried CellProfiler's
object-index columns as features. With those declared metadata, as `examples/published.ini`
does, the same rules give 463 features, then 175, then 97, then the same **95**.
[Reproducing the published run](docs/reproducing.md) sets the two routes out side by side,
says how each check is arranged, and how to run them.

## Documentation

| Page | What is on it |
|---|---|
| [Quickstart](docs/quickstart.md) | the first run annotated: what `check` tells you, what the reports mean, where their numbers come from, and three outputs that look wrong and are not |
| [Configuration](docs/configuration.md) | every key, what it defaults to, and the two patterns worth getting right |
| [Describing your experiment](docs/experiment-design.md) | which keys state your design, what to write when a piece of it is missing, and what the format cannot express |
| [Method notes](docs/method-notes.md) | trimming, pooling replicates, correlation collapsing, the missing-data filter, and the figures |
| [Reproducing the published run](docs/reproducing.md) | what each of the four gates checks, and what the numbers are and are not |
| [Troubleshooting](docs/troubleshooting.md) | what a refusal means, what each summary status means, and where to look when a run retains nothing |

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
path (`git log --oneline --follow -- legacy/plots_emd_model_drc.py`).

**The data is that work's data**, published alongside the article on Zenodo:
<https://doi.org/10.5281/zenodo.17951792> (CC-BY-4.0), 121 GB in total. The single file this
tool starts from is `cell_ID_pooled_median_row_plate_standardization_cid.txt` (3.9 GB), the
untrimmed, row/plate standardized per-cell table, 503,920 cells × 481 columns. It is
deliberately *untrimmed*: trimming lives in this tool, so a run reproduces that step rather
than inheriting it. Put it in `data/`, which is gitignored, and see
[Reproducing the published run](#reproducing-the-published-run).

The quickstart needs none of that. It runs on `tests/fixtures/subset.tsv`, a committed
1,272-cell slice of the same table. See [`tests/fixtures/`](tests/fixtures/README.md).

## Licence and citation

MIT (see [LICENSE](LICENSE)). Copyright National Institute of Biology.

If you use this tool, please cite both it and the article in which the method was first
published; [CITATION.cff](CITATION.cff) has the machine-readable metadata for both.
