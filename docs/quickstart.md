# Quickstart

An annotated first run, to read once before you point the tool at your own data: what
`cdr-fs check` tells you before you start, what the stage reports mean, and where their numbers
come from. The [README](../README.md) carries the commands and the shape of the answer.

Nothing here needs a download. [`examples/quickstart.ini`](../examples/quickstart.ini) runs as
it stands, against `tests/fixtures/subset.tsv`, a committed 1,272-cell slice of the published
dataset. From the repository root, about a minute in total and most of it in `fit`:

```bash
cdr-fs check -c examples/quickstart.ini --scan
cdr-fs run   -c examples/quickstart.ini
```

## Start with `check`

It reads the configuration and the table's header, and prints how the columns resolved. This is
the line to look at:

```
[columns]  30 = 10 metadata + 20 feature(s)
  metadata          Concentration, counts_Cells, counts_Cytoplasm, ...
  features          rp_* 18, counts_* 2
    counts_*        counts_RelateLysoCell, counts_RelateMitoCell
```

If that feature count is not the one you expect, stop and fix `[schema] metadata_patterns`
before running anything else. The patterns name the *metadata*, and everything they fail to
match is a feature, so a pattern that matches too much quietly removes real measurements and
nothing later says so. Both directions of that mistake are set out under
[Configuration](configuration.md#one-warning-worth-reading-before-you-write-metadata_patterns).

`check` prints a second line worth reading, for the same reason. Nothing else can catch this
one either:

```
exposure axis     10 is the LOWEST exposure, 2 the highest   (11.3683 -> 1000)
```

`[design] levels` is the response axis and runs low to high, but the labels tell the tool
nothing; this dataset's own labels count *down* as exposure climbs. Read the line against your
plate map: a list written the wrong way round runs to completion and comes out with a short,
plausible answer. See [Describing your experiment](experiment-design.md#when-a-piece-of-the-design-is-missing).

`--scan` additionally reads the exposure and stratum columns, to confirm that the levels you
declared occur in the data.

## Then the chain

`run` executes the stages in order, printing each one's own report under a rule; every stage is
also a command of its own. Each says what it measured and what it could not:

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

Every count is a product of your design, so check them against it rather than taking them on
trust. The 36 comparisons are 4 strata x 9 exposure levels, and the 714 distances are 36 x 20
features less the 6 cells where a population held no values; the baseline's 24 comparisons are
4 strata x the 6 unordered pairs of 4 replicates. Where a number is short the line above says
why: `select` reports "of 19" rather than "of 20" because `fit` could not fit one feature on any
stratum.

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
retained 9 feature(s)
  feature list      results/final_representatives_retained.txt
  final table       results/final_representatives.tsv
```

## Three things about this run that are normal and look alarming

- **`correlation` removed nothing.** Nine features reach it, and at `|r| >= 0.9` no pair among
  them is redundant. On the reference dataset the stage collapses 175 selected features to 97.
- **`fit` prints a scipy `OptimizeWarning` to stderr.** Fits that genuinely fail are counted on
  stdout, as "20 fit(s) did not converge", and are expected: some shapes do not fit some series.
- **`plot` on its own draws every feature in the distance table.** Here that is 15 figures and
  8 MB, on a full run hundreds of pages. Inside `run` it is given the retained list: 7 figures.

## What this run is not

**This configuration is a smoke test of the tool, not a reproduction of the method.** It turns
trimming off, because the fixture holds one to eight cells per well and a within-well percentile
on two values discards both, and it fits all nine exposure levels rather than withholding the
top one. The method as it was published is
[`examples/published.ini`](../examples/published.ini), and what comes back from it is in
[Reproducing the published run](reproducing.md). For a dataset that is not that one, copy
[`examples/template.ini`](../examples/template.ini): the same ten sections, every switch written
out and commented, with a note on each saying what it means for the experiment.

## See also

- [Troubleshooting](troubleshooting.md): when a stage refuses, and what the summary labels mean
- [Configuration](configuration.md): every key, its default and its meaning
- [Describing your experiment](experiment-design.md): the twelve keys that state facts about the experiment
- [Method notes](method-notes.md): trimming, pooling, collapsing, dropping, and the figures
- [Reproducing the published run](reproducing.md): the four checks against the published outputs
- [README](../README.md): the method, the stages and an example run
