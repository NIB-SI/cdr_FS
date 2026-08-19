# Tests

What the test suite protects, promise by promise, and how to run the two golden regressions
that need the large data files. The suite is deliberately small, so this page is mostly a
list of which claim each group of tests is holding down. Read it before adding a test, or
when one fails and you need to know what it was pinning.

```bash
pytest
```

Forty-six tests, about twenty seconds, no large data needed. The suite is deliberately
small: every test in it either backs one of the four numbers in
[Reproducing the published run](reproducing.md), or pins a boundary that
a one-character edit would move, or is the regression for a bug that actually happened.
Tests that re-checked what another test already covered, or that asserted an implementation
detail rather than a promise, were removed rather than carried into the release.

What that leaves, by what it protects:

- **The reproduction gates**, four of them, described in
  [Reproducing the published run](reproducing.md).
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

## The design comes from the configuration

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

## The golden regressions

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

## See also

- [Configuration](configuration.md) — every key, its default and its meaning
- [Describing your experiment](experiment-design.md) — the twelve keys that state facts about the experiment
- [Method notes](method-notes.md) — pooling, collapsing, dropping, and the figures
- [Reproducing the published run](reproducing.md) — the four checks against the published outputs
- [README](../README.md) — the method, the stages and an example run
