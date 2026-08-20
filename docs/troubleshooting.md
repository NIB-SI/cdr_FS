# Troubleshooting

What to do when a run does not do what you expected: which stages refuse and why, what the
`run` summary's labels mean, and the mistakes a first run on new data actually makes. For a run
that went right, see [Quickstart](quickstart.md).

## Exit codes, and which stages refuse

`0` success, `2` a configuration error, `3` nothing to do.

A stage refuses rather than write an artefact that would read as a result. `fit` will not write
a table when no series was complete, because `select` would read it and retain nothing, which
looks like "no feature responds" rather than "nothing was fitted". `correlation` and
`drop_missing` will not run on an empty feature list; the final table would otherwise hold
metadata and no measurements, which also looks like a result. Called by name, `trim`, `select`
and `correlation` exit 3 when their own `enabled` is false.

Inside `run` a refusal stops the chain: the stage's message is printed, the summary marks the
rest `not reached`, and the run exits 3. A stage a switch turned off is not a refusal. It is a
declared outcome, named in the run's header, and the run still exits 0.

## What the summary labels mean

| Label | What it means |
|---|---|
| `ok` | ran, and wrote its output |
| `off` | a switch. `[select] enabled` takes `emd`, `fit` and `select` together; `[correlation] enabled` takes `correlation`; `plot` is off when the run writes no table any figure could be drawn from |
| `not run` | only ever `trim`, which is never part of a run: every stage trims from `[input] table` itself, so `cdr-fs trim` writes an inspection copy nothing reads |
| `no figures` | `plot` drew nothing, but the table is complete, so the run still exits 0. Reachable with the selection off and one feature reaching `correlation`, where the only figure on the plan is a dendrogram with nothing to cut |
| `stopped` | this stage refused, and the run exits 3 |
| `not reached` | the chain stopped before this stage |

## The `[columns]` feature count is not the one you expected

`[schema] metadata_patterns` matched too much, and nothing downstream can catch it: the split is
by inversion, so whatever the patterns fail to match becomes a feature. On the quickstart
fixture, replacing the three explicit `^counts_...$` patterns with a tidy `^counts_` turns
`30 = 10 metadata + 20 feature(s)` into `30 = 12 metadata + 18 feature(s)` and drops the
`counts_*` group from the features altogether, taking with it a feature the published run
retained. Read that line from `cdr-fs check` before anything else;
[Configuration](configuration.md#one-warning-worth-reading-before-you-write-metadata_patterns)
has both directions of the mistake.

## The run retained nothing

`select_evidence.tsv` first. It carries the linear slope and the winning model for every feature
and stratum, so it says which half of the retention rule the features failed; `select`'s own
report gives the same split as two counts. When `select` retains nothing, `correlation` refuses
on the empty list and the run stops there, so there is no final table to look at.

An empty result is not evidence that the exposure axis is right. A `[design] levels` list
written high to low negates every slope and still produces a short, plausible list: on the
quickstart fixture the correct axis retains nine features and the reversed one four. See
[Describing your experiment](experiment-design.md#when-a-piece-of-the-design-is-missing).

## `plot` skipped a figure

Each figure it could not draw is reported on stderr as `skipped <figure> - <reason>`, and the
rest are still drawn, so one bad table does not lose the others. Drawing nothing at all prints
`error: nothing to draw` and exits 3. Four reasons name a command to run:

```
skipped fits - needs emd and fit (run `cdr-fs fit`)
skipped emd - needs emd.tsv (run `cdr-fs emd`)
skipped baseline - needs emd_baseline.tsv (run `cdr-fs emd`)
skipped dendrogram - needs correlation_linkage (run `cdr-fs correlation`)
```

Two name a configuration setting instead, and running something will not fix them. They report a
figure this run will never have, so naming a command would be wrong:

```
skipped baseline - [emd] baseline is none, so no between-replicate distances were computed
skipped dendrogram - [correlation] enabled is false, so there is no tree to draw
```

Anything else is the figure function's own refusal: an empty feature list, a distance table with
no rows or none of them finite, a linkage table too short to be a tree. What each figure shows is
under [Method notes](method-notes.md#the-figures).

## A stratum came back as `NaN`

Reading a stage table into an analysis of your own loses an unstratified run's stratum label,
which is the empty string: `pandas.read_csv` reads an empty field as `NaN`, and `groupby` then
drops it. Use `cdr_fs.schema.read_stage_table`. The mechanism is at the end of
[Describing your experiment](experiment-design.md#what-the-configuration-cannot-say).

## Where a stage lives in the code

One module per stage, named for it: `emd.py`, `fit.py`, `select.py`, `correlation.py`,
`drop_missing.py`, `trim.py` and `plots.py`, with `config.py` and `schema.py` underneath.

## See also

- [Quickstart](quickstart.md): an annotated first run, and what its numbers are made of
- [Configuration](configuration.md): every key, its default and its meaning
- [Describing your experiment](experiment-design.md): the twelve keys that state facts about the experiment
- [Method notes](method-notes.md): trimming, pooling, collapsing, dropping, and the figures
- [Reproducing the published run](reproducing.md): the four checks against the published outputs
- [README](../README.md): the method, the stages and an example run
