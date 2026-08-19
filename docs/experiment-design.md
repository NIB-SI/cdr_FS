# Describing your experiment

Which configuration keys state a fact about the **experiment** rather than a choice about
the method, and what to write when the design you have is not the one this was built on — no
time course, no replicates, fewer exposure levels. Read this before the first run on a new
dataset. The last section is the one to read second: what the configuration cannot express,
including the parts that fail quietly rather than loudly.

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
| `[correlation] aggregate_by` | the columns identifying one experimental unit for correlation |

`[fit] models` sits between the two: it is a method choice, but the design constrains it,
because a curve needs more points than it has parameters.

## When a piece of the design is missing

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
| **no need to trim** | `[trim] enabled = false` | Nothing, except that `[correlation] aggregate_by` no longer has `[trim] scope` to default to, so set it explicitly |
| **no wish to collapse correlated features** | `[correlation] enabled = false` | The redundancy collapse. `cdr-fs missing_data` then applies `selected.txt` instead of `representatives.txt`, and says which it used |

Two of these deserve spelling out.

**Level labels are opaque strings, and nothing sorts them.** The order in `[design] levels`
is the exposure axis, whatever the labels look like. The reference dataset's own labels run
`11, 10, … 2` with `11` the control and `2` the top dose, so its `levels` line reads
`10,9,8,…,2` — descending as text, ascending in exposure. That is also the quietest way to
get a wrong answer: a `levels` list written high to low runs to completion, produces an
identical distance table, negates every linear slope, and retains nothing.

`[design] dose` guards half of that. When you supply it, the doses must rise along `levels`,
and a vector that falls while the levels are declared to rise is refused — so the two lists
cannot silently disagree. **Supply the doses even when you fit on `rank`**, for that reason
alone. What it cannot catch is both lists being wrong the same way: reverse `levels` and
leave `dose` untouched and the pairing is broken but still monotone, so the check passes.
Nothing can catch that, because only you know which label was your top dose.

The backstop is therefore a line of `cdr-fs check`, which prints the declaration in words:

```
exposure axis     10 is the LOWEST exposure, 2 the highest   (11.3683 -> 1000)
```

Read it against your own plate map. It is the one place the tool states what it believes
about the direction of your experiment.

**`[correlation] aggregate_by` must name a unit that varies along the exposure axis.**
Correlations are computed between unit medians, so if each unit spans the whole dilution
series then every unit looks alike and the stage collapses nothing. In the published plate
layout one well holds one concentration, which is why `Metadata_Day,Metadata_Biorep,
Metadata_Well` works there; if your wells each hold a whole series, add the level column.
The report line `|r| >= 0.9 on median profiles over N unit(s) of …` is where you check that
N is what you expected.

## The smallest configuration that runs

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
[correlation]
enabled = false
[missing_data]
drop_missing = false
[output]
dir = /PATH/TO/results
```

Six exposure levels, because that is what all six models need; one metadata column, one
control, everything else defaulted. `cdr-fs check` lists the twenty-odd keys it filled in.

## What the configuration cannot say

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

And one that bites only afterwards, when you read the stage tables back into an analysis of
your own: **an experiment with no `[schema] group_by` has one stratum, and its label is the
empty string.** Written to a delimited file that is an empty field, and `pandas.read_csv` reads
an empty field as `NaN` whatever dtype it is given; `groupby` then drops it, taking the whole
table with it. Read those tables with `cdr_fs.schema.read_stage_table`, which fills the field
back to the empty string it was, as every stage already does.

## See also

- [Configuration](configuration.md) — every key, its default and its meaning
- [Method notes](method-notes.md) — pooling, collapsing, dropping, and the figures
- [Reproducing the published run](reproducing.md) — the four checks against the published outputs
- [README](../README.md) — the method, the stages and an example run
