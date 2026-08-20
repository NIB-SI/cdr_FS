# Configuration

Every key a `.ini` file can set, what it defaults to and what it means. This is the page to
keep open while writing a configuration; it assumes you have already copied
[`examples/template.ini`](../examples/template.ini) — the annotated starting point, which
carries every section and every switch — and want to know what a particular key does. Two of
its sections run longer than a reference entry because the obvious reading of them is wrong:
the one on `metadata_patterns`, and the one on how the exposure levels are spaced along the
fitted axis.

`.ini` carries no type information, so every value is parsed and cross-checked before any
data is read. A wrong configuration fails in the first second with a message naming the
section, the key and — where there is one to name — the fix, rather than producing
plausible-looking wrong numbers. All ten sections must be present even when every key in them
is left at its default, so that a mistyped section name is an error instead of a silent
fall-back; unknown keys are rejected for the same reason.

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
| | `dose` | *(empty)* | actual doses, index-matched to `levels`. Required for `x_scale = dose`, and worth supplying regardless: it is the only check on the direction of the exposure axis (see [Experiment design](experiment-design.md)) |
| | `exclude_from_fit` | *(empty)* | levels withheld from fitting but still measured |
| `trim` | `enabled` | *required* | trimming changes the data, so it must be stated |
| | `lower_percentile`, `upper_percentile` | *required when enabled* | the interval kept is closed |
| | `scope` | *required when enabled* | columns whose groups the percentiles are computed within |
| `emd` | `contrasts` | `control_vs_each` | or an explicit list, `11v10,11v9,...` |
| | `baseline` | `control_across_replicates` | the between-replicate reproducibility floor, or `none` |
| | `per_replicate` | `false` | see [Pooling replicates](method-notes.md#pooling-replicates) |
| `fit` | `models` | all six | any subset of `BC4,BC5,LL4,WB1.4,Lin,Con` |
| | `x_scale` | `rank` | `rank` (position in the series) or `dose` |
| | `rank_by` | `aic_plus_bic` | `aic_plus_bic`, `aic` or `bic` |
| `select` | `slope_positive` | `any` | positive linear slope on `any` or `all` strata |
| | `nonconstant` | `all` | constant model beaten on `any` or `all` strata |
| | `strata` | *(all present)* | which strata the rule is applied over |
| `correlation` | `enabled` | *required* | |
| | `threshold` | `0.9` | on \|r\|; features are clustered by `1 - |r|` |
| | `linkage` | `average` | `average`, `complete` or `single` |
| | `representative` | `alphabetical` | which member of a cluster to keep — `alphabetical` or `first` |
| | `aggregate_by` | *(`[trim] scope`)* | columns identifying one experimental unit; correlations are computed between unit medians. Explicitly empty correlates the object rows as they are |
| | `fill_missing` | `column_mean` | `column_mean` or `none` — see [Collapsing redundant features](method-notes.md#collapsing-redundant-features) |
| `drop_missing` | `enabled` | *required* | drop features too empty to carry; a filter on the output must be stated |
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

## One warning worth reading before you write `metadata_patterns`

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

Mechanically: one regex per line, so a pattern may itself contain commas, and patterns are
**searched, not fullmatched**. `counts_` alone therefore matches every column with that text
anywhere in it; anchor with `^...$` to name a single column.

The mistake runs the other way too, and it is harder to see. A column can be named like a
measurement and still be a label:

```
rp_norm_Number_Object_Number_*      CellProfiler's object index -> metadata
```

`Number_Object_Number` is the running number of an object within its image. It has no
biological content, but it *does* respond to exposure — the number of objects per image changes
with dose, which shifts the index distribution, which moves the earth mover's distance. Seven of
its eight columns passed the concentration-response gate in the published run and two survived
correlation collapsing; the pipeline removed them by name at the very end. Declaring them metadata
is the same decision made where it belongs, and `examples/published.ini` does. Nothing in the
tool guesses this for you: only a pattern you write keeps a label out of the analysis.

## Spacing the exposure axis

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

## See also

- [Describing your experiment](experiment-design.md) — which of these keys state facts about the experiment
- [Method notes](method-notes.md) — trimming, pooling, collapsing, dropping, and the figures
- [Reproducing the published run](reproducing.md) — the four checks against the published outputs
- [README](../README.md) — the method, the stages and an example run
