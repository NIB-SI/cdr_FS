# Configuration

Every key a `.ini` file can set, what it defaults to and what it means. Copy
[`examples/template.ini`](../examples/template.ini), the annotated starting point that carries
every section and every switch, and use this page to look up what a particular key does.

`.ini` carries no type information, so every value is parsed and cross-checked before any data
is read: a wrong configuration fails in the first second, with a message naming the section,
the key and the fix. All ten sections must be present even when every key in them is left at
its default, so that a mistyped section name is an error rather than a silent fall-back, and
unknown keys are rejected for the same reason.

| Section | Key | Default | Meaning |
|---|---|---|---|
| `input` | `table` | *required* | path to the per-object table |
| | `sep` | `tab` | `tab`, `comma` or `semicolon`, as a keyword: `.ini` has no escapes, so `\t` would arrive as a literal backslash-t |
| `schema` | `metadata_patterns` | *required* | one regex per line; **everything not matched is a feature** |
| | `condition` | *required* | column holding the exposure level |
| | `group_by` | *(empty)* | stratification, e.g. time point; empty means one stratum |
| | `pool_over` | *(empty)* | replicate column merged into one distribution |
| `design` | `control` | *required* | the control level's label |
| | `levels` | *required* | exposure levels, ordered **low to high**; this is the response axis |
| | `dose` | *(empty)* | actual doses, index-matched to `levels`. Required for `x_scale = dose`, and the only automatic check on the direction of the exposure axis, so worth supplying regardless (see [Describing your experiment](experiment-design.md)) |
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
| `select` | `enabled` | *required* | `false` skips `emd`, `fit` and `select`, carrying every feature into the filtering stages |
| | `slope_positive` | `any` | positive linear slope on `any` or `all` strata |
| | `nonconstant` | `all` | constant model beaten on `any` or `all` strata |
| | `strata` | *(all present)* | which strata the rule is applied over |
| `correlation` | `enabled` | *required* | |
| | `threshold` | `0.9` | on \|r\|; features are clustered by `1 - |r|` |
| | `linkage` | `average` | `average`, `complete` or `single` |
| | `representative` | `alphabetical` | which member of a cluster to keep, `alphabetical` or `first` |
| | `aggregate_by` | *(`[trim] scope`)* | columns identifying one experimental unit; correlations are computed between unit medians. Explicitly empty correlates the object rows as they are |
| | `fill_missing` | `column_mean` | `column_mean` or `none`; see [Collapsing redundant features](method-notes.md#collapsing-redundant-features) |
| `drop_missing` | `enabled` | *required* | drop features too empty to carry; a filter on the output must be stated |
| | `max_missing` | `30` | percent; a feature missing this much of the table **or more** is dropped |
| | `exclude` | *(empty)* | exact feature names to leave out of the final table, whatever else says |
| `output` | `dir` | *required* | where results are written |

Two things this configuration makes explicit that the original scripts left implicit. **The
metadata list** lived as a literal in five scripts, in two versions; now one drives every
stage. **The two halves of the retention rule** use *different* quantifiers: a positive slope
on **any** stratum, the constant model beaten on **all** of them, and separate keys let you
choose `all`/`all` for the stricter rule.

## One warning worth reading before you write `metadata_patterns`

Patterns name the *metadata*, so anything you fail to match becomes a feature and anything you
match too broadly stops being one. A new measured channel therefore needs no configuration
change; an over-broad pattern silently removes a real measurement, and nothing later in the
run says so. In the reference dataset five columns begin with `counts_`, and only three of
them are metadata:

```
counts_Cells, counts_Cytoplasm, counts_FilteredNuclei     segmentation QC  -> metadata
counts_RelateLysoCell, counts_RelateMitoCell              organelle counts -> FEATURES
```

The last two are the per-cell lysosomal and mitochondrial counts the article appends to each
cell's profile, and `counts_RelateLysoCell` survives the entire selection into the published
retained list, so a tidy-looking `^counts_` would discard a feature the published run kept.
`cdr-fs check` prints the resolved split grouped by leading name token, so `counts_* 2` is
visible next to `rp_* 461`, and warns about any pattern that matched nothing. One regex per
line, so a pattern may itself contain commas; patterns are **searched, not fullmatched**, so
anchor with `^...$` to name a single column.

The mistake runs the other way too, and it is harder to see: a column can be named like a
measurement and still be a label. `rp_norm_Number_Object_Number_*` is CellProfiler's running
number of an object within its image. It has no biological content, but it *does* respond to
exposure, because the number of objects per image changes with dose, which shifts the index
distribution and moves the earth mover's distance. Seven of its eight columns passed the
published run's gate and two survived correlation collapsing, to be removed by name at the
very end; `examples/published.ini` declares them metadata instead, which is the same decision
made where it belongs. Only a pattern you write keeps a label out of the analysis.

## Spacing the exposure axis

`[fit] x_scale` decides where along the x-axis the exposure levels sit. `rank` puts them at
0, 1, 2, ... - evenly spread, which is what the published run did; `dose` puts them at the
values in `[design] dose`. The catch is that the four sigmoid models evaluate `log(x)`
themselves, so the axis they see is not the axis you supplied:

```
x_scale = rank  ->  log(x) = [-23.03, 0, 0.69, 1.10, 1.39, 1.61, 1.79, 1.95]
x_scale = dose  ->  log(x) = [  2.43, 2.99, 3.55, 4.11, 4.67, 5.23, 5.79, 6.35]
                             spacing 0.559616 throughout - exactly even
```

So **`dose` is what spaces the levels evenly in log-concentration** and `rank` is not: it
flings the lowest exposure out to `log(1e-10)`. Neither choice dominates, because the two
halves of the retention rule want different things:

| the series really is... | `rank` | `dose` |
|---|---|---|
| a logistic in log-dose (textbook DRC) | LL4 AIC -41 | LL4 AIC **-298** |
| a straight line in log-dose | Lin AIC **-573** | Lin AIC -19 |

`rank` suits the **linear** model, whose slope sign is what the retention rule tests, because
for a geometric series the rank *is* log-dose up to an affine shift; `dose` suits the
**sigmoids**, and what changes most is how easily one of them beats the constant model. The
slope test is not strictly invariant either: on arbitrary eight-point series the two axes
disagree about the sign one time in seven, though for a series that rises monotonically with
exposure the sign is positive on either axis.

`rank` stays the default because it is what was published. Switch to `dose` when the exposure
series is not geometrically spaced, since rank is then not a transform of log-dose at all.
There is deliberately no `log_dose`: it would hand `log(dose)` to models that take a logarithm
themselves, and it is NaN for any dose below 1.

## See also

- [Describing your experiment](experiment-design.md): which of these keys state facts about the experiment
- [Method notes](method-notes.md): trimming, pooling, collapsing, dropping, and the figures
- [Reproducing the published run](reproducing.md): the four checks against the published outputs
- [README](../README.md): the method, the stages and an example run
