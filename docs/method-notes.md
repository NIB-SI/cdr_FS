# Method notes

Five decisions inside the pipeline that move the numbers, and the reasoning behind the
defaults. The keys named throughout are documented in [Configuration](configuration.md).

## Trimming

The first data-level step, and the only one that changes the values every later stage reads.
Per feature, within each group of `[trim] scope`, values outside
`[lower_percentile, upper_percentile]` are discarded: quality control against objects the
imaging and segmentation measured badly. Two properties follow from *values* being removed
rather than rows, and both surface later:

- **N differs per feature**, so a trimmed table is ragged. That is what
  [dropping features](#dropping-features-from-the-final-table) cleans up afterwards, and what
  `[correlation] fill_missing` has to fill.
- **Nothing is materialised.** Every stage that needs cell-level data reads `[input] table` and
  applies the trim itself, so `cdr-fs run` never calls it; `cdr-fs trim` exists only to write
  the trimmed copy out for inspection.

`[trim] scope` should name one physical unit of the assay, so that the percentiles are taken
across objects that were handled together; in the published run that is one well of one
replicate on one day. The unit has to hold enough objects for a percentile to mean something.
The committed fixture holds one to eight cells per well, so trimming it the published way to
`[p2.5, p97.5]` removes 55% of its values, which is why `examples/quickstart.ini` turns
trimming off. `[correlation] aggregate_by` defaults
to `scope`, on the reasoning that the unit worth trimming within is the unit worth correlating
across; `scope` is read whether or not trimming is on, so leaving it unset is what makes
`aggregate_by` required once `[correlation] enabled` is true.

## Pooling replicates

`per_replicate = false`, the default, merges all replicates at an exposure level into one
distribution: one distance per feature, stratum and contrast, so the curve is fitted to as
many points as there are exposure levels. This is the published behaviour, and it is
defensible here because the upstream two-step row/plate standardization exists precisely to
remove the batch structure that pooling would otherwise bake in.

`per_replicate = true` instead computes one distance per replicate. That puts batch variation
*between* the points rather than inside the distributions, gives each point a spread, and
multiplies the number of points the curve is fitted to, which matters because a five-parameter
model fitted to eight points is at the edge of what AIC/BIC comparison supports. It changes the
numbers, so it no longer reproduces the article.

## Collapsing redundant features

`cdr-fs correlation` aggregates the objects to one median profile per experimental unit
(`[correlation] aggregate_by`, one well of one replicate on one day in the published run),
correlates the features across those units, and clusters them on `1 - |r|` so that a strong
*negative* correlation counts as the redundancy it is. Average linkage, cut at
`1 - [correlation] threshold`; one member of each cluster survives.

**Missing values are filled by default, and it is not a small effect.** Trimming leaves gaps,
so `fill_missing = column_mean` substitutes the feature's overall mean for each missing object
value before aggregating, which is what the published run did. For a feature with many trimmed
values that pulls its unit medians toward the global mean, shrinking the spread the correlation
sees: on the reference dataset the choice moves the cluster count from 99 to 107. `none`
medians whatever values are actually present.

**Average linkage cuts on cluster means, so a member can sit further from its representative
than the threshold suggests.** On the reference dataset the 83 dropped features sit a median
0.052 from the feature that stands for them, `|r|` = 0.95 as advertised, but the loosest sits
at 0.215, which is `|r|` = 0.785. That is what chaining does, and it is why
`correlation_clusters.tsv` records the distance for every member: it is the column to look at
before trusting one feature to speak for a group.

## Dropping features from the final table

Trimming removes values rather than rows, so a feature can be missing for most objects and
still be present as a column. `[drop_missing] enabled` removes the ones too empty to be worth
carrying: a feature missing `[drop_missing] max_missing` percent of the table **or more** is
left out, 30% by default.

It is applied over the **whole** table, before anything downstream subsamples: one decision for
the whole experiment, so a feature is either in the analysis or out of it. Deciding per
subsample instead lets the same feature be present in one day's file and absent from another.

On the reference dataset the default takes 97 features to 95, dropping two that are missing
47% and 95% of their values. The next thinnest surviving feature is missing 13%, and the rule
drops at the threshold or above, so every threshold in (13%, 47%] gives the same 95: a wide
plateau rather than a knife edge. (The same rule applied to the *published* metadata split
takes 99 to 97, dropping the same two features; see
[Reproducing the published run](reproducing.md) for why there are two chains.)

`[drop_missing] exclude` is the other half: exact feature names to leave out whatever their
quality. Exact names rather than patterns, because a regex that quietly takes a second feature
with it is the wrong tool for a decision made one feature at a time; a name that matches
nothing is reported by `cdr-fs check`. `final_<list>_features.tsv` records which rule removed
each feature, and an explicit exclusion is reported as one even when the missing-data rule
would also have caught it.

Nothing else is filtered: a constant feature, or one whose surviving values all come from one
object, is named in the report and flagged in that table rather than dropped.

## The figures

`cdr-fs plot` draws three things, from the tables the other stages wrote and never by
recomputing anything:

- **Fit panels**: the distance points, the six fitted curves and a legend ordered by
  information criterion, one panel per feature and stratum. This is the figure the article's
  Figure 4 was composed from. `--grid N` sets panels per row and column and `--features FILE`
  restricts which features are drawn; without it a full run is several hundred pages.
- **Distance distributions**: every distance in a table, one column of points per feature,
  features ordered by median, on a broken axis (linear below a split, logarithmic above). On
  `emd_baseline.tsv` this is the between-replicate reproducibility floor; on `emd.tsv` it is the
  treatment distances to read against that floor.
- **Dendrogram**: the tree correlation collapsing cut, with the cut drawn on it and each cluster
  of three or more members coloured from the cluster table.

A curve is the fitted function evaluated at the exposure levels and joined up, as the published
figures show it, and `fit.tsv` carries the parameters at full precision: rounded to six digits
a fitted inflection can settle against the `log(x + 1e-10)` guard and the drawn curve stops
matching the AIC printed beside it, which happens to about one fit in sixteen.

## See also

- [Quickstart](quickstart.md): an annotated first run, and where its numbers come from
- [Troubleshooting](troubleshooting.md): exit codes, the summary labels, and where to look when a run retains nothing
- [Configuration](configuration.md): every key, its default and its meaning
- [Describing your experiment](experiment-design.md): the twelve keys that state facts about the experiment
- [Reproducing the published run](reproducing.md): the four checks against the published outputs
- [README](../README.md): the method, the stages and an example run
