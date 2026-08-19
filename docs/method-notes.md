# Method notes

Four decisions inside the pipeline that move the numbers, and the reasoning behind the
defaults: how replicates are pooled before a distance is computed, how near-redundant
features are collapsed to one representative, which features are dropped from the final
table, and what the three diagnostic figures show. Read these when you are deciding whether
a default suits your data, or when an output is not the one you expected. The keys named
throughout are documented in [Configuration](configuration.md).

## Pooling replicates

When computing the distance between control and exposure level *C*, what counts as "the
population"? `per_replicate = false`, the default, merges all replicates at that level into
one distribution: one distance per feature, stratum and contrast, so the curve is fitted to
as many points as there are exposure levels. This is the published behaviour, and it is
defensible here because the upstream two-step row/plate standardization exists precisely to
remove the batch structure that pooling would otherwise bake in.

`per_replicate = true` instead computes one distance per replicate. That puts batch
variation *between* the points rather than inside the distributions, gives each point a
spread, and multiplies the number of points the curve is fitted to — which matters, because
a five-parameter model fitted to eight points is at the edge of what AIC/BIC comparison
supports. It also changes the numbers, so it no longer reproduces the article. Both are
available; the default reproduces the paper.

## Collapsing redundant features

`cdr-fs prune` aggregates the objects to one median profile per experimental unit -
`[prune] aggregate_by`, one well of one replicate on one day in the published run - correlates
the features across those units, and clusters them on `1 - |r|` so that a strong *negative*
correlation counts as the redundancy it is. Average linkage, cut at `1 - [prune] threshold`.
One member of each cluster survives.

Two things about it are worth knowing before trusting the output.

**Missing values are filled by default, and it is not a small effect.** Trimming leaves gaps,
so `fill_missing = column_mean` substitutes the feature's overall mean for each missing object
value before aggregating - which is what the published run did. For a feature with many
trimmed values that pulls its unit medians toward the global mean and shrinks the spread the
correlation sees. On the reference dataset the choice moves the cluster count by eight, from
99 to 107. `none` medians whatever values are actually present.

**Average linkage cuts on cluster means, so a member can sit further from its representative
than the threshold suggests.** On the reference dataset the 83 dropped features sit a median
0.052 from the feature that stands for them — `|r|` = 0.95, as advertised — but the loosest
sits at 0.215, which is `|r|` = 0.785. That is what chaining does, and it is why
`prune_clusters.tsv` records the distance for every member rather than only the cluster
number: it is the column to look at before trusting one feature to speak for a group.

## Dropping features from the final table

Trimming removes values rather than rows, so a feature can be missing for most objects and
still be present as a column. `[subset] drop_missing` removes the ones too empty to be worth
carrying: a feature missing `[subset] max_missing` percent of the table **or more** is left out,
30% by default.

It is applied over the **whole** table, before anything downstream subsamples. That is one
decision for the whole experiment, so a feature is either in the analysis or out of it —
whereas filtering each subsample separately lets the same feature be present in one day's file
and absent from another, which makes the resulting sets incomparable.

On the reference dataset the default takes 97 features to 95, dropping two that are missing
47% and 95% of their values. The next thinnest surviving feature is missing 13%, and the rule
drops at the threshold or above, so every threshold in (13%, 47%] gives the same 95 — it sits
on a wide plateau rather than on a knife edge, which is the useful thing to know about
choosing one.

(The same rule applied to the *published* metadata split takes 99 to 97, dropping the same two
features. Two chains, two numbers; see [Reproducing the published run](reproducing.md) for why
there are two.)

`[subset] exclude` is the other half: exact feature names to leave out whatever their quality —
the escape hatch for a judgement no rule expresses, such as a feature known to be an artifact of
one assay. Exact names rather than patterns, because a regex that quietly takes a second feature
with it is the wrong tool for a decision made one feature at a time; a name that matches nothing
is reported by `cdr-fs check`, since the failure mode is that the feature stays in and nobody
notices. `subset_<list>_features.tsv` records which rule removed each feature, and an explicit
exclusion is reported as an exclusion even when the missing-data rule would also have caught it.

Nothing else is filtered. A constant feature, or one whose surviving values all come from a
single object, is named in the report and flagged in `subset_<list>_features.tsv`; what to do
about that depends on the analysis, and guessing is worse than saying so.

## The figures

`cdr-fs plot` draws three things, from the tables the other stages wrote and never by
recomputing anything:

- **Fit panels** — one panel per feature and stratum, with the distance points, the six fitted
  curves and a legend ordered by information criterion. This is the figure the article's
  Figure 4 was composed from. `--grid N` sets panels per row and column, `--features FILE`
  restricts which features are drawn; without it a full run is several hundred pages.
- **Distance distributions** — every distance in a table, one column of points per feature,
  features ordered by median, on a broken axis (linear below a split, logarithmic above). On
  `emd_baseline.tsv` this is the between-replicate reproducibility floor; on `emd.tsv` it is
  the treatment distances to read against that floor.
- **Dendrogram** — the tree correlation pruning cut, with the cut drawn on it and each cluster
  of three or more members coloured.

A curve in a fit panel is the fitted function evaluated at the exposure levels and joined up,
as the published figures show it. The `fit.tsv` parameters are written at full precision so a
drawn curve reproduces the AIC printed beside it — which is asserted, and which fails for about
one fit in sixteen if the parameters are rounded to six digits, because a fitted inflection can
settle right against the `log(x + 1e-10)` guard.

The dendrogram colours links and leaf labels from the same cluster table, so a coloured subtree
is a cluster. The original drew the tree with one linkage method and cut it with another, so
its colours were only approximately the clustering.

## See also

- [Configuration](configuration.md) — every key, its default and its meaning
- [Describing your experiment](experiment-design.md) — the twelve keys that state facts about the experiment
- [Reproducing the published run](reproducing.md) — the four checks against the published outputs
- [Tests](testing.md) — what the suite protects, and how to run the golden regressions
- [README](../README.md) — the method, the stages and an example run
