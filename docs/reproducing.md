# Reproducing the published run

What comes back when the tool is run on the dataset the method was published on, and how the
four checks are arranged so that a difference is attributable to the one stage under test.
Read this to see how far the reproduction goes, which numbers belong to the published
metadata split and which do not, and where the route through the last step deliberately
departs from the published one.

`examples/published.ini` is the configuration for the dataset the method was published on.
Four published outputs are checked against, each stage given the published input to the stage
before it rather than this tool's own output, so that a difference is attributable to the
stage under test:

| Stage | Published output | Result |
|---|---|---|
| `emd` | the two EMD tables — 16,946 treatment and 11,292 baseline distances | both population sizes on every one of the 28,238 rows exact, distances agreeing to 8.5e-13 relative, against an asserted tolerance of 1e-9 |
| `select` | the two retained feature lists — **182** across all days, **374** on D5 alone | both reproduced as identical sets |
| `prune` | the all-days list after pruning — **99** features | reproduced, and its composition matches the published categorization across all 4 organelles x 7 measurement families |
| `subset` | the final retained list — **95** features | reproduced: 94 `rp_norm_*` plus `counts_RelateLysoCell` |

The selection gate runs off the published fit table, so it isolates the retention rule from
the curve fitting; `tests/test_golden_selection.py` needs only a 1 MB committed fixture and
runs in a second. The others need the large inputs and are opt-in — see [Tests](testing.md).

Three notes on what those numbers are and are not.

**The pruning check is on composition, not just on the count.** A single wrong merge would
still give a plausible number, so the assertion is against the published feature
categorization, every cell of which is labelled in
[HCS-proc's figure](https://github.com/NIB-SI/HCS-proc). Reproducing the count *and* the
composition takes the same 99 features.

**The 182 and the 99 belong to the published metadata split, and the 95 does not.** The
published run treated the eight `Number_Object_Number` columns as measured features; they are
CellProfiler's within-image object index, so `examples/published.ini` declares them metadata —
see the warning [under Configuration](configuration.md#one-warning-worth-reading-before-you-write-metadata_patterns).
The two gates above are asserted against the published split, because reproducing a published
intermediate means reproducing the choices that produced it. With the object indices out, the
same rules give **463 features → 175 after selection → 97 after pruning → 95**.

**The last step of the published route moved, and it lands in the same place.** The pipeline
dropped features with 30% or more missing values while building the subsets for its
dimension-reduction analyses: *per day-subset*, downstream of the selection, along with two
features excluded by hand. `cdr-fs subset` applies the same 30% rule once, over the whole table,
before anything subsamples — so a feature is either in the analysis or out of it, where a
per-subsample decision lets the same feature be present in one day's file and absent from
another. That takes 97 to **95**, and needs neither hand exclusion: with the object indices
classified correctly, both of those features survive on their own merits. The mechanism for
excluding a feature by name is there, in `[subset] exclude`, and ships empty.

## An aside that is really a check

The published run produced two selections, one per gate: across all days, and on D5 alone.
Neither list contains the other - `counts_RelateLysoCell` is retained across all days but not
on D5, and `counts_RelateMitoCell` the reverse. Only the hybrid retention rule can do that. A
positive slope is required on **any** stratum, so the all-days gate has four chances where D5
has one; if both tests used `all`, the D5 set would necessarily contain the all-days set. The
asymmetry is therefore visible in the published output itself, not merely in a reading of the
code.

## See also

- [Configuration](configuration.md) — every key, its default and its meaning
- [Describing your experiment](experiment-design.md) — the twelve keys that state facts about the experiment
- [Method notes](method-notes.md) — pooling, collapsing, dropping, and the figures
- [Tests](testing.md) — what the suite protects, and how to run the golden regressions
- [README](../README.md) — the method, the stages and an example run
