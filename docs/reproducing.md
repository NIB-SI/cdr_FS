# Reproducing the published run

What comes back when the tool is run on the dataset the method was published on.
`examples/published.ini` is that configuration. Four published outputs are checked against,
each stage given the published input to the stage before it rather than this tool's own output,
so that a difference is attributable to the stage under test:

| Stage | Published output | Result |
|---|---|---|
| `emd` | the two EMD tables, 16,946 treatment and 11,292 baseline distances | both population sizes exact on every one of the 28,238 rows, distances agreeing to 8.5e-13 relative against an asserted tolerance of 1e-9 |
| `select` | the two retained feature lists, **182** across all days and **374** on D5 alone | both reproduced as identical sets |
| `correlation` | the all-days list after the correlation stage, **99** features | reproduced, and its composition matches the published categorization across all 4 organelles x 7 measurement families |
| `drop_missing` | the final retained list, **95** features | reproduced: 94 `rp_norm_*` plus `counts_RelateLysoCell` |

The selection gate runs off the published fit table, so it isolates the retention rule from the
curve fitting: `tests/test_golden_selection.py` needs only a 1 MB committed fixture. The others
read the Zenodo files from `data/` and are opt-in, because reading gigabytes takes minutes:

```bash
pytest                      # everything that needs no data
CDR_FS_GOLDEN=1 pytest      # and the gates that read the published tables
```

The suite is deliberately small: it holds these four numbers and the three shipped
configurations.

**The correlation check is on composition, not just on the count.** A single wrong merge would
still give a plausible number, so the assertion is against the published feature categorization,
every cell of which is labelled in [HCS-proc's figure](https://github.com/NIB-SI/HCS-proc).

**The 182 and the 99 belong to the published metadata split, and the 95 does not.** The
published run treated the eight `Number_Object_Number` columns as measured features; they are
CellProfiler's within-image object index, so `examples/published.ini` declares them metadata
(see the warning [under Configuration](configuration.md#one-warning-worth-reading-before-you-write-metadata_patterns)).
Those two gates are asserted against the published split, because reproducing an intermediate
means reproducing the choices that produced it. With the object indices out, the same rules give
**463 features → 175 after selection → 97 after the correlation stage → 95**.

**The last step of the published route moved, and it lands in the same place.** The pipeline
dropped features with 30% or more missing values while building the subsets for its
dimension-reduction analyses: *per day-subset*, downstream of the selection, along with two
features excluded by hand. `cdr-fs drop_missing` applies the same 30% rule once over the whole
table instead (see [Method notes](method-notes.md#dropping-features-from-the-final-table)),
which takes 97 to **95** and needs neither hand exclusion: with the object indices classified
correctly, both of those features survive on their own merits.

## Reading the reference configuration

Three of `examples/published.ini`'s values need a sentence that does not fit in the file.

**The dose vector is a 1.75-fold serial dilution from 1000 mg/L**, `1000 / 1.75^k` for
k = 8…0, written to six decimals. The article's SI prints the same series truncated
(1000, 571.42, 326.530, 186.588 … 11.36), which is how the factor was pinned to exactly 1.75
rather than the 1.7502 a back-calculation from those truncated endpoints suggests. The run
fits on `rank`, not on these values, but they are what makes the direction of the exposure
axis checkable on a dataset whose labels descend as its doses climb.

**The highest exposure is withheld from curve fitting.** `exclude_from_fit = 2` keeps
concentration–response detection inside the sub-cytotoxic range the article is about. The
distance is still computed for that level and still appears in `emd.tsv`; only the fit skips
it, which is what leaves the eight fitted contrasts 11v10, 11v9 … 11v3.

**`[drop_missing] exclude` ships empty, and that is a result rather than an omission.** The
published run struck out two features by hand at this step,
`rp_norm_AreaShape_Compactness_RelateLysoCell` and
`rp_norm_Texture_AngularSecondMoment_GrayLys_3_00_256_RelateLysoCell`. With the object indices
declared metadata the rules reach the published feature count without them and both survive on
their own merits; naming them would take the final list from 95 to 93.

## An aside that is really a check

The published run produced two selections, one per gate: across all days, and on D5 alone.
Neither list contains the other. `counts_RelateLysoCell` is retained across all days but not on
D5, and `counts_RelateMitoCell` the reverse, and only the hybrid retention rule can do that: a
positive slope is required on **any** stratum, so the all-days gate has four chances where D5
has one, and if both tests used `all` the D5 set would necessarily contain the all-days set.
The asymmetry is visible in the published output itself, not merely in a reading of the code.

## See also

- [Configuration](configuration.md): every key, its default and its meaning
- [Describing your experiment](experiment-design.md): the twelve keys that state facts about the experiment
- [Method notes](method-notes.md): trimming, pooling, collapsing, dropping, and the figures
- [README](../README.md): the method, the stages and an example run
