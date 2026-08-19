# Test fixtures

Small, committed slices of the published dataset, so the suite runs in CI without the
3.7 GB table. Neither file is a copy: `columns_published.txt` is a header, `subset.tsv`
is a row and column subset.

Provenance for both: `cell_ID_pooled_median_row_plate_standardization_cid.txt` from
<https://doi.org/10.5281/zenodo.17951792> (CC-BY-4.0) — the untrimmed, row/plate
standardized per-cell table, 503,920 cells x 481 columns.

## `columns_published.txt` (24 KB)

The 481 column names of the full table, one per line, in table order:

```bash
head -1 cell_ID_pooled_median_row_plate_standardization_cid.txt | tr '\t' '\n'
```

This is what lets `test_schema.py` assert the real arithmetic — 481 columns, 10 metadata,
471 features, 469 `rp_norm_*` plus the two organelle counts — with no data present. It is
also the guard on the `counts_` trap: a `^counts_` metadata pattern resolves to 469
features here, and the test says so.

## `subset.tsv` (440 KB)

1,272 cells x 30 columns: every metadata column, both organelle-count features, and 18
`rp_norm_*` features sampled evenly across the feature list. Up to 8 cells per
(concentration, day, replicate) combination, taken in file order.

The design is preserved, including its one hole: concentration 2 (the highest dose) has no
cells at all on D7/BR4 in the full table, so it has none here either. Concentration 2 is
withheld from curve fitting, so this does not affect the published fit — but it is exactly
the sort of gap a per-replicate EMD would have to handle, so keeping it in the fixture is
deliberate.

Recipe, from the repository root with the full table in `data/`:

```bash
awk 'BEGIN{FS=OFS="\t"}
NR==1{ n=0; for(i=1;i<=12;i++){keep[++n]=i}
       for(k=0;k<18;k++){keep[++n]=13+int(k*468/17+0.5)}; nkeep=n }
{ key=$1"|"$8"|"$9; if(NR>1){ if(++seen[key]>8) next }
  line=""; for(j=1;j<=nkeep;j++){ line = (j==1 ? $(keep[j]) : line OFS $(keep[j])) }
  print line }' data/cell_ID_pooled_median_row_plate_standardization_cid.txt \
  > tests/fixtures/subset.tsv
```

Columns 1-12 of the full table are the 10 metadata columns plus
`counts_RelateLysoCell` and `counts_RelateMitoCell`; 13-481 are the `rp_norm_*` features.

## `published_selected_all_days.txt` (12 KB) and `published_selected_D5.txt` (20 KB)

The two feature lists the published run selected, one per gate: across all four days
(**182** features) and on D5 alone (**374**). Taken from the headers of the published
`all_days_trimmed_features.txt` and `D5_trimmed_features.txt`, which are those lists applied
back to the data and far too large to keep (1.5 GB and 3.0 GB):

```bash
head -1 all_days_trimmed_features.txt | tr '\t' '\n' | grep -vxE \
  "Concentration|counts_Cells|counts_Cytoplasm|counts_FilteredNuclei|Metadata_Well|Metadata_Day|Metadata_Biorep|Tech_replica|Day_Well_BR|cell_ID"
```

`counts_RelateLysoCell` is the first entry of the all-days list - the `counts_` trap seen from
the far end of the pipeline. A `^counts_` metadata pattern would have classified it as
metadata, and a feature the published run retained could never have reached selection.

Both lists predate the BC4 correction, so they pair with `data/model_fit_results.txt` rather
than with a fit produced by the current pipeline. `tests/test_golden_selection.py` uses them,
and skips when that fit table is absent.
