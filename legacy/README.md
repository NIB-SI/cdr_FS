# Inherited scripts — staging area, not part of the package

These are the nine scripts that came through the `git subtree split` of
`NIB-SI/HCS-proc`'s `scripts/feature_selection/`. They are **the record of what was
actually run** for the article, and they are kept here, unmodified, as the reference
against which the rewritten package is checked.

They are deliberately **outside `src/cdr_fs/`**: each one reads `config.ini` from the
current working directory and loads its input table at import time, so none of them can
live inside an importable package. They are also not installed — `pyproject.toml`
packages `src/` only.

This directory is temporary. Each script is deleted once its behaviour is reproduced and
verified in the package; the directory goes away at the end of Phase 5. Where a new module
is substantially a descendant of one script, the migration uses `git mv` so that
`git log --follow` on the module reaches the original author's history.

## Migration map

| Inherited script | Destination | Fate |
|---|---|---|
| `plots_emd_model_drc.py` (model functions + AIC/BIC) | `models.py` | lifted almost verbatim — core maths, do not touch the formulas |
| `plots_emd_model_drc.py` (everything else) | `fit.py` + `plots.py` | rewritten: fitting produces a table, plotting reads it |
| `emd_scores_concs_per_day.py` | `emd.py` | rewritten contrast-driven engine |
| `emd_scores_controls_trimming_well_results.py` | `emd.py` | merged — it is the `baseline` contrast set, not a separate stage |
| `select_features.py` | `select.py` | rewritten, quantifiers made explicit and config-driven |
| `correlation_feature_selection_well_batch.py` | `prune.py` | rewritten as an optional module; maths kept (`1-|r|`, average linkage at 0.1) |
| `parsing_clusters.py` | `prune.py` | **must not survive as-is** — it `eval()`s a line read from a file |
| `trimming_value_include_batch_v1_cid.py` | `subset.py` | folded into one step |
| `trimming_value_include_batch_v2_cid.py` | `subset.py` | folded into one step — near-identical to v1 |
| `plot_emd_controls.py` | `plots.py` | folded in |

The trimming loop these scripts share is also the ancestor of `trim.py`, together with
`scripts/standardization/only_trimming_well.py` in HCS-proc, which produced the trimmed
table the published EMD step consumed. `trim.py` reproduces that step so that `cdr_FS`
can start from the untrimmed standardized table.

`examples/published.ini` is the inherited `config.ini`, moved and rewritten.
