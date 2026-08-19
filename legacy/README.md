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
verified in the package. **Every one of them now is** - the map below has no open rows - so the
directory is ready to go; it is kept only until the figures have been signed off by eye. Where a new module
is substantially a descendant of one script, the migration uses `git mv` so that
`git log --follow` on the module reaches the original author's history.

## Migration map

| Inherited script | Destination | Fate |
|---|---|---|
| `plots_emd_model_drc.py` (model functions + AIC/BIC) | `models.py` | done — lifted almost verbatim; core maths, do not touch the formulas |
| `plots_emd_model_drc.py` (everything else) | `fit.py` + `plots.py` | done — fitting produces a table, plotting reads it |
| `emd_scores_concs_per_day.py` | `emd.py` | done — rewritten as a contrast-driven engine |
| `emd_scores_controls_trimming_well_results.py` | `emd.py` | done — merged; it is the `baseline` contrast set, not a separate stage |
| `select_features.py` | `select.py` | done — rewritten, quantifiers made explicit and config-driven |
| `correlation_feature_selection_well_batch.py` | `prune.py` + `plots.py` | done — maths kept (`1-|r|`, average linkage at 0.1), and the tree it cuts is now the tree it draws |
| `parsing_clusters.py` | `prune.py` | done — clusters are a TSV column now, so nothing `eval()`s a line read from a file |
| `trimming_value_include_batch_v1_cid.py` | `subset.py` | done — folded into one step |
| `trimming_value_include_batch_v2_cid.py` | `subset.py` | done — folded into one step, it was near-identical to v1 |
| `plot_emd_controls.py` | `plots.py` | done — folded in, and it now runs on either distance table |

The trimming loop these scripts share is also the ancestor of `trim.py`, together with
`scripts/standardization/only_trimming_well.py` in HCS-proc, which produced the trimmed
table the published EMD step consumed. `trim.py` reproduces that step so that `cdr_FS`
can start from the untrimmed standardized table.

`examples/published.ini` is the inherited `config.ini`, moved and rewritten.
