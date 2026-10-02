# PERMANOVA configuration and artifact contract (R13)

Normative statistics are SM26–SM30 in the [scientific contract](../../001-downstream-proteomics/contracts/scientific-methods.md#multivariate-permanova). This file fixes the configuration block, outputs and vocabulary.

## Configuration (`multivariate`, optional top-level block of analysis.schema.json)

| Field | Meaning | Default |
|---|---|---|
| `enabled` | Request the `permanova` capability. | — (required) |
| `execution_requirement` | `required` or `optional`; independent of installation state. | `optional` |
| `permanova.distance.metric` | `euclidean` or `manhattan` (vegan::vegdist names). | `euclidean` |
| `permanova.distance.scaling` | `zscore` (per feature, within each test's observations), `center` or `none`. | `zscore` |
| `permanova.permutations` | Permutation count; production profile requires ≥ 999. | 9999 |
| `permanova.seed` | RNG seed (L'Ecuyer-CMRG); every adonis2/permutest call re-seeds with it. | `runtime.seed` |
| `permanova.pairwise` | Run all pairwise post-hoc tests. | `true` |
| `permanova.adjustment` | `holm`, `bonferroni`, `BH` or `none`, within each feature set's pairwise family. | `holm` |
| `permanova.alpha` | Threshold used only for the interpretation label. | 0.05 |
| `permanova.covariates` | `[{column, type: categorical|continuous}]` observation metadata columns. | `[]` |
| `permanova.interaction` | Fit group × each categorical covariate (only when every cell is non-empty). | `false` |
| `permanova.feature_sets` | `[{id, kind: all_complete}]`, `[{id, kind: declared_panel, feature_ids, selection_provenance: independent|same_data|unknown}]`, `[{id, kind: dep_derived, model_id, contrast_id, criterion: family_q|raw_p, threshold}]`. | `[{id: "all", kind: "all_complete"}]` |
| `permanova.random_sets` | Draws for the random equal-size feature-set null. | 2000 |

## Feature universe

Genuinely observed cells of the frozen primary matrix (prior-imputed cells are treated as missing; nothing is imputed). A feature enters a test only when it is observed in every analysed observation of that test and is not constant within it. Excluded features are written with reasons (`incomplete_in_test_observations`, `constant_within_test`, `not_in_primary_matrix`). A set with fewer than two usable features is refused (E_PERMANOVA_FEATURE_SET).

## Permutation schemes

`free` for the group term of independent designs; `within:<subject_column>` for the group term when the design blocks subjects (paired/repeated; covariates are then unsupported, E_PERMANOVA_DESIGN_UNSUPPORTED); `within:<group_column>` for covariate and interaction terms. PERMDISP uses the scheme of the term it accompanies; continuous covariates have no PERMDISP (`not_applicable`).

## Outputs (`<run>/permanova/`)

| File | Content |
|---|---|
| `tests.tsv` | One row per feature set × analysis (`global`, `pairwise`, `group_adjusted`, `covariate`, `interaction`) × term: n, n_features, df, ss, r2, pseudo_f, residual_df, residual_ss, total_ss, p_value, p_floor, p_display, p_adjusted, p_adjusted_display, adjustment, nperm, seed, permutation_scheme, permdisp_f, permdisp_p, permdisp_p_display, permdisp_scheme, interpretation, metric, scaling. |
| `feature_usage.tsv` | Every candidate feature per set and test with `used` and `reason`. |
| `feature_r2.tsv` | Per-feature univariate group R² and group means per set (global groups). No per-feature P values (protein-wise inference belongs to R05). |
| `identity_check.json` | Mean per-feature R² versus multivariate R² per set; `pass`, `fail` or `not_applicable`. |
| `selection_context.tsv`, `random_set_null.tsv` | Random equal-size null (draws, seed), best-possible same-size R², fraction of null ≥ the set's R², optimism flag. |
| `permutation_null.tsv` | Permuted pseudo-F values for each reported term (source data of the null histograms). |
| `ordination.tsv`, `ordination_axes.tsv`, `centroids.tsv` | PCoA coordinates, centroids and distances to centroid (source data of ordination and PERMDISP figures). |
| `figures/*.png` (+ requested pdf/svg) | Ordination with centroids and spider lines, PERMDISP boxplots, permutation-null histograms, R² bars, per-feature R², random-set null. |
| `permanova_result.json` | Settings, versions, seeds, states, interpretation rule and limitations. |

## Interpretation vocabulary

With α from the configuration and the row's P (adjusted P for pairwise rows):

| PERMANOVA | PERMDISP | `interpretation` |
|---|---|---|
| P < α | P ≥ α | `location_shift` |
| P < α | P < α | `location_and_or_dispersion` |
| P ≥ α | P < α | `dispersion_difference_location_not_established` |
| P ≥ α | P ≥ α | `no_separation_detected` |
| any | not applicable | `permdisp_not_applicable` |

The report states the rule: a significant PERMANOVA with a non-significant PERMDISP indicates a location (centroid) shift; when both are significant, location and/or dispersion differences are possible. Non-significance is not evidence of no difference.
