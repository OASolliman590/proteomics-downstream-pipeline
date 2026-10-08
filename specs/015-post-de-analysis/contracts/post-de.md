# Post-differential analysis configuration and artifact contract (R14a–R14f)

Normative statistics are SM31–SM41 in the scientific contract (v1.3.0). This file fixes the configuration block, outputs and vocabulary. All identifiers below are declared by the user and validated against the frozen plan. Nothing is study-specific.

## Configuration: `post_de` (optional top-level block of analysis.schema.json)

Every sub-block has `enabled` and `execution_requirement` (`required` | `optional`, default `optional`). Phase 1–3 configurations cannot enable `post_de`.

### `post_de.sets` (R14a)

| Field | Meaning | Default |
|---|---|---|
| `definitions` | `[{id, rule}]`. `rule` is a JSON tree: `{op: union|intersect|difference|complement_within_tested, args: [...]}`. Leaves are `{family_id, contrast_id, criterion: family_q|treat|exploratory_raw_p, threshold, direction: up|down|any}`. | — |
| `venn` | Set IDs to draw as circle Venn (k ≤ 3). | `[]` |
| `upset` | Set IDs for UpSet (any k). | all defined |
| `concordance_pairs` | `[{left: {family_id, contrast_id}, right: {...}}]` | `[]` |
| `overlap_test` | Request a hypergeometric overlap P (only reported when the contrasts share no biological units). | `false` |

### `post_de.sensitivity` (R14b)

| Field | Meaning | Default |
|---|---|---|
| `covariate_models` | `[{id, add_covariates: [column], interaction_with_group: bool}]` | `[]` |
| `subgroups` | `[{id, filter: {column, in: [levels]}}]`, applied to observation metadata. | `[]` |
| `matched_n` | `{targets: {group: n}, draws: 500, seed}` | off |
| `influence` | Leave-one-biological-unit-out on declared families. | `false` |
| `criteria` | Membership criteria used in comparisons (as in `sets` leaves). | primary family q < 0.05 |

### `post_de.association` (R14c)

| Field | Meaning | Default |
|---|---|---|
| `phenotypes` | `[{column, type: numeric|ordinal, adjust_for: [column], adjust_for_group: true, scope: pooled|within_group}]` | `[]` |
| `method` | `model` (primary engine) or `correlation` (`pearson`/`spearman`/`partial`) | `model` |
| `permutations` | For the correlation method | 9999 |

### `post_de.biomarker` (R14d)

| Field | Meaning | Default |
|---|---|---|
| `contrast` | Binary declared contrast (two groups). | — (required) |
| `min_units_per_class` | Refusal threshold. | 5 |
| `single_feature` | `{features: all_tested|set_id, direction: prespecified|train_only, ci: [delong, bootstrap], bootstrap: 2000}` | as shown |
| `cv` | `{outer: {scheme: repeated_stratified_kfold|loocv, k: 5, repeats: 100}, inner: {k: 5}, group_by_subject: auto}` | as shown |
| `selection` | `{method: top_k_auc|top_k_t|elastic_net|lasso, k: [3,5,10] or penalty grid}` | `top_k_auc`, k ∈ {5, 10} |
| `classifier` | `penalized_logistic` | `svm_linear` | `svm_polynomial` | `svm_radial` | `random_forest` | `penalized_logistic` |
| `threshold_rule` | `youden_train` | `fixed:<value>` | `youden_train` |
| `permutation` | `{B: 1000, seed, compute_guard_hours}` | B = 1000 |
| `compute_policy` | `adapt` (scale repeats, then B, down to `compute_guard_hours`; minimums 10 repeats and B = 99; recorded) or `refuse` (strict) — amended 2026-10-05, D-59 | `adapt` |
| `fixed_panels` | `[{id, feature_ids, provenance: independent|same_data|unknown}]` | `[]` |
| `validation_cohort` | `{matrix_artifact, metadata, subject_column}` (must pass SM23 audit) | none |

### `post_de.networks` (R14e)

| Field | Meaning | Default |
|---|---|---|
| `coabundance` | `{enabled, min_units: 20, rule: wgcna_signed|hclust_correlation, soft_threshold: auto|int, min_module_size: 10, bootstrap: 200, traits: [group|phenotype columns]}` | off |

Amended 2026-10-05 (ADR 0010, D-59): `min_units` is the recommended level; between 4 units and it co-abundance runs as exploratory, and below 4 it is refused. Repeated designs are analysed on subject means (group constant within subject) or after subject-effect removal (group varying within subject). `ppi.null_universe` other than `measured` is refused and the measured-universe null runs. For `association`, `pearson`/`spearman` with an adjustment run in their partial form; `missing` other than `complete_case` is refused and complete case runs. For `biomarker`, a leaky declaration (`threshold_rule`, `permutation.scope`, `imputation`, `validation_cohort.retune`, `cv.group_by_subject: false` in a blocked design) is refused and its leakage-safe form runs. `train_only` (single-feature direction) chooses the direction on the evaluated data and is labelled `in_sample`.
| `ppi` | `{snapshot_id, min_score, sets: [set_id], null_draws: 2000, seed}` | off |

## Claim vocabulary (frozen)

| Label | Meaning |
|---|---|
| `descriptive` | No inferential claim; summary or visual only. |
| `exploratory_raw_p` | Membership by unadjusted P; not a discovery set. |
| `in_sample` | Estimated on the same samples used to choose the direction, features or model. |
| `cross_validated_nested` | Selection and tuning inside training folds; performance on held-out folds. |
| `fixed_panel_cv` | Classifier cross-validated, panel fixed in advance; optimistic if provenance is `same_data`/`unknown`. |
| `independently_validated` | Locked model evaluated once on a disjoint declared cohort. |
| `module_level` | Inference on module eigengenes, not individual proteins. |

## Outputs (`<run>/post_de/<module>/`)

| Module | Files |
|---|---|
| sets | `membership.tsv`, `regions.tsv` (region id, defining sets, protein-group and unique-gene counts, members), `direction_overlap.tsv`, `concordance.tsv`, `overlap_test.tsv` (eligible pairs only), figures (`upset_*.png`, `venn_*.png`, `concordance_*.png`) + source tables |
| sensitivity | `imbalance.tsv`, per-model DEA tables and families, `comparison.tsv`, `matched_n.tsv`, `influence.tsv`, `robustness_summary.tsv`, figures + sources |
| association | `association_<phenotype>.tsv`, `heatmap_source.tsv`, figures |
| biomarker | `single_feature_auc.tsv`, `cv_performance.tsv`, `oof_predictions.tsv`, `selection_stability.tsv`, `permutation_null.tsv`, `calibration.tsv`, `fixed_panels.tsv`, `validation.tsv`, `transform_audit.jsonl`, figures + sources |
| networks | `modules.tsv`, `eigengenes.tsv`, `module_stability.tsv`, `module_trait.tsv`, `ppi_subnetwork.tsv`, `connectivity.tsv`, figures + sources |
| all | `eligibility.json` (per module: state, typed reason, claim label, input hashes, adaptations: analysis, item, requested, used, reason), `stage-result.json` |

Every figure has a source table from which every plotted coordinate can be re-derived (SM25).
