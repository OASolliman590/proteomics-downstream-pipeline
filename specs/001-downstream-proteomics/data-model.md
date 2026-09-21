# Canonical data model v1.2.0

Status: frozen preimplementation contract; all maintained artifacts below are future outputs. [Scientific meaning](contracts/scientific-methods.md), [configuration](contracts/analysis.schema.json) and [CLI/state protocol](contracts/cli-and-artifacts.md) are authoritative in their respective domains.

## Identity and serialization

UTF-8 text; TSV delimiter is TAB and newline LF. Numeric missing values serialize as NA; JSON uses null, never NaN/Infinity. Identifiers remain strings. JSON arrays embedded in TSV encode accessions and membership without splitting ambiguous delimiters. Annotation missingness is separate from numeric missingness. Numeric tables retain full double precision (17 significant digits); formatting belongs to a report copy. Stable sort keys and column order are part of each artifact schema. Unknown flags use unknown, not false.

| Entity | Key | Fields and meaning |
|---|---|---|
| Study | study.id | Explicit assay, organism taxonomy, experimental-unit/provenance statement, optional tissue. |
| Observation | observation_id | biological_unit_id, subject_id (nullable), technical_replicate_id (nullable), group and declared covariates; assay-specific plex/channel/reference fields when needed. |
| Feature | feature_id | Stable protein/group ID, constituent accessions, gene IDs/symbols, ambiguity and tri-state source flags. Gene symbol is not a row key. |
| MatrixArtifact | artifact_id | Ordered feature/observation IDs, matrix path/hash, scale state, numeric-availability mask, original-observed mask or unavailability reason, prior-imputed mask/state and transform lineage. |
| AnalysisPlan | plan_hash | Resolved config, retained IDs/masks, evaluated designs/contrasts, requested capabilities with explicit required/optional intent and separately inspected availability, model/hypothesis/family definitions, input/code/environment/resource hashes and runtime seeds. |
| Design | design_id | Observation order, reversible coefficient map, design/contrast matrices, rank/alias/estimability records, blocking mode and fixed centering values. |
| ModelFit | model_id + plan_hash | Engine/version, design, fitting universe, weights/block correlation, exactness path, featurewise n/df and native fit artifact. |
| HypothesisFamily | family_id + plan_hash | Prespecified endpoint membership, adjustment, universe hash, planned/eligible/finite/failure counts, completeness and rejection count. |
| DifferentialResult | model_id + contrast_id + feature_id + hypothesis_type | Typed endpoint, full-precision estimate/uncertainty, family and estimability. See fields below. |
| ResourceSnapshot | resource_id + version + sha256 | Local source bytes, source/target taxonomy, ID type, source/release/build/terms and orthology evidence. |
| GeneMapping / GeneMatrix | feature_id + gene_id + mapping_hash | Candidate/representative/excluded state, deterministic tie keys, memberships, finite-matrix universe and model lineage. |
| GeneSetResult | model_id + contrast_id + set_id + hypothesis_type + engine | Actual null/statistic, mapped/eligible set membership, universe, native and central adjusted values and diagnostics. |
| DescriptiveResponse | axis_id + model_id + feature_id | d/t/r with covariance reference, RI/class/flags and interval availability; no inferential P/q fields. |
| IndependentScore | score_id + resource_hash | Frozen feature weights/signs, training centers/scales, missing-feature policy, selection/training participant provenance and independent-validation eligibility. |
| StageResult / RunStatus | run_id + stage_id / run_id | Executed state and typed reason, command/exit/timestamps, requested phase/capabilities, inputs/outputs and warnings. |
| ValidationEvidence | acceptance_id + reviewed_tree | Expected fixture/oracle, command/versions, actual observations, artifacts/hashes and review conclusion. Missing artifacts remain an empty array. |

## Canonical tables and masks

Wide matrix: feature_id followed by observation IDs. Long matrix: feature_id, observation_id, abundance; pairs are unique. Feature and observation metadata must align by key, not row order. Minimum observation columns are observation_id, biological_unit_id, subject_id, technical_replicate_id, group, plus every configured covariate. Technical IDs may be missing only when no technical replicate exists. Additional columns are retained as annotations but not silently promoted to covariates. Paired visits are distinct biological specimens but share subject_id; technical injections of a specimen share biological_unit_id and condition.

Minimum feature columns are feature_id, accessions, gene_ids, gene_symbols, is_decoy, is_contaminant, protein_group_ambiguous. Unknown accessions/genes are empty JSON arrays, not invented annotations. Counts are a separate evidence table with feature_id, optional observation_id/plex_id, count_type (peptide or psm), count_value, count_source, count_aggregation and pseudocount_policy. Configured model count aggregation/zero policy must agree with these recorded provenance fields. The table never treats missing count as zero.

All masks have identical matrix keys/order. A genuine original mask is derived only for documented nonimputed input or supplied/validated otherwise. Aggregation retains the source-to-specimen map and coverage rule. Every transform has source_artifact_id, output_artifact_id, input_scale, output_scale, parameters, source/output hashes and changed-mask diagnostics.

## Shared scientific result-row envelope

Every scientific endpoint row (differential, detection, pathway, response, score) includes:

`schema_version`, `run_id`, `plan_hash`, `result_type`, `hypothesis_type`, `family_id`, `engine`, `engine_version`, `model_id`, `contrast_id` (nullable for an axis), `estimable`, `n_obs_by_required_group`, `eligibility`, `reason_code`.

`n_obs_by_required_group` is a JSON object embedded in TSV; counts refer to genuine observed biological units. Unknown counts are null plus a reason, never numeric availability mislabeled observed. `family_id=null` for descriptive results, which cannot carry p_value/q_value/native_q_value. Input/QC/mapping diagnostics are not inferential result rows and have their own keys; they do not acquire a fake family or hypothesis.

`eligibility` is tested, excluded, nonestimable, inapplicable, not_run or numerical_failure. `estimable` states whether the declared endpoint is mathematically identifiable, not whether it was run. Nonexecuted mathematical eligibility not evaluated is null with a reason; tested results have true; rank/coverage failures have false. JSON boolean/null is rendered true/false/NA in TSV.

## Differential outputs

Identifiers additionally include design_id and feature_id. Required numeric/semantic fields: effect, effect_scale, effect_se, ci_lower, ci_upper, ci_level, ci_method, statistic, statistic_type, df_residual, df_inference, p_value, q_value, native_q_value, effect_threshold, null_region, mean_abundance, role. Unavailable numeric values are NA with field-reason entries in the diagnostics sidecar. A zero-null row uses result_type=ProteinZeroNullResult and hypothesis_type=protein_zero_null; TREAT uses ProteinTreatResult/protein_treat; omnibus uses ProteinOmnibusResult/protein_omnibus. Model-native fields are preserved in a namespaced sidecar. An ordinary 95% effect interval does not become a TREAT significance interval.

Every planned feature/contrast endpoint remains in the table, including excluded/nonestimable rows with reasons. A model failure writes a stage failure, not fabricated “tested” rows. Family summaries include family_id, definition_hash, adjustment, dependence_assumption, n_planned, n_eligible, n_finite, n_numerical_failure, completeness, q_cutoff, rejection_count and universe_hash. Do not pool sensitivity/primary or different nulls.

## Pathway outputs

Add set_id, set_name, resource_id/hash, source_taxonomy_id, target_taxonomy_id, mapping_hash, universe_hash/size, original_size, mapped_size, eligible_overlap_size, direction, statistic/statistic_type, p_value, native_q_value, q_value and diagnostic reference. CAMERA uses competitive_enrichment; ROAST separate self_contained_directional and self_contained_mixed rows; fgsea preranked_gene_set; ORA ora_up/ora_down. Preserve ES/NES/log2err/leading-edge for fgsea, correlation for CAMERA, native directional/mixed fields for ROAST, and N/K/n/k plus actual gene IDs for ORA.

## Response, equivalence and score outputs

`DescriptiveResponse`/`descriptive_reversal` stores axis_id, disease/treatment/residual contrast IDs, d/t/r, covariance artifact/hash, dmin, RI, descriptive_class, crossed_control, residual_exceeds_disease and ratio interval kind. It has no P/q fields. A requested `RatioInterval` may be bounded, unbounded, disjoint or unavailable; disjoint limits serialize as an array of intervals, not a misleading single finite CI.

`EquivalenceResult`/`equivalence` is a separate table with margin, alpha, residual estimate/SE/df, p_lower, p_upper, p_value=max(one-sided P), q_value and equivalence CI/level. `FormalRescueResult`/`formal_rescue` is another table with independent-direction resource/hash, named component P values and their maximum, family and eligibility. Neither is a boolean pasted onto a descriptive class.

`DescriptiveScore`/`descriptive_score` has validation-unit ID, value, observed-feature fraction, transform/selection hashes and score_eligibility_reason, but physically no inferential P/q/statistic columns. `IndependentScoreTest`/`independent_score_test` is separate: score_id, admissible scheme/unit/block, exact N or Monte Carlo B, extreme count k, observed/forced allocation count, tie rule, RNG, p_value/q_value and Monte Carlo interval. Selection/train/test overlap or unknown provenance prohibits that inferential object.

## Pending evidence

`traceability.json` separates implementation_targets and proposed_test_target from actual evidence. Prior to execution: status=pending-after-freeze, evidence=[], verified_commit=null. Expected file names in a packet are creation targets, not proof that the maintained software or results exist. State transitions and artifact promotion are defined only in [cli-and-artifacts.md](contracts/cli-and-artifacts.md).
