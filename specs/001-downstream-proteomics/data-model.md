# Canonical data model v1.0.0

This is the normative logical model. JSON schemas and TSV field contracts implement it; R internal fit structures are not the public API. Strings are UTF-8, missing TSV values use `NA`, and actual empty annotations remain distinguishable from missing numeric measurements. Do not coerce identifiers to floating point. Store file paths relative to the configuration/run root when possible; hashes identify content.

## Identity and lineage

| Entity | Key | Required meaning / relationships |
|---|---|---|
| Study | study_id | Assay, taxonomy_id, experimental unit, stated design/provenance, optional tissue and treatment details; absent knowledge is null/unknown |
| Observation | observation_id | One quantified measurement column; references biological_unit_id; optional subject_id, technical_replicate_id, multiplex/channel/batch/acquisition fields |
| BiologicalUnit | biological_unit_id | Independent unit for sampling/replication; may have multiple technical observations and may belong to a repeated-measures subject |
| Feature | feature_id | Unique measured protein/group row; accession list, protein-group ambiguity, source flags and annotations; gene_symbol is never the primary key |
| MatrixArtifact | artifact_id | Numeric feature×observation matrix plus ordered IDs, scale, origin hash, missing mask, prior-imputation mask/state and preprocessing lineage |
| AnalysisPlan | plan_hash | Canonicalized validated config, inclusion/exclusion policies, design matrices, model/contrast definitions, hypothesis families, resources and code/environment identity |
| Design | design_id | Ordered observations, encoded coefficient names, finite evaluated design, term map, rank, block mode and covariate preprocessing |
| Model | model_id | Design ID, engine/version, primary/sensitivity status, coverage policy, fit settings and supported hypothesis |
| Contrast | contrast_id | Named coefficient-weight vector, estimand text, direction, biological groups/covariates involved and role; never inferred from D1/T1 prefix |
| FitArtifact | fit_id | RDS reference, model/design/plan IDs, featurewise n/df/estimability, covariance/uncertainty sufficient for downstream eligible tests |
| HypothesisFamily | family_id | Prespecified endpoint/model/contrast/universe membership, adjustment method, planned/eligible/finite counts and completeness |
| DifferentialResult | model_id+contrast_id+feature_id+hypothesis_id | Effect, SE/CI/df/statistic, P/q, eligibility reason, counts, primary/sensitivity role and family identity |
| ResourceSnapshot | resource_id+version+sha256 | Source/database/release/species/schema/terms, retrieval and mapping evidence; immutable local contents |
| GeneMapping | mapping_id+feature_id+gene_id | Source and target IDs/species, evidence, ambiguity, label-independent representative status and reason |
| GeneSetResult | model_id+contrast_id+set_id+endpoint+method | Null type, eligible universe, set overlaps, native statistic/P/q, centrally adjusted q/family and diagnostics |
| ResponseSummary | axis_id+feature_id+model_id | d/t/r with intervals/covariance, RI eligibility, descriptive category, optional valid equivalence/conjunction evidence |
| IndependentScore | score_id+training_hash | Fixed features/weights/transforms, training biological units/subjects or verifiable disjointness evidence, missing-feature rule and validation eligibility |
| Run | run_id | Frozen plan, execution directory, stage DAG, state, environment, resource checksums, start/end times and overall limitations |
| StageResult | run_id+stage_id | Inputs, outputs, status/reason, warning/error evidence, code/package hashes, RNG and duration/exit |
| ValidationEvidence | acceptance_id+commit | Test/command, independent expected result, fixture hash, actual result, runtime/package versions, reviewer and PASS/FAIL/NOT_RUN |

## Matrix and metadata contracts

Canonical abundance TSV: first column `feature_id`, then exact observation IDs. No duplicate columns or repeated feature IDs. Numeric cells are finite or NA. The feature and observation dimensions are aligned by ID, not incidental row order. Canonical long input contains feature_id, observation_id, abundance and can be losslessly pivoted only when each pair is unique.

Observation metadata TSV fields: observation_id, biological_unit_id, subject_id (nullable), technical_replicate_id (nullable), group (nullable for continuous-only designs), batch (nullable), assay_batch, acquisition_run, multiplex_id, channel, is_reference_channel and arbitrary declared covariate columns. Semantic validation checks each design's required fields and does not require irrelevant TMT fields for LFQ. Biological-unit/subject role is explicitly declared rather than inferred from string prefixes. Technical duplicates require the configured aggregation rule before inferential count calculation.

Feature metadata TSV: feature_id, accessions (JSON array or escaped delimiter with schema), gene_ids, gene_symbols, description, is_decoy, is_contaminant, protein_group_ambiguous, source_database and annotation_version. Flags allow true/false/unknown. Count evidence is separate: feature_id, optional observation_id/plex_id, count_type, count_value, count_source, count_aggregation and pseudocount_policy. Missing counts are not zero.

Masks have identical shape/ordered IDs to the matrix. `observed_mask` records genuinely observed measurements if known; `input_numeric_available_mask` records finite input cells; `imputed_mask` records declared replacements. If original observation information is unavailable for imputed data, observed_mask is unavailable, never all true by assumption. Every transform emits new artifact IDs and a transform record with input scale/output scale/parameters/factors/mask changes.

## Differential export fields

Required identifiers: run_id, plan_hash, model_id, design_id, contrast_id, feature_id, hypothesis_id, engine, engine_version, role, family_id. Required states: eligibility (`tested`, `excluded`, `nonestimable`, `numerical_failure`), reason_code, input_scale, effect_scale and n_biological/observed counts per relevant group (structured companion table allowed).

Numeric fields: effect, effect_se, ci_lower, ci_upper, ci_level, ci_method, statistic, statistic_type, df_residual, df_inference, p_value, q_value, native_q_value, effect_threshold (nullable), null_value/null_region and mean_abundance. Backend-unavailable quantities are NA with field-level reason; unsupported hypotheses fail eligibility. proDA uncertainty may have distinct variance/df meaning and must not be mislabeled limma. TREAT's threshold P is distinct from the zero-null moderated P. Preserve any native engine fields in a namespaced sidecar to avoid losing provenance.

Family summary: family_id, definition_hash, adjustment, dependence_assumption, n_planned, n_eligible, n_finite, n_numerical_failure, completeness, q_cutoff and rejection_count. Numeric failures in eligible planned tests make completeness false and constrain claims. A missing test is never silently removed from the scientific report.

## Pathway exports

Fields: run/plan/model/contrast IDs, method, null_type, endpoint (`competitive_directional`, `rotation_directional`, `rotation_mixed`, `preranked`, `ora_up`, `ora_down`), resource_snapshot/hash, source/target_taxonomy, mapping_hash, universe_hash/universe_size, set_id/name, original_size, mapped_size, eligible_overlap_size, statistic/statistic_type, direction, p_value, native_q_value, q_value, family_id, status/reason, and diagnostic sidecar.

fgsea diagnostics include ES, NES, log2err, leading-edge membership, finite rank/tie counts and seed. CAMERA includes estimated/fixed inter-gene correlation and gene-model details. mroast preserves NGenes, PropDown/PropUp, Direction, directional P and mixed P in unambiguous fields; raw endpoints are adjusted under declared families rather than confusing native FDR with central q. ORA retains N/K/n/k and foreground/background/overlap gene IDs. Every set result links to membership lists used, not only a label.

## Response and score exports

Response includes axis definition, disease/treatment/residual contrast IDs, their effect vectors and covariance hash, d/t/r intervals, denominator minimum, epsilon, RI/RI_method/RI_status/interval status, descriptive_class and overshoot flags. Inference fields identify hypothesis/independent_direction_source/equivalence_margin, both one-sided P values, max-P, interval level, family/q and applicability reason. Do not represent `not tested`, `not significant`, and `equivalent` with the same boolean.

A fixed score lists each feature_id, weight, training_center, training_scale, training_missing_rule and training source identifier/hash. Metadata states biological-unit/subject independence evidence. Score output includes validation-unit IDs, missing-feature fraction, applied fixed transforms, score and eligibility. The test artifact records null/randomization unit/block scheme, all versus sampled allocations, tail definition, tie tolerance, k, N/B, seed, P and Monte Carlo uncertainty. No inferential score artifact is emitted for current-data-selected weights unless a separate approved inference contract exists.

## Versioning and serialization

Schemas reject unknown top-level scientific fields to catch typos. Metadata covariate maps are explicitly extensible with declared types; this does not authorize unknown model methods. Additive fields require schema/minor-version handling; changed meaning requires major version and migration tests. Use deterministic column order and stable sorting for semantic comparison, with proper CSV/TSV escaping and roundtrip tests for Unicode, commas, semicolons and NA strings in annotations. Float serialization preserves meaningful precision; reports format copies, never round underlying tables before inference.

## Artifact state and trust

Only a completed stage with validated output schema/hash can feed downstream analysis. A directory existing or a CSV having rows does not prove completion. Statistics do not read their terminal validation ledger to select a successful result; validation observes computation and may block release, not alter the method to make a gate green. Private source locations are resolved at runtime and omitted/redacted from public release metadata while their content hashes and abstract provenance remain available.
