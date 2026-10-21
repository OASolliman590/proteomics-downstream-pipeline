# Methods extension configuration and artifact contract (R15a–R15d)

Normative statistics are SM42–SM46 in the scientific contract (v1.4.0). This file fixes the configuration block, the study manifest, the error codes, the calibration gates and the outputs. All identifiers are declared by the user and validated against the frozen plan. Nothing here is study-specific. Unknown keys are refused.

## Configuration: `methods_extension` (optional top-level block of analysis.schema.json)

Every sub-block has `enabled` (default `false`) and `execution_requirement` (`required` | `optional`, default `optional`). A configuration without `methods_extension` produces the same outputs as before this slice. Requests are declarative. There is no expression evaluation.

### `methods_extension.welch` (R15a)

| Field | Meaning | Default |
|---|---|---|
| `contrasts` | Contract IDs from the frozen plan. Each must be a pairwise group contrast on log2 data. | all declared pairwise contrasts |
| `role` | `primary` (single primary for the hypothesis, ADR 0002) or `secondary`. | `secondary` |
| `family` | Family ID. A secondary Welch family is always named and separate. | `welch_<hypothesis>` |
| `paired` | `auto` (from the frozen design), `true` or `false`. A declared paired design with `false` is refused as `E_WELCH_PAIRED_UNIT`. | `auto` |
| `ci_level` | Confidence level for the effect interval. | plan CI level |
| `paper_compatible` | `{enabled, raw_p_cutoff}`. Selects features with raw P below the cutoff. Always labelled `not_multiplicity_controlled` (SM43). | `{enabled: false, raw_p_cutoff: 0.05}` |

**Role mapping (F4).** The `role` field above uses the spec's names. The schema's model role admits only `primary` and `sensitivity`, so a secondary Welch model is declared as `models[].role: sensitivity` and `methods_extension.welch.role: secondary`. The mapping is checked in `config.py`. It does not trigger `E_SENSITIVITY_PRIMARY_MODEL`, which concerns only a preprocessing sensitivity matrix assigned to the primary model. A sensitivity matrix assigned to a welch model is refused with `E_WELCH_MATRIX_UNSUPPORTED`, because the welch engine reads the observed primary matrix only. A welch primary writes the `dea/` contract (`zero_null.tsv`, `families.tsv`, `model_status.json`), so it cannot share a run with limma models (`E_PRIMARY_MODEL_CONFLICT`).

The methods extension is Phase 5 only (`runtime.phase` 5). Phase 1 to 4 configurations cannot enable `welch` or `pathways`.

### `methods_extension.pathways` (R15b)

| Field | Meaning | Default |
|---|---|---|
| `libraries` | `[{id, path, sha256, source, version, species, identifier_type, licence_id, redistributable}]`. `path` points to a local file supplied by the user. | `[]` |
| `methods` | `ora` (hypergeometric, measured universe, SM17) and/or `prerank` (exploratory fgsea adapter, FR-075). | `["ora"]` |
| `rank_statistic` | Fixed to the statistic of the primary zero-null table (`dea/zero_null.tsv`): the Welch t of FR-168 for a welch primary. Not configurable in this release. | Welch t |
| `ora_threshold` | Family q at or below which a feature is a foreground hit (direction from the effect sign). | `0.05` |
| `execution_requirement` | `required` or `optional` (a required library stage that fails stops the run). | `optional` |

The stage runs after the primary differential stage in a run. Its inputs are the primary table with one representative gene per feature and contrast (the smallest stable feature ID wins a tie), the measured universe is the mapped genes of the run, and every library is checked by SHA-256 at plan time (`E_RESOURCE_HASH`) and by the licence rule (`E_LIBRARY_REDISTRIBUTION`). The run catalogue draws `path_enrichment_dotplot` for each library and direction with at least one foreground hit.
| `min_size`, `max_size` | Set-size eligibility, applied to set ∩ universe before any test (SM17). | SM17 defaults |

A library whose `licence_id` is KEGG or MSigDB must have `redistributable: false`. Such a library is usable locally and is never copied into an artifact (E_LIBRARY_REDISTRIBUTION).

### `methods_extension.classifiers` (R15c, extends R14d)

R14d `classifier` values gain `lasso` (penalty alpha = 1) and `elastic_net` (`alpha` declared, default 0.5, penalty `lambda` tuned by the inner loop). Both use the R14d leakage-safe procedure, the same outer and inner schemes, and the same permutation rule. Missing glmnet yields NOT_RUN.

The outer scheme also accepts `loocv` (R14d `cv.outer.scheme: loocv`, with `cv.inner.k: 5`) under the same rules. *Amended 2026-10-21 (D-77):* its AUC is the leave-pair-out AUC (`cv_performance.tsv` column `auc_estimator` = `leave_pair_out`); the leave-one-out scores give the threshold metrics. Mixed-class units are adapted to repeated stratified k-fold at plan time (recorded) and refused by the procedure with `E_BIOMARKER_LPO_UNITS`.

### `methods_extension.studies` (R15d)

| Field | Meaning | Default |
|---|---|---|
| `manifest` | Path to the study manifest (below). | none |
| `slots` | Slot IDs to run: `dda_label_free`, `paired_repeated`, `tmt_batched`, `dia`, `multigroup_phenotype_nonhuman`, `reference_welch`. | all declared |
| `public_approval` | `{approved_by, approved_on, accessions, shortlist_sha256}`. Required before any public download (FR-187). | none |

## Study manifest (R15d, FR-185)

A JSON document, `schema_version` `1.0.0`.

| Field | Required | Meaning |
|---|---|---|
| `study_id` | yes | Local identifier. Must not contain a file or person name. |
| `origin` | yes | `private` or `public`. |
| `slot` | yes | One of the six slot IDs above. |
| `design` | yes | `two_group`, `paired`, `repeated`, `multi_group`, with `batch` or `plex` where applicable. |
| `input` | yes | `{format, files: [{name, sha256, bytes}]}`. `format` is one of `maxquant_proteingroups`, `diann`, `spectronaut`, `canonical_matrix`. |
| `organism` | yes | Declared organism resource key. |
| `columns` | yes | `{subject, group, batch, phenotype}` as applicable. |
| `licence` | yes | Licence identifier of the data. |
| `approval_ref` | for `public` | Reference to the `public_approval` record. |

Rules: every input file hash is checked before `validate`. A `private` manifest is never committed; only its status record is (FR-188). Unknown design types are refused (E_MANIFEST_DESIGN).

## Calibration gates

These gates (R15c, FR-183) are fixed before the first run. The gate file and its SHA-256 are written to `docs/validation/016-methods-extension-validation/R15c/gates-preregistration.json` before any run, and the hash is recorded in the R15c receipt. A change after the run fails E_CALIBRATION_GATE_CHANGED.

| Gate | Scenario | Rule | Source |
|---|---|---|---|
| `nested_auc_null` | null, each repeated-stratified family | \|mean − 0.5\| ≤ 0.03 and mean + 1.96·SE ≤ 0.53 | R14f, carried over unchanged |
| `permutation_p_null` | null, each family | binomial P(K ≥ k \| rate 0.05) ≥ 0.01 and chi-square goodness-of-fit P over the B + 1 values ≥ 0.01 | R14f, carried over unchanged |
| `leaky_reference` | leaky selection and scaling on all data | mean AUC > 0.75 | R14f, carried over unchanged |
| `loocv_null` | null, LOOCV outer | same rule as `nested_auc_null` on the pooled out-of-fold AUC | R14f rule, applied to LOOCV (FAIL 2026-10-13; superseded by `lpo_null`, D-77) |
| `lpo_null` | null, LOOCV outer with the leave-pair-out AUC (ridge, lasso) | same rule as `nested_auc_null` on the leave-pair-out AUC | pre-registered 2026-10-21 in `R15c/followup/` (D-77) |
| `rf_permutation_diagnostic` | null, random forest, B = 99, 200 datasets | the `permutation_p_null` rules, plus a 10-bin chi-square; tied permutation AUCs and the mid-P recorded | pre-registered 2026-10-21 in `R15c/followup/` (D-78); diagnostic of the V180 FAIL |
| `blocked_null` | null, subject-blocked | the R14f blocked nested-AUC and blocked permutation rules | R14f, carried over unchanged |
| `blocked_ungrouped_reference` | subject-blocked, ungrouped folds | mean ungrouped AUC − mean grouped AUC > 0.05 | R14f, carried over unchanged |
| `planted_signal_auc` | planted signal, each family | mean nested AUC ≥ 0.70 | **proposed**; requires Maintainer confirmation before the first run |
| `planted_signal_p` | planted signal, each family | permutation P ≤ 0.05 in at least 90 % of replicates | **proposed**; requires Maintainer confirmation before the first run |

Seeds: the null scenarios reuse the R14f seeds (auc 600000, perm 700000, leaky 800000). Planted-signal seeds are declared in the same gate file before the first run. Replicate counts are declared, and a reduced count is recorded as an adaptation (E_CALIBRATION_REPLICATES otherwise).

## Error codes introduced by this slice

Classification follows [ADR 0010](../../../docs/adr/0010-adaptive-post-de-policy.md): `validity` refuses only the invalid piece, `policy` is converted to adapt-with-record, `structural` rejects the item it names.

| Code | Packet | Class | Behaviour |
|---|---|---|---|
| `E_WELCH_NOT_DECLARED` | R15a | validity | No Welch result is produced. A limma result is never relabelled as Welch. |
| `E_WELCH_PAIRED_UNIT` | R15a, R15d | validity | The two-sample Welch form is not run on a declared paired design. |
| `W_WELCH_ADAPTED` | R15a | warning | Recorded adaptation for a covariate or block that Welch cannot honour (ADR 0010). |
| `E_LIBRARY_REDISTRIBUTION` | R15b | structural | The artifact or repository would contain a non-redistributable library. Nothing is written. |
| `E_ORA_UNIVERSE` | R15b | validity | ORA on a whole-genome or library-wide background is refused. |
| `E_RANK_DUPLICATE`, `E_RANK_TIE` | R15b | validity | Pre-ranked enrichment on an invalid ranking. |
| `E_ENRICHMENT_FAMILY_POOLED` | R15b | validity | A BH adjustment spans libraries. |
| `E_CALIBRATION_GATE_CHANGED` | R15c | structural | Gate values or hash differ from the pre-registration. |
| `E_CALIBRATION_REPLICATES` | R15c | validity | A reduced replicate count without a recorded adaptation. |
| `E_MANIFEST_HASH` | R15d | structural | Input file hash mismatch. The study is not run. |
| `E_MANIFEST_DESIGN` | R15d | structural | Unknown or unsupported design type. |
| `E_PUBLIC_DATA_UNAPPROVED` | R15d | structural | Public download without a matching approval record. No file is written. |
| `E_ARTIFACT_HASH` | R15d | verify FAIL | An artifact changed after the run. |
| `E_PRIVATE_IDENTIFIER` | R15d | structural | A committed record contains a private identifier. The commit is refused. |
| `E_MAPPING_UNMAPPED` | R15d | validity | An unmapped required column. Only that column is refused. |
| `E_WELCH_BIOLOGICAL_UNIT` | R15a | validity | An unblocked Welch model whose observations repeat a biological unit (technical replicates are not independent, SM03). The stage stops with this reason. |
| `E_WELCH_CONTRAST_NOT_PAIRWISE` | R15a | validity | A declared family member whose contrast is not a pairwise group contrast (its weights are not a +w and -w pair on two group levels). It is kept as a `not_run` endpoint with this reason, so the family stays complete (ADR 0010). |
| `E_WELCH_FAMILY_MIXED` | R15a | validity | A Welch family that also holds endpoints of another engine. Refused, because it would pool Welch and limma P values. |
| `E_WELCH_WEIGHTS_UNSUPPORTED` | R15a | validity | Precision weights on a welch model. The Welch engine takes no weights; the model is INAPPLICABLE. |
| `E_WELCH_MATRIX_UNSUPPORTED` | R15a | validity | A preprocessing sensitivity matrix assigned to a welch model. Refused at configuration time, because the welch engine reads the observed primary matrix only. |
| `E_LIBRARY_FORMAT` | R15b | validity | A malformed GMT file: a row with fewer than three fields, an empty set, or a duplicated set name. The file is not analysed. |
| `E_RANK_INVALID` | R15b | validity | A ranking with a non-finite statistic. The pre-ranked enrichment is refused for that contrast. |
| `E_BIOMARKER_LPO_UNITS` | R15c | validity | A loocv request whose units are not all of one class (within-subject designs): the leave-pair-out AUC is undefined (D-77). The plan adapts the scheme to repeated stratified k-fold and records it; the procedure refuses if called directly. |
| `W_WELCH_ADAPTED` | R15a | warning | A recorded adaptation for a covariate or block that the Welch test cannot honour (run unadjusted, ADR 0010). |

The welch primary and limma co-existence refusal reuses `E_PRIMARY_MODEL_CONFLICT`, and the single-primary rule is enforced by the schema as `E_CONFIG_SCHEMA` ("Too many items").

Existing codes are reused and are not redefined: `E_RESOURCE_HASH`, `E_RESOURCE_DOWNLOAD`, `E_RESOURCE_TAXONOMY`, `E_BIOMARKER_LEAKAGE`, `E_BIOMARKER_GROUP_LEAKAGE`, `E_PHENOTYPE_IMPUTATION`, `E_CONTRAST_NONESTIMABLE`, `E_DEPENDENCY_UNAVAILABLE`. A missing optional package produces NOT_RUN with `E_DEPENDENCY_UNAVAILABLE`, never INAPPLICABLE.

## Claim and method labels

| Label | Meaning |
|---|---|
| `welch` | The Welch engine produced this result (Welch effect, CI and P). |
| `not_multiplicity_controlled` | A paper-compatible selection by raw P. It has no FDR or family control. |
| `exploratory` | Pre-ranked enrichment (exploratory fgsea adapter) or a result from a non-primary library. |
| `calibration_scope` | Calibration applies only to the listed families and designs. |
| `private_status_only` | A private study with only its PASS/FAIL status committed. |

## Outputs

| Module | Location | Files |
|---|---|---|
| Welch (R15a) | `<run>/methods_extension/welch/` | `welch_results.tsv` (effect, CI, t, df, P, BH q, label), `paper_compatible.tsv` (raw-P rows labelled `not_multiplicity_controlled`), `eligibility.json` (adaptations). |
| Pathways (R15b) | `<run>/methods_extension/pathways/` | `libraries.json` (manifest records with hashes), `ora_<library>.tsv`, `prerank_<library>.tsv` (exploratory). |
| Calibration (R15c) | `docs/validation/016-methods-extension-validation/R15c/calibration/` | `gates-preregistration.json`, raw and summary JSON per family and scenario, run log, seeds and versions. |
| Studies (R15d) | `<run>/studies/<study_id>/` (local) and `docs/validation/016-methods-extension-validation/R15d/status/<study_id>.json` (committed) | Local: stage results and artifacts. Committed: slot, design, state and check counts only. |

Every table that feeds a figure has a source table from which each plotted value can be re-derived (SM25).
