# Tasks: Preprocessing, missingness and quality control

**Input:** [spec.md](spec.md), [plan.md](plan.md), parent shared contracts.
**Roadmap:** R03. **Dependencies:** R02.

## Phase A: Preconditions

Read the live constitution, frozen packet contract, prerequisite evidence/commits and existing interfaces. Maintainer issues one packet starter only after audit, freeze and explicit authorization; parallel work is limited to a frozen disjoint-file group.

## Phase B: Independently testable implementation

- [ ] T021 [US1] Implement and independently test **Preserve processed abundance by default** in `r/proteomicsCore/R/preprocess.R`. Acceptance V021: An already normalized log2 fixture remains unchanged; every optional transform has a new artifact and recorded policy. Evidence: command/output/version/hash in `docs/validation/004-preprocessing-qc/V021.json`. Covers FR-021.
- [ ] T022 [US1] Implement and independently test **LFQ normalization policies** in `r/proteomicsCore/R/normalization.R`. Acceptance V022: Median/reference normalization factors match analytic fixtures; quantile is explicit sensitivity; missing reference coverage is not silently ignored. Evidence: command/output/version/hash in `docs/validation/004-preprocessing-qc/V022.json`. Covers FR-022.
- [ ] T023 [US1] Implement and independently test **TMT plex and bridge policies** in `r/proteomicsCore/R/tmt.R`. Acceptance V023: Two-plex synthetic bridges align loading offsets; missing bridge blocks bridge mode, while identifiable no-bridge plex-covariate inference is separately supported and confounding fails. Evidence: command/output/version/hash in `docs/validation/004-preprocessing-qc/V023.json`. Covers FR-023.
- [ ] T024 [US1] Implement and independently test **Contrast-aware coverage masks** in `r/proteomicsCore/R/filtering.R`. Acceptance V024: Original observed coverage is computed at biological-unit grain; ordinary available-case all-missing groups are nonestimable while eligible proDA follows native dropout rules with explicit uncertainty. Evidence: command/output/version/hash in `docs/validation/004-preprocessing-qc/V024.json`. Covers FR-024.
- [ ] T025 [US1] Implement and independently test **Missingness mechanism diagnostics** in `r/proteomicsCore/R/missingness.R`. Acceptance V025: Observation/group/abundance summaries distinguish observed versus previously imputed values and do not assert MCAR/MNAR from a threshold. Evidence: command/output/version/hash in `docs/validation/004-preprocessing-qc/V025.json`. Covers FR-025.
- [ ] T026 [US1] Implement and independently test **PCA and correlation diagnostics** in `r/proteomicsCore/R/qc.R`. Acceptance V026: PCA-only fill is isolated; constant/empty/small matrices return valid diagnostic states; variance and feature/sample identities match direct calculations. Evidence: command/output/version/hash in `docs/validation/004-preprocessing-qc/V026.json`. Covers FR-026.
- [ ] T027 [US1] Implement and independently test **Prespecified exclusions and watchlists** in `src/proteomics_pipeline/planning.py; r/proteomicsCore/R/qc.R`. Acceptance V027: Outlier flags never auto-delete observations; a changed exclusion reason/policy creates a new plan and sensitivity lineage. Evidence: command/output/version/hash in `docs/validation/004-preprocessing-qc/V027.json`. Covers FR-027.
- [ ] T028 [US1] Implement and independently test **Explicit imputation sensitivities** in `r/proteomicsCore/R/sensitivity.R`. Acceptance V028: Deterministic minimum, stochastic Gaussian and KNN have correct names/seeds/masks, remain secondary and never alter the primary matrix. Evidence: command/output/version/hash in `docs/validation/004-preprocessing-qc/V028.json`. Covers FR-028.
- [ ] T029 [US1] Implement and independently test **Detection-only exploratory endpoint** in `r/proteomicsCore/R/detection.R`. Acceptance V029: Independent two-group detection uses exact 2x2 counts and a separate BH family; paired/complex detection requests are marked unsupported. Evidence: command/output/version/hash in `docs/validation/004-preprocessing-qc/V029.json`. Covers FR-029.
- [ ] T030 [US1] Implement and independently test **QC source data and stage report** in `src/proteomics_pipeline/reporting/qc_data.py; tests/integration/test_qc.py`. Acceptance V030: All plots have exact source tables; QC-only run completes without R modeling or any significant features. Evidence: command/output/version/hash in `docs/validation/004-preprocessing-qc/V030.json`. Covers FR-030.

## Phase C: Review checkpoint

All checkboxes start unchecked. Complete a checkbox only when Maintainer has read the diff, rerun the relevant tests and recorded acceptance evidence plus commit. A task can span multiple Packet Implementer correction tasks. If a runtime/resource is unavailable, keep its verification pending and continue independent work; do not write a successful fixture-shaped substitute.

## Dependencies & Execution Order

The order above follows shared contract → implementation → integration → adversarial verification. Preserve data-model and method boundaries. Resolve packet prerequisites before dependent tasks and obey the exact allowlist in the frozen packet index. A thin vertical demonstration is a checkpoint, not the end of the parent roadmap.
