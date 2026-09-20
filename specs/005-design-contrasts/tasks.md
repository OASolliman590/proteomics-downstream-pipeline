# Tasks: Design validation, exact contrasts and blocking

**Input:** [spec.md](spec.md), [plan.md](plan.md), parent shared contracts.
**Roadmap:** R04. **Dependencies:** R03.

## Phase A: Preconditions

Read the live constitution, frozen packet contract, prerequisite evidence/commits and existing interfaces. Maintainer issues one packet starter only after audit, freeze and explicit authorization; parallel work is limited to a frozen disjoint-file group.

## Phase B: Independently testable implementation

- [ ] T031 [US2] Implement and independently test **Safe declarative design grammar** in `src/proteomics_pipeline/planning.py; r/proteomicsCore/R/design.R`. Acceptance V031: Groups, continuous covariates and interactions generate expected matrices; arbitrary function calls/code strings fail before R evaluation. Evidence: command/output/version/hash in `docs/validation/005-design-contrasts/V031.json`. Covers FR-031.
- [ ] T032 [US2] Implement and independently test **Stable levels and aligned matrices** in `r/proteomicsCore/R/design.R`. Acceptance V032: Sample order shuffling with metadata alignment preserves coefficients; absent references and missing covariates give explicit errors. Evidence: command/output/version/hash in `docs/validation/005-design-contrasts/V032.json`. Covers FR-032.
- [ ] T033 [US2] Implement and independently test **Rank and confounding diagnostics** in `r/proteomicsCore/R/estimability.R`. Acceptance V033: Batch=treatment fixture is detected, with aliased terms; no pseudoinverse silently manufactures an interpretable treatment effect. Evidence: command/output/version/hash in `docs/validation/005-design-contrasts/V033.json`. Covers FR-033.
- [ ] T034 [US2] Implement and independently test **Independent, paired and repeated plans** in `r/proteomicsCore/R/design.R`. Acceptance V034: Subject-fixed and eligible duplicateCorrelation designs preserve the experimental unit and reject invalid/unidentified blocks. Evidence: command/output/version/hash in `docs/validation/005-design-contrasts/V034.json`. Covers FR-034.
- [ ] T035 [US2] Implement and independently test **Numeric contrasts and interaction semantics** in `r/proteomicsCore/R/contrasts.R`. Acceptance V035: Known disease/treatment/residual/interaction coefficients match analytic expectations; significance difference is never used as interaction. Evidence: command/output/version/hash in `docs/validation/005-design-contrasts/V035.json`. Covers FR-035.
- [ ] T036 [US2] Implement and independently test **Featurewise estimability and df** in `r/proteomicsCore/R/estimability.R`. Acceptance V036: Sparse features get engine-specific n/df/reason records; ordinary available-case all-missing groups are nonestimable and eligible native-dropout estimates retain their prior-dependence caveat. Evidence: command/output/version/hash in `docs/validation/005-design-contrasts/V036.json`. Covers FR-036.
- [ ] T037 [US2] Implement and independently test **Exact contrast covariance strategy** in `r/proteomicsCore/R/contrasts.R`. Acceptance V037: Sparse nonorthogonal weighted designs match direct coefficient-refit SE/P; exact shortcut is restricted to proven eligible designs. Evidence: command/output/version/hash in `docs/validation/005-design-contrasts/V037.json`. Covers FR-037.
- [ ] T038 [US2] Implement and independently test **Method-design applicability registry** in `r/proteomicsCore/R/capabilities.R; src/proteomics_pipeline/planning.py`. Acceptance V038: Unsupported proDA/DEqMS block or threshold combinations fail explicitly and no engine is silently substituted. Evidence: command/output/version/hash in `docs/validation/005-design-contrasts/V038.json`. Covers FR-038.
- [ ] T039 [US2] Implement and independently test **Frozen analysis plan and matrices** in `src/proteomics_pipeline/planning.py`. Acceptance V039: Plan saves configuration, inclusion masks, matrices, hypotheses, families and hashes before fitting; config changes invalidate the plan. Evidence: command/output/version/hash in `docs/validation/005-design-contrasts/V039.json`. Covers FR-039.
- [ ] T040 [US2] Implement and independently test **Design adversarial verification** in `tests/scientific/test_design_contract.py; r/proteomicsCore/tests/testthat/test-design.R`. Acceptance V040: Unequal n, singletons, missing factor levels, repeated-unit imbalance and reordered IDs have reference-backed expected behavior. Evidence: command/output/version/hash in `docs/validation/005-design-contrasts/V040.json`. Covers FR-040.

## Phase C: Review checkpoint

All checkboxes start unchecked. Complete a checkbox only when Maintainer has read the diff, rerun the relevant tests and recorded acceptance evidence plus commit. A task can span multiple Packet Implementer correction tasks. If a runtime/resource is unavailable, keep its verification pending and continue independent work; do not write a successful fixture-shaped substitute.

## Dependencies & Execution Order

The order above follows shared contract → implementation → integration → adversarial verification. Preserve data-model and method boundaries. Resolve packet prerequisites before dependent tasks and obey the exact allowlist in the frozen packet index. A thin vertical demonstration is a checkpoint, not the end of the parent roadmap.
