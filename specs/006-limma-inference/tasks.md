# Tasks: Core limma inference and multiplicity

**Input:** [spec.md](spec.md), [plan.md](plan.md), parent shared contracts.
**Roadmap:** R05. **Dependencies:** R04.

## Phase A: Preconditions

Read the live constitution, frozen packet contract, prerequisite evidence/commits and existing interfaces. Maintainer issues one packet starter only after audit, freeze and explicit authorization; parallel work is limited to a frozen disjoint-file group.

## Phase B: Independently testable implementation

- [ ] T041 [US2] Implement and independently test **Observed-data limma adapter** in `r/proteomicsCore/R/model_limma.R`. Acceptance V041: Eligible fixtures match direct lmFit/eBayes calls for effects/SE/df/P; primary input retains NA and declared normalization. Evidence: command/output/version/hash in `docs/validation/006-limma-inference/V041.json`. Covers FR-041.
- [ ] T042 [US2] Implement and independently test **Exact sparse and weighted contrasts** in `r/proteomicsCore/R/model_limma.R`. Acceptance V042: General designs use validated contrast reparametrization/direct covariance; nonorthogonal sparse example matches a separate oracle. Evidence: command/output/version/hash in `docs/validation/006-limma-inference/V042.json`. Covers FR-042.
- [ ] T043 [US2] Implement and independently test **Robust/trend diagnostics** in `r/proteomicsCore/R/diagnostics.R`. Acceptance V043: Settings, prior/posterior variance and trend diagnostics are saved; robust eBayes is not mislabeled robust sample regression. Evidence: command/output/version/hash in `docs/validation/006-limma-inference/V043.json`. Covers FR-043.
- [ ] T044 [US2] Implement and independently test **Effect-threshold TREAT endpoint** in `r/proteomicsCore/R/inference.R`. Acceptance V044: Boundary effects compare to direct TREAT; true threshold-null evidence is distinct from observed effect filtering and zero-null results. Evidence: command/output/version/hash in `docs/validation/006-limma-inference/V044.json`. Covers FR-044.
- [ ] T045 [US2] Implement and independently test **Backend-correct confidence intervals** in `r/proteomicsCore/R/inference.R`. Acceptance V045: Intervals use correct moderated variance/df, match independent calculations and preserve unbounded/unavailable cases with reason. Evidence: command/output/version/hash in `docs/validation/006-limma-inference/V045.json`. Covers FR-045.
- [ ] T046 [US2] Implement and independently test **Declared multiple-testing families** in `r/proteomicsCore/R/multiplicity.R`. Acceptance V046: Finite-P BH/BY matches analytic and R references; separate endpoints/primary/sensitivity sets do not accidentally pool or omit rows. Evidence: command/output/version/hash in `docs/validation/006-limma-inference/V046.json`. Covers FR-046.
- [ ] T047 [US2] Implement and independently test **Omnibus and interaction outputs** in `r/proteomicsCore/R/inference.R`. Acceptance V047: F tests have distinct hypotheses/families, complete planned contrasts are exported and no unadjusted post-hoc claims appear. Evidence: command/output/version/hash in `docs/validation/006-limma-inference/V047.json`. Covers FR-047.
- [ ] T048 [US2] Implement and independently test **Complete typed differential tables** in `schemas/differential-result.schema.json; r/proteomicsCore/R/export.R`. Acceptance V048: Every planned feature-contrast row is tested or explicitly excluded/nonestimable; no zero/NA ambiguity or mislabeled q fields. Evidence: command/output/version/hash in `docs/validation/006-limma-inference/V048.json`. Covers FR-048.
- [ ] T049 [US2] Implement and independently test **Biological-unit influence diagnostics** in `r/proteomicsCore/R/influence.R`. Acceptance V049: Whole subjects/technical sets are omitted appropriately; relevant omission counts are correct and influence stays descriptive. Evidence: command/output/version/hash in `docs/validation/006-limma-inference/V049.json`. Covers FR-049.
- [ ] T050 [US2] Implement and independently test **Core inference golden suite** in `tests/scientific/; r/proteomicsCore/tests/testthat/test-limma.R`. Acceptance V050: Complete/sparse/paired/blocked fixtures match pinned direct references and no-null/no-discovery result paths pass. Evidence: command/output/version/hash in `docs/validation/006-limma-inference/V050.json`. Covers FR-050.

## Phase C: Review checkpoint

All checkboxes start unchecked. Complete a checkbox only when Maintainer has read the diff, rerun the relevant tests and recorded acceptance evidence plus commit. A task can span multiple Packet Implementer correction tasks. If a runtime/resource is unavailable, keep its verification pending and continue independent work; do not write a successful fixture-shaped substitute.

## Dependencies & Execution Order

The order above follows shared contract → implementation → integration → adversarial verification. Preserve data-model and method boundaries. Resolve packet prerequisites before dependent tasks and obey the exact allowlist in the frozen packet index. A thin vertical demonstration is a checkpoint, not the end of the parent roadmap.
