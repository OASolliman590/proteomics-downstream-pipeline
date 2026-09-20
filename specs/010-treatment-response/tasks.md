# Tasks: Treatment response, equivalence and independent scores

**Input:** [spec.md](spec.md), [plan.md](plan.md), parent shared contracts.
**Roadmap:** R09. **Dependencies:** R05.

## Phase A: Preconditions

Read the live constitution, frozen packet contract, prerequisite evidence/commits and existing interfaces. Maintainer issues one packet starter only after audit, freeze and explicit authorization; parallel work is limited to a frozen disjoint-file group.

## Phase B: Independently testable implementation

- [ ] T081 [US4] Implement and independently test **Explicit disease-treatment-residual axes** in `r/proteomicsCore/R/response.R`. Acceptance V081: d, t, r and interaction estimates align to planned contrasts with covariance; r=d+t matches analytic fixtures. Evidence: command/output/version/hash in `docs/validation/010-treatment-response/V081.json`. Covers FR-081.
- [ ] T082 [US4] Implement and independently test **Descriptive reversal without double dipping** in `r/proteomicsCore/R/response.R`. Acceptance V082: Current-data selected features can be displayed but produce no confirmatory score P or rescue label; selection provenance is explicit. Evidence: command/output/version/hash in `docs/validation/010-treatment-response/V082.json`. Covers FR-082.
- [ ] T083 [US4] Implement and independently test **Bounded response categories** in `r/proteomicsCore/R/response.R`. Acceptance V083: Configured epsilon classes handle zero/small d, RI=1, RI>1+epsilon and RI=3; old unlimited full-reversal label is absent. Evidence: command/output/version/hash in `docs/validation/010-treatment-response/V083.json`. Covers FR-083.
- [ ] T084 [US4] Implement and independently test **Ratio uncertainty contract** in `r/proteomicsCore/R/response_uncertainty.R`. Acceptance V084: Supported joint-covariance method handles unstable/unbounded RI intervals; dividing CI endpoints is rejected and descriptive-only mode is explicit. Evidence: command/output/version/hash in `docs/validation/010-treatment-response/V084.json`. Covers FR-084.
- [ ] T085 [US4] Implement and independently test **Model-correct residual equivalence** in `r/proteomicsCore/R/equivalence.R`. Acceptance V085: Analytic TOST examples match both one-sided tests/max-P and 90% equivalence intervals; a nonsignificant zero-null test alone cannot pass. Evidence: command/output/version/hash in `docs/validation/010-treatment-response/V085.json`. Covers FR-085.
- [ ] T086 [US4] Implement and independently test **Conjunction claim eligibility** in `r/proteomicsCore/R/equivalence.R`. Acceptance V086: Only independently defined direction and valid model-supported conjunction yield adjusted restoration endpoints; FDR-list intersection is not labeled rescue FDR. Evidence: command/output/version/hash in `docs/validation/010-treatment-response/V086.json`. Covers FR-086.
- [ ] T087 [US4] Implement and independently test **Independent fixed-score artifacts** in `schemas/independent-score.schema.json; r/proteomicsCore/R/scores.R`. Acceptance V087: Training identity/weights/transforms are frozen and checked; validation-outcome-derived selection or refitting produces descriptive-only eligibility. Evidence: command/output/version/hash in `docs/validation/010-treatment-response/V087.json`. Covers FR-087.
- [ ] T088 [US4] Implement and independently test **Exact and Monte Carlo randomization** in `r/proteomicsCore/R/randomization.R`. Acceptance V088: All 70 four-versus-four allocations use k/N; Monte Carlo +1 and interval are distinct; paired/restricted permutations honor blocks and ties. Evidence: command/output/version/hash in `docs/validation/010-treatment-response/V088.json`. Covers FR-088.
- [ ] T089 [US4] Implement and independently test **Score/equivalence reporting semantics** in `src/proteomics_pipeline/reporting/response_data.py`. Acceptance V089: Report distinguishes movement, overshoot, equivalence, conjunction evidence and unknown training overlap; no causal rescue prose is generated. Evidence: command/output/version/hash in `docs/validation/010-treatment-response/V089.json`. Covers FR-089.
- [ ] T090 [US4] Implement and independently test **Circularity/null/overshoot tests** in `tests/scientific/test_response.py; r/proteomicsCore/tests/testthat/test-response.R`. Acceptance V090: Shared-untreated null exposes negative effect covariance but cannot yield a spurious inferential score path; exact old examples and RI overshoot regressions pass. Evidence: command/output/version/hash in `docs/validation/010-treatment-response/V090.json`. Covers FR-090.

## Phase C: Review checkpoint

All checkboxes start unchecked. Complete a checkbox only when Maintainer has read the diff, rerun the relevant tests and recorded acceptance evidence plus commit. A task can span multiple Packet Implementer correction tasks. If a runtime/resource is unavailable, keep its verification pending and continue independent work; do not write a successful fixture-shaped substitute.

## Dependencies & Execution Order

The order above follows shared contract → implementation → integration → adversarial verification. Preserve data-model and method boundaries. Resolve packet prerequisites before dependent tasks and obey the exact allowlist in the frozen packet index. A thin vertical demonstration is a checkpoint, not the end of the parent roadmap.
