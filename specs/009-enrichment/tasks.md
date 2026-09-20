# Tasks: Design-compatible pathways and enrichment

**Input:** [spec.md](spec.md), [plan.md](plan.md), parent shared contracts.
**Roadmap:** R08. **Dependencies:** R06, R07.

## Phase A: Preconditions

Read the live constitution, frozen packet contract, prerequisite evidence/commits and existing interfaces. Maintainer issues one packet starter only after audit, freeze and explicit authorization; parallel work is limited to a frozen disjoint-file group.

## Phase B: Independently testable implementation

- [ ] T071 [US3] Implement and independently test **Pathway null and method dispatcher** in `r/proteomicsCore/R/pathways.R`. Acceptance V071: Plan distinguishes competitive, rotation self-contained, preranked and ORA nulls and rejects unsupported method/design combinations. Evidence: command/output/version/hash in `docs/validation/009-enrichment/V071.json`. Covers FR-071.
- [ ] T072 [US3] Implement and independently test **Independent-design CAMERA** in `r/proteomicsCore/R/pathway_camera.R`. Acceptance V072: Finite gene matrix/design/contrast match official CAMERA with estimated correlation; block/correlation arguments cannot be ignored. Evidence: command/output/version/hash in `docs/validation/009-enrichment/V072.json`. Covers FR-072.
- [ ] T073 [US3] Implement and independently test **Blocked-design ROAST** in `r/proteomicsCore/R/pathway_roast.R`. Acceptance V073: Fixed-subject pairing passes its design only; duplicateCorrelation passes block/correlation without duplicate subject encoding; midpFALSE/rotation settings and reference outputs match. Evidence: command/output/version/hash in `docs/validation/009-enrichment/V073.json`. Covers FR-073.
- [ ] T074 [US3] Implement and independently test **Directional and mixed rotation endpoints** in `r/proteomicsCore/R/pathway_roast.R`. Acceptance V074: Directional and mixed raw P/statistics/families remain distinct; central adjustment uses correct columns rather than package-native FDR mislabeled global. Evidence: command/output/version/hash in `docs/validation/009-enrichment/V074.json`. Covers FR-074.
- [ ] T075 [US3] Implement and independently test **Exploratory fgsea adapter** in `r/proteomicsCore/R/pathway_fgsea.R`. Acceptance V075: Declared native ranks, ties, seeds and eps settings are honored; NES/log2err/warnings/leading edges and nonfinite tests are retained. Evidence: command/output/version/hash in `docs/validation/009-enrichment/V075.json`. Covers FR-075.
- [ ] T076 [US3] Implement and independently test **Exact background-aware ORA** in `r/proteomicsCore/R/pathway_ora.R`. Acceptance V076: Hand-enumerated hypergeometric examples match; all eligible sets enter the family and empty foreground gives a valid explanation. Evidence: command/output/version/hash in `docs/validation/009-enrichment/V076.json`. Covers FR-076.
- [ ] T077 [US3] Implement and independently test **Pathway family adjustment** in `r/proteomicsCore/R/multiplicity.R`. Acceptance V077: Collections/contrasts/null types/primary status obey the frozen families; original and central q fields remain distinguishable. Evidence: command/output/version/hash in `docs/validation/009-enrichment/V077.json`. Covers FR-077.
- [ ] T078 [US3] Implement and independently test **Alternative-engine pathway sensitivity** in `r/proteomicsCore/R/pathways.R`. Acceptance V078: proDA/DEqMS ranks and linear gene-model sensitivity are labeled correctly and never sold as equivalent correlation-adjusted primary-engine inference. Evidence: command/output/version/hash in `docs/validation/009-enrichment/V078.json`. Covers FR-078.
- [ ] T079 [US3] Implement and independently test **Leading-edge redundancy and diagnostics** in `r/proteomicsCore/R/pathway_summary.R`. Acceptance V079: Overlap relationships preserve original set IDs/P values and do not merge independent discoveries or invent causal mechanisms. Evidence: command/output/version/hash in `docs/validation/009-enrichment/V079.json`. Covers FR-079.
- [ ] T080 [US3] Implement and independently test **Enrichment golden and failure suite** in `tests/scientific/test_pathways.py; r/proteomicsCore/tests/testthat/test-pathways.R`. Acceptance V080: Exact ORA, reference CAMERA/ROAST/fgsea, missing snapshots, NA/Inf matrix, all nonfinite tests and no discoveries exercise honest statuses. Evidence: command/output/version/hash in `docs/validation/009-enrichment/V080.json`. Covers FR-080.

## Phase C: Review checkpoint

All checkboxes start unchecked. Complete a checkbox only when Maintainer has read the diff, rerun the relevant tests and recorded acceptance evidence plus commit. A task can span multiple Packet Implementer correction tasks. If a runtime/resource is unavailable, keep its verification pending and continue independent work; do not write a successful fixture-shaped substitute.

## Dependencies & Execution Order

The order above follows shared contract → implementation → integration → adversarial verification. Preserve data-model and method boundaries. Resolve packet prerequisites before dependent tasks and obey the exact allowlist in the frozen packet index. A thin vertical demonstration is a checkpoint, not the end of the parent roadmap.
