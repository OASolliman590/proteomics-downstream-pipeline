# Packet index — v1.2.0-frozen

Kit frozen 2026-09-20. R01 accepted 2026-09-21. On 2026-10-01 the operator gave GO for R02 onward in this order and authorized R13 (PERMANOVA, [ADR 0006](../../docs/adr/0006-permanova-scope-amendment.md)) with the Claude route ([ADR 0007](../../docs/adr/0007-claude-route-and-integration-amendments.md)). R02–R11 and R13 pass their gates (V107: CI run 37087074374 on fa4d948) and are merged into `main`. R14a–R14f (ADR 0009) pass their gates on branch `claude/post-de`, independent audit NOT_RUN. R12: V115–V117 PASS; the private gates (V111, V112, V114), V113, the release/publication gate V118, V119 and V120 are NOT_RUN (Maintainer-only or operator decisions).

[PHASES.md](PHASES.md) defines release scope. [packet-ownership.json](packet-ownership.json) is the authoritative exact write allowlist and dispatch-dependency record. All listed source/test paths are pending creation after freeze. All paths not assigned to a packet are forbidden to that implementer. There are twelve original slices plus the amended R13 slice, and fourteen dispatch units because R10 is split without changing any FR/T/V identity.

| Dispatch unit | Phase | Dependency | Acceptance | Risk | Slice |
|---|---|---|---|---|---|
| R01 | 1 | Frozen baseline/contracts | V001–V010 | high | [Runtime, package and contract foundation](../002-foundation/spec.md) |
| R02 | 1 | R01 | V011–V020 | high | [Canonical protein input and biological identity](../003-intake/spec.md) |
| R03 | 1 | R02 | V021–V030 | high | [Preprocessing, missingness and quality control](../004-preprocessing-qc/spec.md) |
| R04 | 1 | R03 | V031–V040 | high | [Design validation, exact contrasts and blocking](../005-design-contrasts/spec.md) |
| R05 | 1 | R04 | V041–V050 | high | [Core limma inference and multiplicity](../006-limma-inference/spec.md) |
| R10a | 1 | R05 | V091–V094 | high | [Unified offline report and streamlined workflow](../011-reporting/spec.md) |
| R13 | 2 | R05, R10a | V121–V130 | high | [Multivariate PERMANOVA and dispersion testing](../014-multivariate-permanova/spec.md) (amendment) |
| R06 | 2 | R05, R10a | V051–V060 | high | [Assay-qualified DEqMS and proDA backends](../007-assay-engines/spec.md) |
| R07 | 2 | R05, R10a | V061–V070 | high | [Versioned annotation, protein groups and gene sets](../008-resources-mapping/spec.md) |
| R09 | 2 | R05, R10a | V081–V090 | high | [Treatment response, equivalence and independent scores](../010-treatment-response/spec.md) |
| R08 | 2 | R06, R07, R09 | V071–V080 | high | [Design-compatible pathways and enrichment](../009-enrichment/spec.md) |
| R10b | 3 | R10a, R08, R09 | V095–V100 | high | [Unified offline report and streamlined workflow](../011-reporting/spec.md) |
| R11 | 3 | R10b | V101–V110 | high | [Reproducibility, calibration and continuous validation](../012-validation/spec.md) |
| R12 | 3 | R11 | V111–V120 | high | [Archived-study regression, documentation and versioned successor](../013-release/spec.md) |
| R14a | 4 | R11, R13 | V131–V136 | high | [Post-DE: result structure and concordance](../015-post-de-analysis/spec.md) (amendment) |
| R14b | 4 | R14a | V137–V143 | high | [Post-DE: robustness and sensitivity](../015-post-de-analysis/spec.md) (amendment) |
| R14c | 4 | R14b | V144–V148 | high | [Post-DE: protein–phenotype association](../015-post-de-analysis/spec.md) (amendment) |
| R14d | 4 | R14c, R09 | V149–V160 | high | [Post-DE: biomarker discrimination evaluation](../015-post-de-analysis/spec.md) (amendment) |
| R14e | 4 | R14d, R07 | V161–V165 | high | [Post-DE: co-abundance modules and networks](../015-post-de-analysis/spec.md) (amendment) |
| R14f | 4 | R14a–R14e, R10b | V166–V167 | high | [Post-DE: eligibility reporting and generality matrix](../015-post-de-analysis/spec.md) (amendment) |

## Ownership and scheduling

Every packet forbids pipeline/**, legacy/**, docs/audit/**, private_evidence/**, original archive/manifest, credentials and unrelated/global configuration. Implementers also cannot edit .specify/**, specs/**, commit, push, select a license or receive private study inputs/logs. The Maintainer alone edits kit/progress/traceability and runs private V111/V114; those actions are not concealed in a packet worker allowlist.

All write allowlists are disjoint, including serial packets. There is no R01 R/** wildcard, shared planning.py ownership or later reopening of runtime.py/multiplicity.R. The fixed module/file handoff is specified in the [runtime extension seam](contracts/cli-and-artifacts.md#runtime-extension-seam). Each packet uses its own fixture/test/evidence subdirectories. A needed interface change returns to the Maintainer before dispatch; it is not a silent allowlist exception.

The only permitted concurrent group is **R06 ∥ R07 ∥ R09 after accepted R05 and accepted R10a**. The machine-readable dependency record includes both prerequisites for every member. R08 runs serially only after all three group members are accepted; R09 is therefore an explicit dispatch dependency even though R08's scientific inputs come from R06/R07. R10b follows R10a/R08/R09; R11 and R12 follow serially. No other packet, including reporting/documentation, runs concurrently with this group. All twelve slices remain necessary for SYS-03/full v1.0; the slice-completion graph in roadmap.json is not a substitute for this subpacket dispatch order.

## Minimal handoff

A packet receipt identifies its frozen tree/contract versions, changed paths, acceptance IDs, executed commands and actual outputs/failures. The Maintainer reviews the diff and independently runs the gates before recording completion. Evidence status remains pending-after-freeze or NOT_RUN when no verified execution exists. Independent Reviewer audit does not authorize its author to mark their own candidate frozen.

Current next action: the Maintainer's review of `claude/post-de`, an independent audit of R14a–R14f, CI on the pushed branch, and the Maintainer-only R12 steps in [the runbook](../../docs/validation/013-release/MAINTAINER_RUNBOOK.md). Shared-interface amendments made by the Maintainer route are listed in `maintainer_amendments` of [packet-ownership.json](packet-ownership.json). No lockfile or validated general platform is implied.

Phase 4 amendment (2026-10-02, [ADR 0009](../../docs/adr/0009-post-de-scope-amendment.md)): slice 015 adds dispatch units R14a–R14f (FR/T/V-131–167), serial after accepted R11 and the R13 audit fixes. Implemented 2026-10-04 on `claude/post-de` (see the receipts). The shared-interface amendment (capability maps, optional `post_de` schema block, R Suggests) is made by the Maintainer before R14a, as described in [the slice plan](../015-post-de-analysis/plan.md).
