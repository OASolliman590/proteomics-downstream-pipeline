# Release phases — normative scope

Kit: **1.2.0-frozen**. No maintained phase has shipped. R01 is the only next packet eligible for explicit authorization.

| Phase | Packets | Release name | Meaning |
|---|---|---|---|
| 0 | Recovered baseline | baseline-recovered | Inspectable study scripts and immutable audit, not a validated general platform. |
| 1 | R01–R05 + R10a | v0.1-limma-core | CLI, all assigned intake/QC/design/limma requirements and a thin offline HTML report. |
| 2 | R06, R07, R08, R09 | v0.2-qualified-methods | Qualified DEqMS/proDA, frozen mapping/resources, design-valid enrichment and response. |
| 3 | R10b, R11, R12 | v1.0-defensible | Full report, environment locks, calibration, offline reproduction and legacy reconciliation. |

R10 remains one of twelve slices. R10a owns FR-091–FR-094 / T091–T094 / V091–V094. R10b owns FR-095–FR-100 / T095–T100 / V095–V100. No identifier is added, retired or renumbered. Phase 1 does not mark R10 as a whole complete; Phase 3 re-exercises the R10a interface without changing its ownership.

Phase 1 requires all Phase 1 acceptance cases and the cross-cutting scenarios in [validation-strategy.md](validation-strategy.md). This includes qualified TMT/vendor intake and other scope already assigned to R02–R05; “limma core” is not permission to silently remove those requirements. R06-specific numerical behavior remains Phase 2. No Phase 1 score module or score P-value is exposed. An explicit QC-only execution may accept unknown scale; an inference run may not.

Phase 1 records the actual solved environment and runs its direct numerical oracles. Fully reproducible lock restoration and 1,000-dataset release calibration remain R11, not requirements for calling Phase 1 a milestone. Phase 1 MUST NOT claim completed R11 validation. Production in every implemented phase still uses local hashed resources only.

The only permitted concurrent group is R06 ∥ R07 ∥ R09 after R05 and R10a have both been accepted. R10a runs serially after R05. R08 waits for acceptance of R06, R07 and R09 and is scheduled serially after that whole-group barrier. R10b waits for R10a, R08 and R09. R11 follows R10b; R12 follows R11. Exact ownership and dependencies are in [packet-index.md](packet-index.md).

A phase gate records executed scope, reviewed tree hash and evidence IDs. Missing optional study prerequisites are INAPPLICABLE or NOT_RUN with reasons; absent software is not scientific inapplicability. A failed required gate blocks the corresponding milestone. License/disclosure authorization is a separate unresolved public-release gate, not permission to invent terms. Private study validation cannot be claimed until the Maintainer actually performs it.
