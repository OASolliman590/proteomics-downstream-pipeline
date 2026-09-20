# Roadmap: Streamlined downstream proteomics

Twelve Spec Kit packets turn the recovered study code into a maintained pipeline. Each retains its own spec.md, plan.md and tasks.md; the corrected ownership contract is in [packet-index.md](packet-index.md). All packets are in scope. Status reflects verified implementation, not specification completeness.

| ID | Slice | Intent / boundary | Dependencies | Status | Specification |
|---|---|---|---|---|---|
| R01 | Runtime, package and contract foundation | Assigned capabilities only; shared contracts bind downstream interfaces | none | draft / unverified prior-worktree candidate | [002-foundation](../002-foundation/spec.md) |
| R02 | Canonical protein input and biological identity | Assigned capabilities only; shared contracts bind downstream interfaces | R01 | draft | [003-intake](../003-intake/spec.md) |
| R03 | Preprocessing, missingness and quality control | Assigned capabilities only; shared contracts bind downstream interfaces | R02 | draft | [004-preprocessing-qc](../004-preprocessing-qc/spec.md) |
| R04 | Design validation, exact contrasts and blocking | Assigned capabilities only; shared contracts bind downstream interfaces | R03 | draft | [005-design-contrasts](../005-design-contrasts/spec.md) |
| R05 | Core limma inference and multiplicity | Assigned capabilities only; shared contracts bind downstream interfaces | R04 | draft | [006-limma-inference](../006-limma-inference/spec.md) |
| R06 | Assay-qualified DEqMS and proDA backends | Assigned capabilities only; shared contracts bind downstream interfaces | R05 | draft | [007-assay-engines](../007-assay-engines/spec.md) |
| R07 | Versioned annotation, protein groups and gene sets | Assigned capabilities only; shared contracts bind downstream interfaces | R05 | draft | [008-resources-mapping](../008-resources-mapping/spec.md) |
| R08 | Design-compatible pathways and enrichment | Assigned capabilities only; shared contracts bind downstream interfaces | R06, R07 | draft | [009-enrichment](../009-enrichment/spec.md) |
| R09 | Treatment response, equivalence and independent scores | Assigned capabilities only; shared contracts bind downstream interfaces | R05 | draft | [010-treatment-response](../010-treatment-response/spec.md) |
| R10 | Unified offline report and streamlined workflow | Assigned capabilities only; shared contracts bind downstream interfaces | R03, R06, R08, R09 | draft | [011-reporting](../011-reporting/spec.md) |
| R11 | Reproducibility, calibration and continuous validation | Assigned capabilities only; shared contracts bind downstream interfaces | R10 | draft | [012-validation](../012-validation/spec.md) |
| R12 | Archived-study regression, documentation and versioned successor | Assigned capabilities only; shared contracts bind downstream interfaces | R11 | draft | [013-release](../013-release/spec.md) |

The parent specification defines full downstream scope. Each sub-spec points back to its immutable roadmap ID. A boundary change updates this roadmap, dependent specs, contracts, packet index and traceability before dependent work. The only draft parallel group is R06 + R07 + R09 after R05; every other edge is serial. No packet is dispatchable before Independent Reviewer audit, Maintainer freeze and explicit authorization.
