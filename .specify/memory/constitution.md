<!-- Sync impact: 1.1.0 introduces the reviewed Maintainer/Independent Reviewer/Packet Implementer workflow; scientific and evidence-preservation principles are unchanged. All feature plans must read this live. -->
# Downstream Proteomics Constitution

Version: 1.1.0. Ratified: 2026-09-12. Last amended: 2026-09-20.

## I. Experimental evidence precedes statistical inference

The system MUST preserve the experimental unit, input scale, assay type, sample/feature identity and provenance. It MUST distinguish biological specimens, technical injections and multiplex reference channels. Missing tissue or preprocessing details MUST remain missing and constrain claims. Input evidence cannot be upgraded by a hardcoded PASS or a plausible inferred label. Raw-MS identification, peptide identification FDR and tissue mechanisms are outside the evidence supplied by a processed abundance table.

## II. An analysis plan is an executable scientific contract

Every run MUST freeze validated configuration, model/design, contrasts, missingness/filter policy, primary method, gene mapping, resource hashes and multiplicity families before fitting. Parameter changes create a new plan hash. No automatic choice may maximize significant hits. Method eligibility MUST be checked explicitly. A statistically inapplicable method MUST fail with a useful reason rather than silently changing engines or filling data. A null or inconclusive result is a successful analysis outcome.

## III. Biological uncertainty must survive computation

Effect estimates, uncertainty intervals, raw P values, adjusted values, estimability and missing-data diagnostics MUST retain their meaning and provenance. Zero-null, effect-threshold, equivalence, competitive enrichment and self-contained enrichment hypotheses MUST be different result types. Sign opposition, nonsignificance and pathway labels cannot substitute for evidence of restoration, equivalence or causality. Protein selection and testing MUST not reuse labels in an invalid inferential procedure.

## IV. Reproducibility is verified, not asserted

Each stage MUST emit content hashes, software/resource versions, RNG state/seed, input/output lineage, warnings and a status derived from executed work. Failures cannot become empty healthy outputs. Production runs MUST not auto-install dependencies or silently fetch unversioned resources. Scientific output tables MUST be reproducible under a pinned environment to a documented tolerance; cosmetic timestamps are excluded from semantic comparisons. Clean local and CI executions, independent numerical oracles, synthetic calibration and the archived study regression are separate gates.

## V. The historical baseline remains inspectable

The recovered scripts, audit and 165 source evidence files are an immutable reference. New implementation MUST live in a maintained package rather than rewriting the historical scripts. Corrected results MUST identify intentional changes and preserve the old result for comparison. The auxiliary Negr1 result and old permutation discrepancy are regression cases, not constants to bake into production logic. Existing tests may change only with a documented contract correction and replacement coverage.

## VI. Spec Kit is the development record

Every development slice MUST have a specification, implementation plan, uniquely numbered tasks, acceptance tests and an implementation brief linked from the roadmap. Requirement → task → test → evidence links MUST remain complete. Tasks become complete only after the orchestrator verifies code and gates. Material design discoveries update the affected contracts/spec/plan before the next dependent dispatch, with an ADR explaining the change. Routine implementation decisions within this scope need no repeated user confirmation.

## VII. Small, coherent interfaces and real implementations

The user experience MUST be a validated configuration and a single downstream run command with inspectable intermediate results. The implementation MUST use a small Python orchestration layer and a cohesive R statistics package with explicit file contracts. Avoid a web service, database, workflow cluster or plugin framework unless a demonstrated requirement needs it. No placeholder backend, fixture-driven production result, unsupported optional flag or silent skip may count as implementation. Optional means scientifically inapplicable for some studies, not unimplemented.

## VIII. Delegated implementation has an accountable reviewer

The Maintainer owns the kit, freeze corrections, scientific review, independent gate execution, commits, private regression and loop closure. The Independent Reviewer performs specification-freeze and risky-diff audits; that role never implements, receives private test logs or closes work. The Packet Implementer implements one frozen packet at a time without specification edits and MUST NOT commit or push. The Maintainer MUST inspect staged, unstaged and untracked changes, reject weakened tests and verify every claimed result. External private evidence is not given to the Packet Implementer, Independent Reviewer or documentation contributors; the Maintainer runs private regression locally. `pipeline/`, `legacy/`, `docs/audit/`, the historical archive package, its original manifest and `private_evidence/original/` remain immutable evidence. No license choice or private-data upload is authorized by these specifications.

## Governance and completion

Normative words MUST/MUST NOT are release requirements. SHOULD states a recommended choice and requires a recorded reason if changed. A justified methodological change requires a versioned ADR and updates to schema, specification, tests and report language. Breaking scientific/output meaning increments a major contract version; additive capabilities increment minor; compatible corrections increment patch. The initial version is approved for implementation by the user's one-shot request, not scientifically validated by its existence.

The entire roadmap is the delivery objective. A first working CLI, a passed Python-only suite or a partial backend is not completion. A missing runtime, resource, permission or external source metadata may block a specific gate, but Maintainer MUST continue independent work, report NOT_RUN accurately and never substitute fabricated evidence. Public-release ownership/license and study-level scientific provenance can remain explicitly unresolved; core software implementation and executable validation cannot be silently deferred.
