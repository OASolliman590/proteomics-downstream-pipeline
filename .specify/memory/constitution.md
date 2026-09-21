# Downstream Proteomics Constitution

Version: **1.2.0**. Status: **frozen**. Amended and frozen: 2026-09-20.
The prior constitution was ratified on 2026-09-12; this amendment is not yet Maintainer-frozen.

## I. Experimental evidence precedes statistical inference
The system MUST retain assay, source scale, observation/feature identity, original missingness evidence and biological versus technical units. Unknown provenance stays unknown. Quantified protein tables do not establish raw-MS identification confidence, tissue identity or a drug mechanism. See [measurement rules](../../specs/001-downstream-proteomics/contracts/scientific-methods.md#measurement-and-preprocessing).

## II. An analysis plan is an executable scientific contract
The system MUST freeze configuration, inclusion rules, design, contrasts, models, hypotheses, families and local resource hashes before fitting. Limma is the default primary engine. DEqMS/proDA are eligibility-gated adapters, never hit-count-selected competing defaults. Changes create a new plan hash. A null result is not execution failure. See [methods](../../specs/001-downstream-proteomics/contracts/scientific-methods.md).

## III. Biological uncertainty must survive computation
Zero-null, TREAT, equivalence, competitive enrichment, self-contained enrichment and descriptive reversal MUST retain separate result types and families. Confirmatory score testing is disabled unless its protein set, signs and transforms were frozen using data excluding every tested biological unit and subject. In-sample selected scores have no P-value column. Nonsignificance is not equivalence; opposition is not independent corroboration. See [response rules](../../specs/001-downstream-proteomics/contracts/scientific-methods.md#treatment-response-and-scores).

## IV. Reproducibility is verified, not asserted
Stages MUST derive status from executed work, publish atomically and retain hashes, software versions, RNG state and warnings. Production consumes hashed local resources and does not install dependencies or fetch snapshots. Semantic reproduction, independent numerical reference checks and calibration are separate gates. Missing runtime/private evidence is NOT_RUN, not PASS. See [validation](../../specs/001-downstream-proteomics/validation-strategy.md).

## V. The historical baseline remains inspectable
`pipeline/`, `legacy/`, `docs/audit/`, recovered source hashes and external original evidence MUST remain unchanged. The Negr1 auxiliary contrast and the old permutation fractions are regression fixtures, never production constants. Corrections belong in the maintained successor and a separate reconciliation. The protected baseline for this candidate is `3aefdd95c46a1c56dae89e79ddf441b6cc978148`.

## VI. Spec Kit is the development record
FR-001–FR-120, T001–T120 and V001–V120 MUST retain their identity. A merge, split or retirement needs a redirect in traceability and an ADR. Evidence is empty/pending until execution; a proposed path is not evidence. Normative science lives in one shared contract; other documents link to it. This is a file-based kit, not an assumed slash-command installation.

## VII. Small, coherent interfaces and real implementations
The successor MUST use a small Python orchestration layer, a cohesive R statistics package and explicit file contracts. No placeholder backend, fixture-derived production answer or silent method substitution counts as implementation. Scope excludes spectra search, peptide quantification, PTM localization, single-cell, classifier training, web UI, databases and deployed networks. Optional per-study methods still require real software by their release phase.

## VIII. Delegated implementation has an accountable reviewer
The Maintainer owns freeze, contract changes, independent gate execution, private regression and closure. An Independent Reviewer audits the candidate and risky diffs but does not implement engines or close work. A Packet Implementer edits only its frozen allowlist, never specifications, private data, commits or remote state. Private regression and its private logs remain Maintainer-only. No license selection or disclosure is authorized.

## Governance and completion
[PHASES.md](../../specs/001-downstream-proteomics/PHASES.md) defines milestone scope. Phase 1 may be versioned as `v0.1-limma-core`; it MUST NOT be called completed `v1.0-defensible` or imply Phases 2–3 shipped. SYS-03 retains the full twelve-packet objective. Public README and GitHub description MUST match the actually released phase and kit status; release checks verify this.

This pass authorizes specification editing only. Independent audit and explicit Maintainer freeze of an identified tree are still required before any implementation dispatch. The first authorized implementation packet, once frozen, is R01 only. No packet is authorized by the existence of this document. Versioned scientific/output meaning changes require an ADR and synchronized schema, tests and report terms. No waiver permits fabricated PASS, a weakened eligibility rule or an invented license.
