# Feature specification — defensible downstream protein abundance

**Version:** 1.2.0-frozen. **Status:** maintained implementation NOT_STARTED; R01 requires explicit GO before dispatch.

## Product and user journeys

The target consumes quantified protein-level matrices and metadata, then performs intake/identity validation → QC → frozen design → eligible inference → design-valid enrichment → treatment-response description → honest offline reporting. It does not perform raw instrument conversion/search, peptide quantification/protein inference from spectra, PTM localization, single-cell analysis, classifier training, web UI, database or network deployment, or causal drug-mechanism inference.

US1 validates scale, assay, missingness, biological units and a frozen analysis plan. US2 obtains exact, qualified abundance inference with separate zero-null/TREAT/omnibus endpoints and multiplicity. US3 obtains species-aware mapping and enrichment with the correct null. US4 describes response without circular rescue claims and performs only eligible independent-score/equivalence/conjunction inference. US5 reviews QC, all contrasts and failures in one offline report, including zero discoveries. US6 reproduces semantic tables under qualified environments/resources and records validation and legacy limits.

## Parent requirements

- **SYS-01:** One maintained CLI MUST expose coherent validation/planning/run/verification services and the full downstream workflow by its specified phases; phase scope is visible in outputs.
- **SYS-02:** Every stage MUST preserve explicit scale, assay, observed/imputed evidence, unit hierarchy, feature identity, model/contrast meaning and provenance.
- **SYS-03:** All twelve specified slices MUST be implemented and accepted for full v1.0-defensible. R01–R05 plus R10a may be accepted only as the v0.1-limma-core milestone; it does not complete R10 or imply Phases 2–3 shipped.
- **SYS-04:** Statistical selection MUST use prespecified eligibility and primary-method declarations, not hit counts; inapplicable/failed methods retain typed reasons and never silently fall back.
- **SYS-05:** The recovered baseline and audit MUST remain byte-preserved; private regression is Maintainer-only and unexecuted reproduction remains NOT_RUN.
- **SYS-06:** Routine analysis MUST require a declared configuration and one run command; deterministic design work precedes plan freeze and all fitting follows it.
- **SYS-07:** Release claims MUST be supported by actual tests, model/resource/environment qualification and reviewed evidence for the claimed phase.
- **SYS-08:** Automated output MUST NOT assert tissue, identity, mechanism, raw-MS confidence, study-wide FDR or restoration beyond the actual evidence/hypothesis/family.

Detailed requirements retain FR-001–FR-120 with one-to-one T001–T120 and V001–V120 identities in [traceability.json](traceability.json) and the linked [slice specs](packet-index.md). New dispatch subdivisions do not create or retire requirements. The [scientific contract](contracts/scientific-methods.md) is the vocabulary/eligibility authority; other documents cannot weaken it.

## Success and phase acceptance

SC-01: Actual end-to-end examples produce typed traceable outputs through the claimed phase. SC-02: All 120 requirements retain reviewed evidence or explicit pending/external status; no fabricated evidence path. SC-03: Golden numerical tests meet [fixed tolerances](validation-strategy.md). SC-04: R11 calibration meets its prespecified error/coverage limits before a full validated-method claim. SC-05: Pinned offline semantic reproduction is verified by R11. SC-06: Maintainer-run legacy reconciliation explicitly handles historical exceptions and intentional changes. SC-07: The R11 20,000×100/eight-contrast/2,000-set benchmark is measured against ≤8 GiB and ≤30 minutes on documented four-core hardware. SC-08: A versioned successor is exported only after the relevant gates and release authorization, without overwriting originals or an unauthorized remote push.

## Edge cases

Primary acceptance covers log2 zero/second log, technical pseudoreplication, confounded batch, featurewise missing-data estimability, weighted/blocked exact SE, no discoveries and failed stages. Later phases cover count/dropout evidence, ambiguous mappings, zero-overlap ORA, finite pathway matrices, independent versus gene-set nulls, unstable RI, overshoot, equivalence versus nonsignificance, selection overlap and exact/Monte Carlo enumeration. Full cases and negative oracles are in every slice, not generic boilerplate.

## Current reality

The recovered study remains the only existing analysis implementation. R was not rerun during recovery; environment.yml is not a lockfile; gene-set resources and a repository license remain unresolved. [PHASES.md](PHASES.md) and [progress.md](progress.md) distinguish this baseline from future milestones. No implementation packet is authorized before independent audit and Maintainer freeze.
