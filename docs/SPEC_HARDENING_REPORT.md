# Specification hardening report

**Frozen kit:** v1.2.0. **Constitution:** 1.2.0. **Scientific contract:** 1.1.0. **Freeze date:** September 20, 2026. **Inspected baseline:** `3aefdd95c46a1c56dae89e79ddf441b6cc978148`.

This deliverable began as a guarded documentation/schema replacement overlay: **80 files (61 replacements, 19 additions)**. It is now applied to an exact-baseline local checkout as uncommitted files, together with the supplied root `AGENTS.md`; it has not been committed, pushed or frozen. **No successor engine was implemented. No implementation packet is authorized by this report.**

## Changes and tightened rules

The constitution retains its eight principles and removes contradictory pre-freeze implementation approval. The roadmap separates Phase 1 (`v0.1-limma-core`), Phase 2 (`v0.2-qualified-methods`) and full Phase 3 (`v1.0-defensible`). All 120 FR/T/V identities and twelve slices remain; R10a covers 091–094 and R10b covers 095–100 as dispatch subdivisions, not new or renumbered requirements. `redirects=[]` is intentional and explained in ADR 0001. Thirteen dispatch units have disjoint write ownership. The only parallel group is R06/R07/R09 after accepted R05 and R10a; R08 waits for all three group members.

Scientific rules specify triggers, required/forbidden behavior, typed results, families and diagnostics. They cover scale/mask/biological-unit integrity, no new primary imputation, exact general-design contrast SEs, limma as default, qualified DEqMS/proDA with no fallback, separate zero-null/TREAT/equivalence/pathway endpoints, primary-versus-secondary multiplicity, offline hashed resources, label-independent mapping, design-compatible enrichment, and covariance-aware descriptive response. Selected-on-tested-data scores have no inferential columns. Exact allocation counting and Monte Carlo sampling are distinct contracts. Historical Negr1/fraction cases are regression fixtures, never production constants.

The Draft 2020-12 configuration now has explicit primary design/engine/hypothesis, per-model execution requiredness, source scale, assay, organism, families, resources and response/score defaults. Primary models are always required; sensitivity models explicitly declare required or optional without inspecting the environment. Scientific eligibility and operational availability are separate frozen-plan fields. Three protocol schemas and the R01 brief close foundation field/interface ambiguities. Three runnable-input examples use only declared tiny synthetic data (8 features, at most 12 observations); no private study matrix or fabricated production analysis output is included. Every V001–V120 has a fixture, independent oracle, exact assertion and negative case. Future evidence remains empty and pending-after-freeze.

R10a retains a real minimal offline HTML target, rather than weakening FR-093 to a Markdown-only stub. Phase 1 requires 54 cases (V001–V050/V091–V094), Phase 2 40 (V051–V090), Phase 3 26 (V095–V120), with report-compatibility rechecks. The 1,000-dataset calibration remains a Phase 3/R11 obligation, not a Phase 1 blocker.

## Contradictions resolved

| Conflict | Resolution |
|---|---|
| Constitution implied approval; parent documents required freeze | All candidate documents require independent audit and explicit Maintainer freeze. |
| Twelve full packets versus useful early delivery | Named milestones preserve full v1.0 requirements without claiming later phases shipped. |
| Shared/wildcard writable files across packets | Exact ownership and fixed lazy service/R-dispatch interfaces replace reopening shared implementation files. |
| R10's full-slice dependency versus Phase 1 report | Separate R10a/R10b dispatch ownership, with explicit full-slice dependency semantics. |
| Paired CAMERA exception versus block-aware dispatch | No paired/subject-blocked CAMERA route in this contract; use the qualified ROAST route. |
| Continuous response bands versus old unlimited reversal | Literal bounded classes and separate crossing-control flags; formal inference is not a class label. |
| Selected-score permutation versus confirmatory independence | Overlap/unknown provenance permits only a descriptive score, physically without P/q columns. |
| Exhaustive and sampled permutation corrections | Exhaustive k/N; separately specified full-space Monte Carlo plus-one sampling. |
| Missing artifacts described as evidence | Empty actual-evidence arrays; future filenames are explicitly pending creation targets. |
| Requested report inside an immutable audit tree | Report written as `docs/SPEC_HARDENING_REPORT.md`; ADR 0001 records this choice. |
| Narrative phase barrier versus machine dispatch dependencies | R06/R07/R09 now explicitly require R05 and R10a; R08 explicitly requires R06/R07/R09. |
| Required/optional engine state depended on installation wording | `model.execution_requirement` is mandatory; scientific eligibility precedes installation inspection and V059 fixes exact FAILED/PARTIAL outcomes. |
| R03 V027/V029 appeared to depend on future R04/R05 services | V027 uses the R01 canonical hash seam with a complete synthetic plan envelope and later V039 regression; R03 owns narrow detection Fisher/BH and R05 later reruns family-separation regression. |

## Verification performed here

| Check | Result and scope |
|---|---|
| JSON parsing and schema validation | PASS: all overlay JSON parses; four Draft 2020-12 schemas; all three examples valid; 16 invalid-configuration cases rejected. |
| Traceability and concrete cases | PASS: 120 unique FR, T and V identities; 12 slices; all cases have four required test sections; no implementation evidence or verified commit asserted. |
| Ownership and dependencies | PASS after correction and independent re-audit: 13 dispatch units, zero allowlist overlaps, R06/R07/R09 require R05 and R10a, and R08 requires the full group barrier. |
| Documentation links | PASS within the overlay; `docs/PIPELINE.md` is the only unchanged baseline-only target, separately confirmed through its GitHub blob `984e71dd77e488c7821447ab95bfc379c2214a07`. This is not a local-byte claim. |
| Tiny input fixtures | PASS: dimensions/IDs/flag spelling and synthetic provenance hashes checked. These hashes concern newly authored fixtures, not historical evidence. |
| Small mathematical oracles | PASS: 2/70 and 6/70 enumeration/counting, shared-control covariance, hand-counted ORA and distinct technical-aggregation estimands. These do not close future implementation V cases. |
| External apply-helper guards | Packaged record: nine tests passed. Actual Windows rerun: eight passed; the symlink-redirection case was NOT_RUN because the account lacks symlink privilege. The helper's own target-symlink check still ran during the successful actual-checkout preflight. |
| Forbidden paths in deliverable | PASS: no overlay file under pipeline/, legacy/, docs/audit/ or maintained production package paths. |
| Original `scripts/check_spec_kit.py` on full checkout | **PASS on the applied local checkout:** 12 slices; 120 unique requirements, tasks and acceptance cases; links and registered source hashes checked. |
| Registered historical SHA-256 and full protected-tree bytes | **PASS on the applied local checkout:** all 9 registered destinations match; `pipeline`, `legacy` and `docs/audit` have no working-tree diff and retain baseline tree IDs. |
| R/statistical backend execution, scientific calibration, private archive regression | **NOT_RUN.** None is evidence supplied by this documentation pass. |
| Independent freeze audit / Maintainer freeze | Initial Astra audit against manifest `e65317b7873d63361399fe8b55f7b90c893f0ae92b5193ff119f6c642a9ee90d` produced three P0 and two P1 findings. Corrections were independently re-audited against manifest `be73544c0999d85adb46c55d5d9e56219b5c90c2e50ec07f80e6be96369c0e85`; no P0 remained and the final wording-only P1 was resolved by clarifying that R handler existence is registration. Maintainer freeze was recorded on 2026-09-20. |
| Repository license / remote description correction | **UNRESOLVED / OUTSTANDING.** No license invented; remote settings not modified. |

The delivery-time checker record and helper-test log are in `verification/` at the delivery-bundle root. The in-repository `checklists/structure-validation.json` now also records the actual-checkout verification. Static success is not software acceptance, R execution, scientific calibration or private regression.

## Applying and finishing the review

The bundle's default application mode is a read-only preflight against a clean checkout at the exact baseline. It verifies registered historical SHA-256 values, validates the payload, exports tracked public source into a temporary directory, overlays the candidate there and invokes the original repository checker. It refuses writes on failure. `--apply` writes only manifest-listed documentation files after that preflight. The full real-checkout preflight and application have now passed; no commit, push, license choice or freeze occurred.

The Maintainer ran the preflight, inspected staged/unstaged/untracked changes, confirmed protected bytes and resolved the independent review findings before recording the frozen contract identity in progress.md. The public GitHub description should be corrected separately to the status text in the release checklist. Licensing remains an explicit owner decision; no software or publication gate may fabricate it.

**Next implementation action, only after successful checkout verification and freeze: R01, from `specs/001-downstream-proteomics/IMPLEMENTATION_BRIEF_R01.md`.** Do not start R02 or any production statistics adapter in this pass.

## File inventory

“Replace” means replacement bytes supplied for a path observed in the baseline; it does not imply those remote bytes have already been changed. “Add” means a new documentation/contract/fixture path. The external `MANIFEST.json` provides each delivered file's SHA-256 and byte count.

| Action | Repository-relative path |
|---|---|
| Replace | `.specify/memory/constitution.md` |
| Replace | `README.md` |
| Replace | `docs/RELEASE_CHECKLIST.md` |
| Add | `docs/SPEC_HARDENING_REPORT.md` |
| Add | `docs/adr/0001-phase-milestones.md` |
| Add | `docs/adr/0002-primary-engine-limma.md` |
| Add | `docs/adr/0003-response-descriptive-default.md` |
| Add | `docs/adr/0004-enrichment-dispatch.md` |
| Add | `docs/adr/0005-no-silent-method-fallback.md` |
| Add | `specs/001-downstream-proteomics/IMPLEMENTATION_BRIEF_R01.md` |
| Add | `specs/001-downstream-proteomics/PHASES.md` |
| Replace | `specs/001-downstream-proteomics/START_HERE.md` |
| Replace | `specs/001-downstream-proteomics/checklists/spec-quality.md` |
| Replace | `specs/001-downstream-proteomics/checklists/structure-validation.json` |
| Replace | `specs/001-downstream-proteomics/contracts/analysis.schema.json` |
| Replace | `specs/001-downstream-proteomics/contracts/cli-and-artifacts.md` |
| Replace | `specs/001-downstream-proteomics/contracts/example-effect-threshold.json` |
| Replace | `specs/001-downstream-proteomics/contracts/example-independent.json` |
| Replace | `specs/001-downstream-proteomics/contracts/example-paired.json` |
| Add | `specs/001-downstream-proteomics/contracts/fixtures/README.md` |
| Add | `specs/001-downstream-proteomics/contracts/fixtures/features.tsv` |
| Add | `specs/001-downstream-proteomics/contracts/fixtures/independent-abundance.tsv` |
| Add | `specs/001-downstream-proteomics/contracts/fixtures/independent-observations.tsv` |
| Add | `specs/001-downstream-proteomics/contracts/fixtures/paired-abundance.tsv` |
| Add | `specs/001-downstream-proteomics/contracts/fixtures/paired-observations.tsv` |
| Add | `specs/001-downstream-proteomics/contracts/fixtures/provenance.json` |
| Add | `specs/001-downstream-proteomics/contracts/run-status.schema.json` |
| Replace | `specs/001-downstream-proteomics/contracts/scientific-methods.md` |
| Replace | `specs/001-downstream-proteomics/contracts/semantic-validation.md` |
| Add | `specs/001-downstream-proteomics/contracts/stage-request.schema.json` |
| Add | `specs/001-downstream-proteomics/contracts/stage-result.schema.json` |
| Replace | `specs/001-downstream-proteomics/data-model.md` |
| Replace | `specs/001-downstream-proteomics/packet-index.md` |
| Add | `specs/001-downstream-proteomics/packet-ownership.json` |
| Replace | `specs/001-downstream-proteomics/plan.md` |
| Replace | `specs/001-downstream-proteomics/progress.md` |
| Replace | `specs/001-downstream-proteomics/quickstart.md` |
| Replace | `specs/001-downstream-proteomics/research.md` |
| Replace | `specs/001-downstream-proteomics/roadmap.json` |
| Replace | `specs/001-downstream-proteomics/roadmap.md` |
| Replace | `specs/001-downstream-proteomics/spec.md` |
| Replace | `specs/001-downstream-proteomics/tasks.md` |
| Replace | `specs/001-downstream-proteomics/traceability.json` |
| Replace | `specs/001-downstream-proteomics/validation-strategy.md` |
| Replace | `specs/002-foundation/plan.md` |
| Replace | `specs/002-foundation/spec.md` |
| Replace | `specs/002-foundation/tasks.md` |
| Replace | `specs/003-intake/plan.md` |
| Replace | `specs/003-intake/spec.md` |
| Replace | `specs/003-intake/tasks.md` |
| Replace | `specs/004-preprocessing-qc/plan.md` |
| Replace | `specs/004-preprocessing-qc/spec.md` |
| Replace | `specs/004-preprocessing-qc/tasks.md` |
| Replace | `specs/005-design-contrasts/plan.md` |
| Replace | `specs/005-design-contrasts/spec.md` |
| Replace | `specs/005-design-contrasts/tasks.md` |
| Replace | `specs/006-limma-inference/plan.md` |
| Replace | `specs/006-limma-inference/spec.md` |
| Replace | `specs/006-limma-inference/tasks.md` |
| Replace | `specs/007-assay-engines/plan.md` |
| Replace | `specs/007-assay-engines/spec.md` |
| Replace | `specs/007-assay-engines/tasks.md` |
| Replace | `specs/008-resources-mapping/plan.md` |
| Replace | `specs/008-resources-mapping/spec.md` |
| Replace | `specs/008-resources-mapping/tasks.md` |
| Replace | `specs/009-enrichment/plan.md` |
| Replace | `specs/009-enrichment/spec.md` |
| Replace | `specs/009-enrichment/tasks.md` |
| Replace | `specs/010-treatment-response/plan.md` |
| Replace | `specs/010-treatment-response/spec.md` |
| Replace | `specs/010-treatment-response/tasks.md` |
| Replace | `specs/011-reporting/plan.md` |
| Replace | `specs/011-reporting/spec.md` |
| Replace | `specs/011-reporting/tasks.md` |
| Replace | `specs/012-validation/plan.md` |
| Replace | `specs/012-validation/spec.md` |
| Replace | `specs/012-validation/tasks.md` |
| Replace | `specs/013-release/plan.md` |
| Replace | `specs/013-release/spec.md` |
| Replace | `specs/013-release/tasks.md` |
