# Corrected packet index — draft pending audit and freeze

**Kit version:** 0.2-draft. **Revision date:** 2026-09-20. **State:** corrected draft; independent audit and maintainer freeze are pending; no implementation dispatch is authorized.

This index is the orchestration source of truth for the existing twelve slices. The detailed scientific requirements remain the 120 `FR-*` / `T*` / `V*` entries in the numbered slice specifications and `traceability.json`; this index changes ownership, ordering and review boundaries without adding or weakening product scope.

## Decisions proposed for freeze

1. The scientific scope and acceptance strength in FR-001–FR-120 and V001–V120 are unchanged.
2. `pipeline/`, `legacy/` and `docs/audit/` are immutable historical evidence. Maintained work is confined to the packet allowlists below.
3. The maintainer owns kit changes, private-data regression, verification, commits, progress/traceability updates, delivery and closure. An independent reviewer audits the specification freeze and risky diffs without implementing or closing work. A packet implementer works on one frozen packet without editing the specification.
4. Private study data and the historical archive package never enter a worker task. Only synthetic or explicitly redistributable fixtures may be committed.
5. A packet is complete only when every acceptance ID in its range has inspected evidence on a reviewed commit. `NOT_RUN`, `SKIPPED`, `INAPPLICABLE`, missing runtime and missing dependency are never `PASS`.
6. R01 begins from the recovered baseline and is not credited as complete until all V001–V010 evidence is reviewed on a commit.
7. Sequential packets may reopen named shared files only after their dependencies are complete. Parallel packets may not share any writable path.
8. Documentation work in R12 follows the same audit, freeze, verification and closure requirements as every other packet.

## Global forbidden paths and actions

Every packet forbids `pipeline/**`, `legacy/**`, `docs/audit/**`, `private_evidence/**`, the original archive package and manifest, credentials, global configuration, unrelated files, and every maintained path not in that packet's allowlist. Packet implementers must not edit `.specify/**`, `specs/**` or the project-level progress and traceability records. They must not commit, push, publish, upload data, choose a license or change branches.

## Packet table

| ID | Title | Goal | Dependencies | Acceptance | Risk |
|---|---|---|---|---|---|
| R01 | Runtime, package and contract foundation | Establish executable package, schemas, state/runtime, R bridge, baseline registry and real test harness. | none | V001–V010 | high |
| R02 | Canonical protein input and biological identity | Canonicalize protein matrices, scale/missingness, identities, technical units, vendor and legacy imports. | R01 | V011–V020 | high |
| R03 | Preprocessing, missingness and quality control | Preserve declared abundance semantics while implementing qualified normalization, filtering, missingness, QC and sensitivities. | R02 | V021–V030 | high |
| R04 | Design validation, exact contrasts and blocking | Freeze safe designs, contrasts, estimability, blocking and method applicability before fitting. | R03 | V031–V040 | high |
| R05 | Core limma inference and multiplicity | Produce validated limma estimates, exact contrasts, threshold tests, uncertainty, families and complete typed exports. | R04 | V041–V050 | high |
| R06 | Assay-qualified DEqMS and proDA backends | Add evidence-qualified alternative engines without silent substitution or borrowed uncertainty. | R05 | V051–V060 | high |
| R07 | Versioned annotation, protein groups and gene sets | Create offline versioned resources, defensible mapping/gene aggregation and measured pathway universes. | R05 | V061–V070 | high |
| R08 | Design-compatible pathways and enrichment | Dispatch CAMERA/ROAST/fgsea/ORA by their valid nulls and preserve complete pathway evidence. | R06, R07 | V071–V080 | high |
| R09 | Treatment response, equivalence and independent scores | Separate descriptive movement, bounded response, equivalence, conjunction and independent-score inference. | R05 | V081–V090 | high |
| R10 | Unified offline report and streamlined workflow | Integrate the verified DAG into one command and render honest offline reports and figure-source data. | R03, R06, R08, R09 | V091–V100 | high |
| R11 | Reproducibility, calibration and continuous validation | Pin environments, verify caching/offline behavior, scientific oracles, calibration, CI, performance and evidence ledger. | R10 | V101–V110 | high |
| R12 | Archived-study regression, documentation and versioned successor | Run private local regression, reconcile differences, complete documentation/traceability and build a versioned SSD successor. | R11 | V111–V120 | high |

## Packet contracts

### R01 — Runtime, package and contract foundation

- **Goal:** Establish executable package, schemas, state/runtime, R bridge, baseline registry and real test harness.
- **Files allowed:** `pyproject.toml`; `src/proteomics_pipeline/{__init__.py,__main__.py,cli.py,config.py,doctor.py,errors.py,paths.py,provenance.py,runtime.py,py.typed}`; `src/proteomics_pipeline/schemas/**`; `schemas/{analysis.schema.json,stage-request.schema.json,stage-result.schema.json}`; `r/proteomicsCore/{DESCRIPTION,NAMESPACE,.Rbuildignore,R/**,tests/testthat.R,tests/testthat/test-io-roundtrip.R}`; `scripts/maintained/{run_stage.R,test_r.R,check_baseline.py,write_solved_versions.py,bootstrap.ps1,bootstrap.sh}`; `configs/examples/{synthetic-independent.yml,synthetic-paired.yml,synthetic-effect-threshold.yml}`; `tests/contract/{baseline_test.py,bootstrap_test.py,cli_install_test.py,config_schema_test.py,harness_test.py,r_package_test.py}`; `tests/unit/{cli_validate_test.py,doctor_test.py,provenance_test.py,r_bridge_test.py,runtime_state_test.py}`; `tests/fixtures/echo_argv.py`; `docs/development.md`; `docs/validation/{baseline-hashes.json,002-foundation/**}`.
- **Files forbidden:** global forbidden paths plus `.gitignore`, all later-packet source/tests/evidence, and any R analysis implementation beyond the minimal non-fabricating I/O/runtime skeleton.
- **Dependencies:** none; recovered baseline and approved parent contracts are read-only inputs.
- **Acceptance lines:** all V001–V010 criteria pass exactly as written; historical hashes are verified and mutation-negative-tested; both Python and real R gates report missing prerequisites honestly.
- **Risk:** high — defines shared schemas, runtime states, package layout and evidence baseline used by every later packet.

### R02 — Canonical protein input and biological identity

- **Goal:** Canonicalize protein matrices, scale/missingness, identities, technical units, vendor and legacy imports.
- **Files allowed:** `src/proteomics_pipeline/intake/**`; `r/proteomicsCore/R/intake.R`; `configs/mappings/**`; `tests/contract/test_intake.py`; `tests/scientific/test_technical_units.py`; `tests/fixtures/synthetic/intake/**`; `docs/validation/003-intake/**`.
- **Files forbidden:** global forbidden paths plus R03–R12 implementation files and R01 shared schemas/runtime except through their frozen interfaces.
- **Dependencies:** R01.
- **Acceptance lines:** all V011–V020 criteria pass exactly as written; aligned wide/long roundtrips preserve values, masks and identities; scale/zero/missing encodings and biological-versus-technical units are explicit; unsupported grain/vendor inputs fail without inference.
- **Risk:** high — canonical identity, scale and observation grain determine every downstream estimand.

### R03 — Preprocessing, missingness and quality control

- **Goal:** Preserve declared abundance semantics while implementing qualified normalization, filtering, missingness, QC and sensitivities.
- **Files allowed:** `r/proteomicsCore/R/{preprocess.R,normalization.R,tmt.R,filtering.R,missingness.R,qc.R,sensitivity.R,detection.R}`; `src/proteomics_pipeline/planning.py`; `src/proteomics_pipeline/reporting/qc_data.py`; `tests/integration/test_qc.py`; packet-specific testthat files named `test-preprocessing*.R` or `test-qc*.R`; `docs/validation/004-preprocessing-qc/**`.
- **Files forbidden:** global forbidden paths plus R04–R12 implementation files and edits to R02 canonical artifacts.
- **Dependencies:** R02.
- **Acceptance lines:** all V021–V030 criteria pass exactly as written; processed log2 data are not transformed again; primary data remain distinct from imputation sensitivities; TMT bridge/no-bridge eligibility, contrast-aware coverage, outlier policy and QC source tables are testable and honest.
- **Risk:** high — transformations, filters and missingness policy can silently change the feature universe and all later results.

### R04 — Design validation, exact contrasts and blocking

- **Goal:** Freeze safe designs, contrasts, estimability, blocking and method applicability before fitting.
- **Files allowed:** `src/proteomics_pipeline/planning.py`; `r/proteomicsCore/R/{design.R,estimability.R,contrasts.R,capabilities.R}`; `tests/scientific/test_design_contract.py`; `r/proteomicsCore/tests/testthat/test-design.R`; `docs/validation/005-design-contrasts/**`.
- **Files forbidden:** global forbidden paths plus R05–R12 implementation files and arbitrary R evaluation from user configuration.
- **Dependencies:** R03.
- **Acceptance lines:** all V031–V040 criteria pass exactly as written; confounding/rank failures block before fit; pairing/repeated-unit and interaction semantics preserve the experimental unit; sparse/nonorthogonal contrast covariance matches an independent oracle; plan changes invalidate hashes.
- **Risk:** high — shared design, contrast and capability contracts govern every inferential engine.

### R05 — Core limma inference and multiplicity

- **Goal:** Produce validated limma estimates, exact contrasts, threshold tests, uncertainty, families and complete typed exports.
- **Files allowed:** `r/proteomicsCore/R/{model_limma.R,diagnostics.R,inference.R,multiplicity.R,export.R,influence.R}`; `schemas/differential-result.schema.json`; `tests/scientific/test_limma.py`; `tests/scientific/fixtures/limma/**`; `r/proteomicsCore/tests/testthat/test-limma.R`; `docs/validation/006-limma-inference/**`.
- **Files forbidden:** global forbidden paths plus alternative engines, resources, pathways, response and report implementation files.
- **Dependencies:** R04.
- **Acceptance lines:** all V041–V050 criteria pass exactly as written; official direct limma calls and separate sparse/nonorthogonal oracles agree within frozen tolerances; TREAT, zero-null, intervals, F tests and multiplicity families remain distinct; every planned row is tested or carries an explicit reason.
- **Risk:** high — primary statistical estimates, uncertainty and multiplicity semantics feed all downstream claims.

### R06 — Assay-qualified DEqMS and proDA backends

- **Goal:** Add evidence-qualified alternative engines without silent substitution or borrowed uncertainty.
- **Files allowed:** `schemas/count-evidence.schema.json`; `src/proteomics_pipeline/intake/counts.py`; `src/proteomics_pipeline/planning.py`; `src/proteomics_pipeline/doctor.py`; `src/proteomics_pipeline/reporting/model_data.py`; `r/proteomicsCore/R/{model_deqms.R,model_proda.R,capabilities.R,inference.R,sensitivity.R}`; `tests/scientific/test_assay_engines.py`; `r/proteomicsCore/tests/testthat/test-assay-engines.R`; `docs/validation/007-assay-engines/**`.
- **Files forbidden:** global forbidden paths plus every R07 and R09 allowed path, all pathway/response algorithms, and edits to R05's primary limma adapter/multiplicity implementation.
- **Dependencies:** R05.
- **Acceptance lines:** all V051–V060 criteria pass exactly as written; count provenance and LFQ dropout eligibility are enforced; direct official DEqMS/proDA references agree; unsupported block/weight/threshold cases fail explicitly; alternative-engine failures never alter the primary engine silently.
- **Risk:** high — shared planning/capability/inference files and method-specific statistical assumptions can break later dispatch.

### R07 — Versioned annotation, protein groups and gene sets

- **Goal:** Create offline versioned resources, defensible mapping/gene aggregation and measured pathway universes.
- **Files allowed:** `src/proteomics_pipeline/resources.py`; `schemas/resource.schema.json`; `configs/resources/**`; `r/proteomicsCore/R/{mapping.R,gene_matrix.R,genesets.R}`; `tests/contract/test_resources.py`; `r/proteomicsCore/tests/testthat/test-mapping.R`; `docs/validation/008-resources-mapping/**`; small redistributable `resources/**` with provenance.
- **Files forbidden:** global forbidden paths plus every R06 and R09 allowed path, auto-fetch behavior during analysis, and private/unlicensed resource content.
- **Dependencies:** R05.
- **Acceptance lines:** all V061–V070 criteria pass exactly as written; every snapshot is versioned/hashed/offline; mappings and orthology retain ambiguity/evidence; representative selection is label-independent; measured mapped universes and overlap loss are explicit.
- **Risk:** high — resource identity, mapping and universe definitions directly change pathway hypotheses and multiplicity.

### R08 — Design-compatible pathways and enrichment

- **Goal:** Dispatch CAMERA/ROAST/fgsea/ORA by their valid nulls and preserve complete pathway evidence.
- **Files allowed:** `r/proteomicsCore/R/{pathways.R,pathway_camera.R,pathway_roast.R,pathway_fgsea.R,pathway_ora.R,pathway_summary.R,multiplicity.R}`; `tests/scientific/test_pathways.py`; `tests/scientific/fixtures/pathways/**`; `r/proteomicsCore/tests/testthat/test-pathways.R`; `docs/validation/009-enrichment/**`.
- **Files forbidden:** global forbidden paths plus resource-building/mapping files, assay-engine fit files and treatment-response files.
- **Dependencies:** R06 and R07.
- **Acceptance lines:** all V071–V080 criteria pass exactly as written; method/design/null compatibility is enforced; CAMERA, ROAST, fgsea and ORA match independent references; directional/mixed and central/native q fields stay distinct; all tested sets and honest empty/failure states are exported.
- **Risk:** high — null choice, correlation/block handling and family adjustment determine the validity of pathway claims.

### R09 — Treatment response, equivalence and independent scores

- **Goal:** Separate descriptive movement, bounded response, equivalence, conjunction and independent-score inference.
- **Files allowed:** `r/proteomicsCore/R/{response.R,response_uncertainty.R,equivalence.R,scores.R,randomization.R}`; `schemas/independent-score.schema.json`; `src/proteomics_pipeline/reporting/response_data.py`; `tests/scientific/test_response.py`; `tests/scientific/fixtures/response/**`; `r/proteomicsCore/tests/testthat/test-response.R`; `docs/validation/010-treatment-response/**`.
- **Files forbidden:** global forbidden paths plus every R06 and R07 allowed path, pathway files and primary model-fitting files.
- **Dependencies:** R05.
- **Acceptance lines:** all V081–V090 criteria pass exactly as written; covariance-aware disease/treatment/residual axes and RI boundary cases match analytic examples; nonsignificance never becomes equivalence; selected-on-current-data scores remain descriptive; exact/randomized tests count allocations and ties correctly.
- **Risk:** high — circular selection, ratio instability and invalid rescue/equivalence language can create false biological claims.

### R10 — Unified offline report and streamlined workflow

- **Goal:** Integrate the verified DAG into one command and render honest offline reports and figure-source data.
- **Files allowed:** `src/proteomics_pipeline/{cli.py,runtime.py}`; `src/proteomics_pipeline/reporting/{data.py,render.py,qc_data.py,model_data.py,pathway_data.py,response_data.py,methods.py,templates/**}`; `schemas/report.schema.json`; `r/proteomicsCore/R/plots.R`; `tests/integration/test_reports.py`; `docs/user-guide/**`; `configs/examples/README.md`; `docs/validation/011-reporting/**`.
- **Files forbidden:** global forbidden paths plus changes to scientific estimators, mapping, pathway tests or response tests; R01 example data/config files are read-only inputs.
- **Dependencies:** R03, R06, R08 and R09.
- **Acceptance lines:** all V091–V100 criteria pass exactly as written; one command follows the valid DAG and stops required failures; reports are offline, accessible and status-derived; every table/figure links exact source data; null, empty, partial and failed runs remain renderable without false PASS language.
- **Risk:** high — reopens shared CLI/runtime/report adapters and controls the language users see.

### R11 — Reproducibility, calibration and continuous validation

- **Goal:** Pin environments, verify caching/offline behavior, scientific oracles, calibration, CI, performance and evidence ledger.
- **Files allowed:** `renv.lock`; the selected Python lockfile; `containers/**`; `src/proteomics_pipeline/{provenance.py,runtime.py}`; `tests/integration/{test_offline.py,test_resume.py,test_report_artifacts.py}`; `tests/scientific/reference/**`; `tests/scientific/calibration_pathways.R`; `scripts/maintained/{run_calibration.py,benchmark.py}`; packet-specific CI helpers; `.github/workflows/maintained.yml`; `r/proteomicsCore/R/simulate.R`; `schemas/validation-evidence.schema.json`; `docs/validation/{reference-matrix.json,acceptance.json,benchmarks/**,012-validation/**}`.
- **Files forbidden:** global forbidden paths plus product estimators/report logic except the named runtime/provenance cache seams; no reduced fixtures or relaxed thresholds may replace failed targets.
- **Dependencies:** R10.
- **Acceptance lines:** all V101–V110 criteria pass exactly as written; clean compatible Python/R/Bioconductor restoration and offline rerun are demonstrated; cache invalidation is semantic; independent references and prespecified calibration retain failures; Windows/Linux, performance and rendered-report evidence distinguish PASS from NOT_RUN.
- **Risk:** high — environment locks, shared cache semantics and validation criteria determine whether any release claim is defensible.

### R12 — Archived-study regression, documentation and versioned successor

- **Goal:** Run private local regression, reconcile differences, complete documentation/traceability and build a versioned successor.
- **Files allowed:** `scripts/maintained/{legacy_regression.py,build_release.py}`; `tests/regression/**`; local-only config template with no study data; `docs/validation/{legacy-comparison.md,delivery.json,final-acceptance.md,013-release/**}`; `docs/user-guide/**`; `docs/methods/**`; `README.md`; `CHANGELOG.md`; `.gitignore`; maintainer-owned `specs/001-downstream-proteomics/{traceability.json,progress.md}`; a new versioned delivery directory under `<release-root>/releases/` only after verification.
- **Files forbidden:** global forbidden paths; the original archive package, manifest and `<private-evidence-root>/original/**` are read-only even for the maintainer; packet implementers, independent reviewers and documentation contributors may not receive private inputs or private-derived contents.
- **Dependencies:** R11.
- **Acceptance lines:** all V111–V120 criteria pass exactly as written, with V111/V114 executed locally by Maintainer; legacy non-reproduction and intentional scientific differences remain explicit; docs match implemented/verified behavior; final source whitelist excludes private content; delivery is versioned and the original hashes remain unchanged.
- **Risk:** high — private evidence, historical reconciliation, documentation claims, release hygiene and external-storage writes require serial maintainer control.

## Dependency graph

```text
R01 → R02 → R03 → R04 → R05
                         ├─→ R06 ─┐
                         ├─→ R07 ─┼─→ R08 ─┐
                         └─→ R09 ───────────┼─→ R10 → R11 → R12
R03 ────────────────────────────────────────┘
R06 ────────────────────────────────────────┘
```

Direct dependency edges are: `R01→R02`, `R02→R03`, `R03→R04`, `R04→R05`, `R05→R06`, `R05→R07`, `R05→R09`, `R06→R08`, `R07→R08`, `R03→R10`, `R06→R10`, `R08→R10`, `R09→R10`, `R10→R11`, `R11→R12`.

## Parallel groups

Only one group is legal in this draft:

- **PG-01:** R06, R07 and R09, after R05 is verified. Their allowlists and named Python/R test files are disjoint. They have no dependency on one another. R08 waits for both R06 and R07; R10 waits for R06, R08 and R09.

All other packets are serial. Tests under a shared directory are not shared ownership: only the exact named test files in each allowlist are writable.

## Risky pile

All twelve packets are high-risk at slice granularity. The specifically serial/high-coupling seams are:

- R01 shared package/schema/runtime/evidence foundation.
- R02 canonical identifiers, scale, missing encodings and biological-unit grain.
- R03–R06 sequential handoffs through `planning.py`, `capabilities.R`, `inference.R` and `sensitivity.R`.
- R05/R08 sequential handoff through `multiplicity.R`.
- R10 shared `cli.py`, `runtime.py` and report adapters that can misstate upstream results.
- R11 shared cache/provenance semantics, environment locks, calibration and CI truthfulness.
- R12 private regression, documentation claims, release whitelist and external delivery.

No authentication or database migration packet exists in product scope. Delegate authentication is not a scientific acceptance substitute and is not part of the frozen product.

## Questions blocking freeze

None identified in the current revision. Independent specification audit remains a mandatory gate before freeze; any unresolved disagreement on a product requirement must be escalated to the project owner.
