# Feature Specification: Streamlined downstream proteomics

**Feature Branch**: `spec/001-downstream-proteomics` (logical epic; the app may allocate the actual worktree branch).
**Created**: 2026-09-12. **Status**: Corrected draft awaiting Independent Reviewer audit and Maintainer freeze; functionality not yet verified in this checkout.
**Input**: Mature the recovered dapagliflozin code into a scientifically defensible, streamlined, complete downstream protein-abundance pipeline. Develop through detailed Spec Kit specifications and frozen packet contracts. The maintainer owns the kit, review and closure; an independent reviewer audits specification freeze and risky diffs; a packet implementer implements one frozen packet at a time.

## User Scenarios & Testing

### US1 — Start a scientifically interpretable run (Priority P1)

An analyst supplies a protein abundance matrix, sample metadata and a configuration. They receive actionable validation before any model runs, a frozen plan and an immutable run directory. The pipeline handles general sample IDs, group sizes and species rather than the old fixed twenty-rat layout.

**Independent test:** execute validate/plan on synthetic independent, paired, multiplex and invalid fixtures. No R model is required to demonstrate intake correctness.

**Acceptance:** Given repeated technical injections, when a biological contrast is planned, then injections are explicitly aggregated or rejected as unmodeled replication. Given processed log2 values, when the run begins, then no second log transform is applied. Given a sample-group/batch confounding, then the affected contrast is blocked before fitting with a rank/alias explanation.

### US2 — Obtain valid abundance inference (Priority P1)

An analyst estimates planned disease, treatment, continuous-covariate, interaction or paired effects with an eligible model. The output separates zero-null evidence, practical effect thresholds, nonestimable features and exploratory sensitivity results.

**Independent test:** supported engines match pinned official reference calls on complete, sparse and nonorthogonal synthetic designs; row/sample permutations do not change aligned estimates.

**Acceptance:** Given a threshold-null test, results identify the declared effect threshold and never substitute post-hoc filtering of a zero-null P value. Given partially missing nonorthogonal designs, the model uses exact contrast parametrization or featurewise covariance rather than approximate contrasts.fit standard errors. Given a method that cannot model the requested correlation structure, validation rejects it rather than ignoring arguments.

### US3 — Investigate pathways with the right null (Priority P1)

The analyst selects versioned gene sets and explicit protein-to-gene mapping, sees coverage/universe loss and gets pathway results compatible with the design. Reported FDR is tied to a declared test family, and dependent pathways are not presented as independent discoveries.

**Independent test:** independent linear designs exercise CAMERA; blocked designs exercise ROAST; a small hand-counted universe exercises ORA; fixed rank/seed exercises fgsea. Missing resources and all-empty overlaps remain visible.

**Acceptance:** A blocked model cannot pass block/correlation to CAMERA and claim those were used. ROAST directional and mixed hypotheses have distinct rows/families. Preranked fgsea is labeled exploratory with a gene-set null. All significant and nonsignificant tested sets are exported, with full diagnostics.

### US4 — Describe treatment response without circular rescue claims (Priority P1)

An analyst can inspect disease, treatment, residual and interaction effects with uncertainty, directional opposition and bounded descriptive restoration classes. They can test a genuinely independent fixed score or prespecified practical equivalence when eligibility is satisfied.

**Independent test:** null shared-control simulation, exact enumeration, overshoot examples and analytic TOST cases expose the old false interpretations.

**Acceptance:** Selecting proteins on the current disease contrast disables confirmatory fixed-score testing. RI=3 is overshoot, never full restoration. Nonsignificance of treated-versus-control cannot produce an equivalence label. All allocations in a small exact test are counted once with tie tolerance.

### US5 — Review QC and limitations in one report (Priority P1)

An analyst runs one command and reads an HTML report with QC, eligible/failed methods, complete contrast coverage, results, pathway caveats and methods/provenance. Null results, empty candidate sets and optional unavailable metadata do not break report construction.

**Independent test:** complete, null, partial-failure and no-pathway synthetic runs render offline with validated numerical contents, original-scale/missingness notes and traceable figure-source data.

### US6 — Reproduce and extend the pipeline (Priority P1)

A computational scientist installs a pinned environment, executes documented examples, reproduces semantic outputs, audits scientific calibration and safely resumes compatible stages. New datasets use schema-driven configuration. Private study data stay local.

**Independent test:** clean environment, offline cached rerun, semantic table hashes, null/mixture simulations and local archive regression. A code/config/resource change invalidates the relevant stage cache.

## Requirements

The normative detailed requirements are in numbered slice specifications 002–013, the live constitution and shared contracts. Their globally unique FR-001…FR-120 identifiers map one-to-one to implementable tasks and then acceptance evidence. The parent requirements below bind all slices:

- **SYS-01:** One maintained `proteomics` CLI MUST drive intake, validation, planning, QC, differential abundance, pathways, treatment-response summaries and offline HTML reporting; subcommands expose the same stages without divergent logic.
- **SYS-02:** All analyses MUST preserve explicit source scale, imputation history, sample hierarchy, feature identity, filter lineage, model/contrast identity and provenance. Unknowns remain visible.
- **SYS-03:** All twelve roadmap slices MUST be implemented; method adapters cannot remain placeholders or substitute silently for one another.
- **SYS-04:** Statistical defaults MUST be justified by applicability and primary references, not significance counts. Missing scientific prerequisites MUST block only dependent claims/stages with a typed status.
- **SYS-05:** The historical baseline MUST be byte-preserved and regression-tested where executable. New decisions, output changes and remaining uncertainties MUST be recorded separately.
- **SYS-06:** Routine successful operation MUST need only a configuration and one run command. Complexity belongs in the frozen plan and diagnostics, not in manually invoking six R scripts.
- **SYS-07:** The implementation MUST deliver real tests, pinned runtime/resource manifests, executable examples, coherent error/status handling and an independently verified acceptance report.
- **SYS-08:** Automated output MUST not assert species/tissue/rat identity, causal effects, global FDR, raw-MS quality or restoration beyond the data and tested hypotheses.

## Key Entities

Study, Specimen, Observation, Feature, Assay, MatrixArtifact, AnalysisPlan, Design, Contrast, ModelFit, DifferentialResult, ResourceSnapshot, GeneMapping, GeneSetResult, ResponseSummary, IndependentScore, Run, StageResult and ValidationEvidence. Exact fields, keys and relationships are specified in data-model.md and contracts.

## Edge Cases

Zero/negative linear values with explicit missing encodings; zero on log2 scale; duplicate feature IDs; ambiguous protein groups; duplicate sample rows; technical replicate imbalance; absent factor levels; all missing values in a relevant group; a single biological replicate; aliased batch and treatment; unequal group sizes; repeated observations; count covariates with wrong grain; reimputed data offered to a dropout model; missing or negative assay weights; no mapped genes; all pathway tests nonfinite; no significant proteins; zero disease denominator; extreme RI overshoot; nonexchangeable permutation blocks; truncated R output; stage crash; interrupted writes; changed plan on resume; cached database mismatch; unavailable optional engine; font/render failure; private data appearing in a release manifest.

## Success Criteria

- **SC-01:** The documented synthetic end-to-end example completes from one run command with every expected table, figure-source dataset and report section traceable to executed code.
- **SC-02:** All 120 requirements/tasks have a verified acceptance record or an explicit external prerequisite restriction; no unimplemented core capability is described as completed.
- **SC-03:** Golden numerical tests match independently invoked pinned reference packages at the tolerances in validation-strategy.md; design/contrast permutation invariance and exact enumeration tests pass.
- **SC-04:** Null/mixture calibration meets prespecified bounds without tuning seeds, thresholds or simulations after observing outcomes. Violations are investigated and block the affected validated-method claim.
- **SC-05:** A clean offline rerun with frozen resources reproduces semantic tables and figure data. Timestamps and cosmetic image metadata are distinguished from scientific nondeterminism.
- **SC-06:** The local archived-study regression reproduces source matrix/fold-change arithmetic, records legacy result exceptions and enumerates intentional differences in filtering, model, FDR and reversal reporting.
- **SC-07:** The medium benchmark (20,000 features × 100 observations, 8 contrasts, 2,000 gene sets) uses ≤8 GiB peak memory and completes QC+core limma+CAMERA+report within 30 minutes on documented 4-core hardware. Specialized engines/calibration are timed separately; inability to meet a target is reported, not hidden by smaller fixtures.
- **SC-08:** The verified successor is packaged on the SSD with a source manifest, changelog, tests and Spec Kit evidence; original historical files remain unchanged. No remote push is performed.

## Assumptions and Scope

Input is quantified protein-level abundance or a source export that can be losslessly mapped to it. LFQ DDA/DIA and multiplex TMT protein tables are supported when their declared contracts are met. The system does not establish peptide identification confidence from such matrices. Current mature statistical libraries remain the mathematical engines; the project adds scientifically constrained orchestration, validation and reporting rather than reimplementing empirical Bayes or enrichment algorithms. Execution begins only after independent specification audit, maintainer correction/freeze and explicit authorization.
