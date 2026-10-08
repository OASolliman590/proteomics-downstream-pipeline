# ADR 0011 — Methods extension and multi-study validation scope amendment

Status: scope decided by the Maintainer on 2026-10-08. Recorded by the Claude route acting under [ADR 0007](0007-claude-route-and-integration-amendments.md). The specification drafted under this ADR is pending the Maintainer's review; no packet is dispatched by this ADR. Decisions D-64 to D-67, D-69, D-70 and D-73 in [decisions.md](../../specs/001-downstream-proteomics/decisions.md).

## Context

- The paper-compatible analysis (a Welch t-test on log2 abundance, a raw-P rule with BH reported alongside) can only be reproduced inside the maintained pipeline if the Welch engine is a declared, labelled method. The primary engine stays limma ([ADR 0002](0002-primary-engine-limma.md)).
- Pathway enrichment exists (SM16–SM18) but cannot consume user-supplied offline libraries under the licence constraints of the commonly used collections.
- Calibration evidence (R14f) covers penalised logistic regression under one resampling scheme. The other classifiers, LOOCV and subject-blocked resampling that users request are not yet calibrated.
- The pipeline has been exercised on one private study. Cross-design evidence needs a multi-study framework with explicit confidentiality rules, and public datasets must not be fetched without approval.

## Decision

- Add slice **016-methods-extension-validation**, Phase 5 (v1.2-methods-validation), with dispatch units **R15a–R15d** and **FR-168–FR-194, T168–T194, V168–V194**, appended after the existing ranges. No identity is renumbered or retired.
- Bump the scientific methods contract to **v1.4.0** with **SM42–SM46**. SM01–SM41 are unchanged.
- **Welch engine (R15a).** `welch` is an optional engine for pairwise group contrasts on log2 data. It may be declared as the primary or a secondary engine. The default remains limma, a single primary model is still required, and the Welch result carries its own family. Every paper-compatible DEP rule (raw P below a cutoff) is labelled "not multiplicity-controlled". BH is always reported next to it.
- **Offline pathway libraries (R15b).** User-supplied local GMT libraries are read only from hashed files whose version and licence are recorded. KEGG and MSigDB files are never redistributed in the repository or in artifacts. Reactome, GO and WikiPathways files may be used under their recorded open licences. The public repository ships only synthetic fixtures.
- **Calibration extension (R15c).** The R14d calibration is extended to lasso, elastic net with inner-CV tuned lambda, linear and polynomial SVM, optional random forest, LOOCV outer resampling and subject-blocked units, under the null and under a planted signal. Gates are recorded before the first run and are not revised afterwards.
- **Multi-study validation (R15d).** A study-manifest format, five design slots (two-group label-free DDA; paired or repeated measures; TMT with batches; DIA; multi-group with a continuous phenotype in a non-human organism) and one explicit reference re-run slot for a three-group serum study analysed with Welch t. Studies are a mix of Maintainer-private local studies and public ProteomeXchange/PRIDE datasets. Public data are downloaded only after the Maintainer approves an accession shortlist; the approval, URL and SHA-256 are recorded and the data are stored outside the repository. Private studies stay local and only a non-identifying PASS/FAIL status may be committed.
- Phase 5 runs serially: R15a → R15b → R15c → R15d.
- Phase 6 ([ADR 0012](0012-visualization-outputs.md)) is independent and may be built before Phase 5; see PHASES.md.

## Consequences

- Shared-interface amendments are made by the Maintainer before R15a, not by allowlist exceptions. They are: the `welch` value in the primary-engine declaration, applicable only when declared (the single-primary rule of ADR 0002 is kept); an optional `methods_extension` block in `analysis.schema.json`; and R Suggests entries for e1071 (SVM), randomForest (optional) and glmnet (already required by R14d). Missing optional packages give NOT_RUN, never INAPPLICABLE.
- Every refusal introduced here follows [ADR 0010](0010-adaptive-post-de-policy.md): the invalid piece is refused as narrowly as possible, the valid alternative runs, and every adaptation is recorded.
- Welch results are a secondary analysis of the same data unless the Maintainer declares otherwise. A Welch P value is never presented as a limma P value, and vice versa.
- Pathway libraries are supplied by the user. The repository does not provide or download them.
- Private study outcomes are Maintainer-only. The repository holds no private data, identifiers or study numbers.

## Verification

V168–V194. All fixtures are synthetic. The R15c calibration gates are evaluated mechanically from recorded seeds. The Maintainer's private reference re-run and the private study runs are NOT_RUN for the implementer and recorded only as non-identifying PASS/FAIL statuses.
