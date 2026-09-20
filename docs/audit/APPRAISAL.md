# Appraisal of the recovered proteomics pipeline

Audit date: 2026-09-12. Decision: **retain as an auditable exploratory study pipeline; do not describe it as a validated general proteomics platform or confirmatory rescue analysis.** The recoverable code is useful, and the core numerical exports are internally consistent. Scientific inference and full computational reproducibility require further work.

## What the evidence supports

The workbook-derived matrix has 3,714 proteins and 20 samples, four per recorded group. The independent audit matched workbook values to the exported abundance matrix, matched Methods A and B, verified the outer archive manifest and recomputed all stored refit fold changes and BH protein adjustments to numerical precision. No protein passes FDR 0.05 in the four principal disease/direct-treatment contrasts. No formal rescue is reported under the historical rule.

The archive's universal “no FDR-supported proteins under the tested contrasts” statement is too broad. Joint-model **Negr1 / A0A8I6AJV3** is significant in PreDM+dapa versus CTRL: log2FC −1.07835, P=4.5723×10⁻⁶, q=0.0169815. This is an auxiliary sensitivity result and does not demonstrate drug rescue. The joint direct PreDM-treatment contrast has one protein at q<0.10 (minimum 0.07936); the corresponding axis model has none (minimum 0.11490).

| Dimension | Appraisal | Basis |
|---|---|---|
| Data preservation | Strong within this archive | Original ZIP retained; all 164 entries in the outer manifest match. |
| Core downstream statistics | Defensible with limitations | Explicit group contrasts, empirical-Bayes variance moderation, no extra log transform, no blanket primary imputation; exported arithmetic verified. |
| Experimental provenance | Incomplete | Rat independence and treatment identity are asserted in later files; original confirmation, tissue and important design metadata absent. |
| Reversal claims | Exploratory only | Selection/testing reuse, shared untreated-group noise, ratio instability and overshoot labels. |
| Pathway claims | Hypothesis-generating | Preranked enrichment has a gene-set null; correlation sensitivity, mapping provenance and numerical diagnostics are incomplete. |
| Exact computational reproduction | Not established | Broken historical export layout, incomplete/unpinned dependencies, mutable gene sets; no local R rerun. |
| General reuse | Requires redesign | Fixed groups, sample counts, workbook names, thresholds and study assertions are embedded in code. |
| GitHub preparation | Isolated source candidate prepared | Code/data boundary, provenance, docs, Python tests and portable launcher prepared. License/attribution and release scope unresolved. |

## Priority findings

1. **Correct the biological headline.** Scope the no-DEP statement to the principal contrasts and show the Negr1 sensitivity exception. Do not equate lack of protein significance with absence of an effect, or differences in significance with a stage interaction.
2. **Constrain reversal interpretation.** Existing selected-score P values are not valid confirmatory tests. Fixing the exhaustive permutation denominator alone does not address selection. Classify overshoot separately and define restoration to control before making rescue claims.
3. **Recover experimental and computational provenance.** Obtain original sample confirmation and upstream processing details. Capture an actual R dependency lock and exact gene-set snapshot. Separate archival reconstruction from a scientifically revised implementation.
4. **Validate pathway robustness.** State the collection/contrast testing family; assess consequential findings with correlation-aware analysis, mapping alternatives and the chosen primary model. Selected-pathway leave-one-out is descriptive stability, not independent validation.
5. **Complete a clean R rerun before claiming reproducibility.** Compare all key outputs, save per-stage loaded package versions and numerical diagnostics, and document expected changes. Resolve downstream reporting failures when optional results are absent.

## Validation achieved and its limits

Method A executed successfully from the workbook. Its reproduced abundance, annotations, standardized DEA and validation tables were compared with the archive. The new preparation layer preserved the Method B values and missingness, and all ten Python validation/failure-path tests passed. Python sources compiled. Preflight validated the study input and correctly reported Rscript unavailable. The generated CI configuration has not run on GitHub.

Independent checks verify arithmetic and integrity; they do not re-estimate limma moderated variances/P values, rerun fgsea, establish rat identity, or validate the original MS measurements. No statistical repair was silently applied to the historical R stages. The [scientific audit](SCIENTIFIC_AUDIT.md), [technical audit](TECHNICAL_AUDIT.md), [validation record](VALIDATION.md) and [release checklist](../RELEASE_CHECKLIST.md) give the next steps.
