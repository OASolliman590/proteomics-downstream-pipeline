# ADR 0009 — Post-differential analysis scope amendment

Status: accepted by the operator (Omar, Maintainer) on 2026-10-02, with decisions: include R14e, penalised logistic as the default classifier, Phase 4 placement. Recorded by the Claude (Opus) route acting for the Maintainer role (ADR 0007). This is an explicit scope amendment, not a silent expansion. Implementation has not started.

## Context

The frozen kit stops after differential abundance, enrichment, treatment response and multivariate PERMANOVA. Real analyses need several further steps:
- structuring results across contrasts (declared set logic, concordance);
- stress-testing them against covariates, subgroups, sample size and influential units;
- relating proteins to continuous phenotypes;
- evaluating candidate biomarkers without selection leakage;
- co-abundance and interaction-network context.

The operator's own reviewer-driven work showed how fragile these are when done ad hoc:
- an undocumented set rule (108 vs 133);
- fixed-panel LOOCV presented as validation;
- sex confounding answered only by subgroup significance counts.

The spec currently excludes "classifier training" (spec 001, line 7). The operator asked that the post-DE layer fit any dataset, not one study.

## Decision

- Add slice **015-post-de-analysis**, Phase 4 (v1.1-post-de), with dispatch units **R14a–R14f** and **FR-131–FR-167, T131–T167, V131–V167**, appended after the existing ranges. No identity is renumbered.
- Bump the scientific methods contract to **v1.3.0** with **SM31–SM41**. SM01–SM30 are unchanged.
- Amend spec 001's scope sentence:
  - Replace "classifier training" with "deployment of clinical classifiers or diagnostic/clinical-utility claims; leakage-safe biomarker discrimination evaluation with explicit claim labels is in scope (US8/US9)".
  - Add "interaction-network analysis only from offline versioned snapshots". Network *deployment* stays excluded.
- Add US8 (result structure, robustness, phenotype association) and US9 (biomarker evaluation, co-abundance and networks) to the user journeys.
- Defaults (operator-confirmed or drafted defaults left unchanged):
  - classifier: penalised logistic regression;
  - minimum 5 biological units per class for biomarker evaluation;
  - minimum 20 units for co-abundance;
  - compute guard: 6 hours on the reference machine.
- Generality is normative (SM41). Configuration is declarative, every module declares its eligibility, and a synthetic design matrix (V167) qualifies every module across design types, sample sizes, imbalance and species.

## Consequences

- The kit grows to 167 requirements across 14 slices. `scripts/check_spec_kit.py` must check the new counts and that FR/T/V 001–130 are unchanged.
- Biomarker outputs are performance *estimates* with claim labels, not discoveries. They have no FDR family.
- All six packets, including R14e (co-abundance and PPI), are in v1.1 by operator decision. R14e's modules are opt-in at run time and refused below the declared sample-size minimum.
- Compute-heavy procedures (nested CV × whole-procedure permutation) refuse beyond a declared guard rather than silently cutting work.
- Missing optional packages (glmnet, randomForest, WGCNA) are NOT_RUN, never INAPPLICABLE.

## Verification

V131–V167. The R11 calibration evidence is extended with null calibration of nested CV and permutation, a leaky-reference inflation demonstration, and calibration of the network-connectivity null. All fixtures are synthetic.
