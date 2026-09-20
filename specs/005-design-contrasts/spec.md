# Feature Specification: Design validation, exact contrasts and blocking

**Feature Branch:** `spec/005-design-contrasts` (logical slice; stay on the current worktree branch unless a branch change is safe and needed).
**Created:** 2026-09-12. **Status:** Corrected draft awaiting Independent Reviewer audit and Maintainer freeze; not verified.
**Input:** Parent roadmap `specs/001-downstream-proteomics/roadmap.md` → **R04**. Deliver design validation, exact contrasts and blocking for user journey US2.

## User Scenarios & Testing

### US2 — Design validation, exact contrasts and blocking (Priority P1)

This slice delivers a usable and independently verifiable part of the parent user journey. It consumes the shared canonical contracts and exposes the behavior below through maintained package interfaces. It must not silently change the historical baseline or scientific interpretation.

**Why this priority:** dependent stages cannot reliably interpret their input without this contract being enforced.
**Independent Test:** run the acceptance cases below against actual implementations, including the negative cases; downstream slices may use controlled canonical fixtures rather than requiring an unfinished upstream feature.
**Dependencies:** R03 (`specs/004-preprocessing-qc/`).

## Requirements

- **FR-031**: The system MUST provide safe declarative design grammar. Groups, continuous covariates and interactions generate expected matrices; arbitrary function calls/code strings fail before R evaluation.
- **FR-032**: The system MUST provide stable levels and aligned matrices. Sample order shuffling with metadata alignment preserves coefficients; absent references and missing covariates give explicit errors.
- **FR-033**: The system MUST provide rank and confounding diagnostics. Batch=treatment fixture is detected, with aliased terms; no pseudoinverse silently manufactures an interpretable treatment effect.
- **FR-034**: The system MUST provide independent, paired and repeated plans. Subject-fixed and eligible duplicateCorrelation designs preserve the experimental unit and reject invalid/unidentified blocks.
- **FR-035**: The system MUST provide numeric contrasts and interaction semantics. Known disease/treatment/residual/interaction coefficients match analytic expectations; significance difference is never used as interaction.
- **FR-036**: The system MUST provide featurewise estimability and df. Sparse features get engine-specific n/df/reason records; ordinary available-case all-missing groups are nonestimable and eligible native-dropout estimates retain their prior-dependence caveat.
- **FR-037**: The system MUST provide exact contrast covariance strategy. Sparse nonorthogonal weighted designs match direct coefficient-refit SE/P; exact shortcut is restricted to proven eligible designs.
- **FR-038**: The system MUST provide method-design applicability registry. Unsupported proDA/DEqMS block or threshold combinations fail explicitly and no engine is silently substituted.
- **FR-039**: The system MUST provide frozen analysis plan and matrices. Plan saves configuration, inclusion masks, matrices, hypotheses, families and hashes before fitting; config changes invalidate the plan.
- **FR-040**: The system MUST provide design adversarial verification. Unequal n, singletons, missing factor levels, repeated-unit imbalance and reordered IDs have reference-backed expected behavior.

## Key Entities

Use entities/keys in `specs/001-downstream-proteomics/data-model.md`. All artifacts carry schema_version, run_id/plan_hash where applicable, source lineage and execution status. No alternate implicit identity or inference representation is permitted.

## Acceptance Scenarios

### V031: Safe declarative design grammar

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** groups, continuous covariates and interactions generate expected matrices; arbitrary function calls/code strings fail before R evaluation.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V032: Stable levels and aligned matrices

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** sample order shuffling with metadata alignment preserves coefficients; absent references and missing covariates give explicit errors.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V033: Rank and confounding diagnostics

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** batch=treatment fixture is detected, with aliased terms; no pseudoinverse silently manufactures an interpretable treatment effect.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V034: Independent, paired and repeated plans

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** subject-fixed and eligible duplicateCorrelation designs preserve the experimental unit and reject invalid/unidentified blocks.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V035: Numeric contrasts and interaction semantics

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** known disease/treatment/residual/interaction coefficients match analytic expectations; significance difference is never used as interaction.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V036: Featurewise estimability and df

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** sparse features get engine-specific n/df/reason records; ordinary available-case all-missing groups are nonestimable and eligible native-dropout estimates retain their prior-dependence caveat.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V037: Exact contrast covariance strategy

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** sparse nonorthogonal weighted designs match direct coefficient-refit SE/P; exact shortcut is restricted to proven eligible designs.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V038: Method-design applicability registry

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** unsupported proDA/DEqMS block or threshold combinations fail explicitly and no engine is silently substituted.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V039: Frozen analysis plan and matrices

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** plan saves configuration, inclusion masks, matrices, hypotheses, families and hashes before fitting; config changes invalidate the plan.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V040: Design adversarial verification

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** unequal n, singletons, missing factor levels, repeated-unit imbalance and reordered IDs have reference-backed expected behavior.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

## Edge Cases and Success Criteria

Test the named invalid/missing/empty/reordered/boundary cases without weakening the shared scientific contract. Every requirement must have a passing eligible-case test and a relevant explicit failure or boundary case where applicable. Preserve finite/NA distinctions and scientific eligibility reasons. Success requires real behavior, schema-valid outputs, independently rerun gates, reviewed diff and evidence linked to a commit; process exit alone is insufficient.

## Assumptions and Scope Boundary

Implement only the capabilities assigned above and the narrow helpers they require. Read the live constitution and parent contracts. An optional per-study analysis is still a mandatory implemented adapter when assigned here. Do not implement unrelated services, raw-MS analysis, private-data transmission or speculative interfaces. Use synthetic/public-permitted fixtures. Maintainer handles local private regression in the final slice.
