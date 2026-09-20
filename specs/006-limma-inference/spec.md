# Feature Specification: Core limma inference and multiplicity

**Feature Branch:** `spec/006-limma-inference` (logical slice; stay on the current worktree branch unless a branch change is safe and needed).
**Created:** 2026-09-12. **Status:** Corrected draft awaiting Independent Reviewer audit and Maintainer freeze; not verified.
**Input:** Parent roadmap `specs/001-downstream-proteomics/roadmap.md` → **R05**. Deliver core limma inference and multiplicity for user journey US2.

## User Scenarios & Testing

### US2 — Core limma inference and multiplicity (Priority P1)

This slice delivers a usable and independently verifiable part of the parent user journey. It consumes the shared canonical contracts and exposes the behavior below through maintained package interfaces. It must not silently change the historical baseline or scientific interpretation.

**Why this priority:** dependent stages cannot reliably interpret their input without this contract being enforced.
**Independent Test:** run the acceptance cases below against actual implementations, including the negative cases; downstream slices may use controlled canonical fixtures rather than requiring an unfinished upstream feature.
**Dependencies:** R04 (`specs/005-design-contrasts/`).

## Requirements

- **FR-041**: The system MUST provide observed-data limma adapter. Eligible fixtures match direct lmFit/eBayes calls for effects/SE/df/P; primary input retains NA and declared normalization.
- **FR-042**: The system MUST provide exact sparse and weighted contrasts. General designs use validated contrast reparametrization/direct covariance; nonorthogonal sparse example matches a separate oracle.
- **FR-043**: The system MUST provide robust/trend diagnostics. Settings, prior/posterior variance and trend diagnostics are saved; robust eBayes is not mislabeled robust sample regression.
- **FR-044**: The system MUST provide effect-threshold treat endpoint. Boundary effects compare to direct TREAT; true threshold-null evidence is distinct from observed effect filtering and zero-null results.
- **FR-045**: The system MUST provide backend-correct confidence intervals. Intervals use correct moderated variance/df, match independent calculations and preserve unbounded/unavailable cases with reason.
- **FR-046**: The system MUST provide declared multiple-testing families. Finite-P BH/BY matches analytic and R references; separate endpoints/primary/sensitivity sets do not accidentally pool or omit rows.
- **FR-047**: The system MUST provide omnibus and interaction outputs. F tests have distinct hypotheses/families, complete planned contrasts are exported and no unadjusted post-hoc claims appear.
- **FR-048**: The system MUST provide complete typed differential tables. Every planned feature-contrast row is tested or explicitly excluded/nonestimable; no zero/NA ambiguity or mislabeled q fields.
- **FR-049**: The system MUST provide biological-unit influence diagnostics. Whole subjects/technical sets are omitted appropriately; relevant omission counts are correct and influence stays descriptive.
- **FR-050**: The system MUST provide core inference golden suite. Complete/sparse/paired/blocked fixtures match pinned direct references and no-null/no-discovery result paths pass.

## Key Entities

Use entities/keys in `specs/001-downstream-proteomics/data-model.md`. All artifacts carry schema_version, run_id/plan_hash where applicable, source lineage and execution status. No alternate implicit identity or inference representation is permitted.

## Acceptance Scenarios

### V041: Observed-data limma adapter

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** eligible fixtures match direct lmFit/eBayes calls for effects/SE/df/P; primary input retains NA and declared normalization.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V042: Exact sparse and weighted contrasts

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** general designs use validated contrast reparametrization/direct covariance; nonorthogonal sparse example matches a separate oracle.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V043: Robust/trend diagnostics

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** settings, prior/posterior variance and trend diagnostics are saved; robust eBayes is not mislabeled robust sample regression.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V044: Effect-threshold TREAT endpoint

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** boundary effects compare to direct TREAT; true threshold-null evidence is distinct from observed effect filtering and zero-null results.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V045: Backend-correct confidence intervals

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** intervals use correct moderated variance/df, match independent calculations and preserve unbounded/unavailable cases with reason.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V046: Declared multiple-testing families

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** finite-P BH/BY matches analytic and R references; separate endpoints/primary/sensitivity sets do not accidentally pool or omit rows.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V047: Omnibus and interaction outputs

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** f tests have distinct hypotheses/families, complete planned contrasts are exported and no unadjusted post-hoc claims appear.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V048: Complete typed differential tables

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** every planned feature-contrast row is tested or explicitly excluded/nonestimable; no zero/NA ambiguity or mislabeled q fields.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V049: Biological-unit influence diagnostics

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** whole subjects/technical sets are omitted appropriately; relevant omission counts are correct and influence stays descriptive.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V050: Core inference golden suite

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** complete/sparse/paired/blocked fixtures match pinned direct references and no-null/no-discovery result paths pass.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

## Edge Cases and Success Criteria

Test the named invalid/missing/empty/reordered/boundary cases without weakening the shared scientific contract. Every requirement must have a passing eligible-case test and a relevant explicit failure or boundary case where applicable. Preserve finite/NA distinctions and scientific eligibility reasons. Success requires real behavior, schema-valid outputs, independently rerun gates, reviewed diff and evidence linked to a commit; process exit alone is insufficient.

## Assumptions and Scope Boundary

Implement only the capabilities assigned above and the narrow helpers they require. Read the live constitution and parent contracts. An optional per-study analysis is still a mandatory implemented adapter when assigned here. Do not implement unrelated services, raw-MS analysis, private-data transmission or speculative interfaces. Use synthetic/public-permitted fixtures. Maintainer handles local private regression in the final slice.
