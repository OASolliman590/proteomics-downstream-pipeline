# Feature Specification: Treatment response, equivalence and independent scores

**Feature Branch:** `spec/010-treatment-response` (logical slice; stay on the current worktree branch unless a branch change is safe and needed).
**Created:** 2026-09-12. **Status:** Corrected draft awaiting Independent Reviewer audit and Maintainer freeze; not verified.
**Input:** Parent roadmap `specs/001-downstream-proteomics/roadmap.md` → **R09**. Deliver treatment response, equivalence and independent scores for user journey US4.

## User Scenarios & Testing

### US4 — Treatment response, equivalence and independent scores (Priority P1)

This slice delivers a usable and independently verifiable part of the parent user journey. It consumes the shared canonical contracts and exposes the behavior below through maintained package interfaces. It must not silently change the historical baseline or scientific interpretation.

**Why this priority:** dependent stages cannot reliably interpret their input without this contract being enforced.
**Independent Test:** run the acceptance cases below against actual implementations, including the negative cases; downstream slices may use controlled canonical fixtures rather than requiring an unfinished upstream feature.
**Dependencies:** R05 (`specs/006-limma-inference/`).

## Requirements

- **FR-081**: The system MUST provide explicit disease-treatment-residual axes. d, t, r and interaction estimates align to planned contrasts with covariance; r=d+t matches analytic fixtures.
- **FR-082**: The system MUST provide descriptive reversal without double dipping. Current-data selected features can be displayed but produce no confirmatory score P or rescue label; selection provenance is explicit.
- **FR-083**: The system MUST provide bounded response categories. Configured epsilon classes handle zero/small d, RI=1, RI>1+epsilon and RI=3; old unlimited full-reversal label is absent.
- **FR-084**: The system MUST provide ratio uncertainty contract. Supported joint-covariance method handles unstable/unbounded RI intervals; dividing CI endpoints is rejected and descriptive-only mode is explicit.
- **FR-085**: The system MUST provide model-correct residual equivalence. Analytic TOST examples match both one-sided tests/max-P and 90% equivalence intervals; a nonsignificant zero-null test alone cannot pass.
- **FR-086**: The system MUST provide conjunction claim eligibility. Only independently defined direction and valid model-supported conjunction yield adjusted restoration endpoints; FDR-list intersection is not labeled rescue FDR.
- **FR-087**: The system MUST provide independent fixed-score artifacts. Training identity/weights/transforms are frozen and checked; validation-outcome-derived selection or refitting produces descriptive-only eligibility.
- **FR-088**: The system MUST provide exact and monte carlo randomization. All 70 four-versus-four allocations use k/N; Monte Carlo +1 and interval are distinct; paired/restricted permutations honor blocks and ties.
- **FR-089**: The system MUST provide score/equivalence reporting semantics. Report distinguishes movement, overshoot, equivalence, conjunction evidence and unknown training overlap; no causal rescue prose is generated.
- **FR-090**: The system MUST provide circularity/null/overshoot tests. Shared-untreated null exposes negative effect covariance but cannot yield a spurious inferential score path; exact old examples and RI overshoot regressions pass.

## Key Entities

Use entities/keys in `specs/001-downstream-proteomics/data-model.md`. All artifacts carry schema_version, run_id/plan_hash where applicable, source lineage and execution status. No alternate implicit identity or inference representation is permitted.

## Acceptance Scenarios

### V081: Explicit disease-treatment-residual axes

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** d, t, r and interaction estimates align to planned contrasts with covariance; r=d+t matches analytic fixtures.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V082: Descriptive reversal without double dipping

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** current-data selected features can be displayed but produce no confirmatory score P or rescue label; selection provenance is explicit.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V083: Bounded response categories

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** configured epsilon classes handle zero/small d, RI=1, RI>1+epsilon and RI=3; old unlimited full-reversal label is absent.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V084: Ratio uncertainty contract

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** supported joint-covariance method handles unstable/unbounded RI intervals; dividing CI endpoints is rejected and descriptive-only mode is explicit.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V085: Model-correct residual equivalence

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** analytic TOST examples match both one-sided tests/max-P and 90% equivalence intervals; a nonsignificant zero-null test alone cannot pass.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V086: Conjunction claim eligibility

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** only independently defined direction and valid model-supported conjunction yield adjusted restoration endpoints; FDR-list intersection is not labeled rescue FDR.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V087: Independent fixed-score artifacts

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** training identity/weights/transforms are frozen and checked; validation-outcome-derived selection or refitting produces descriptive-only eligibility.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V088: Exact and Monte Carlo randomization

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** all 70 four-versus-four allocations use k/N; Monte Carlo +1 and interval are distinct; paired/restricted permutations honor blocks and ties.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V089: Score/equivalence reporting semantics

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** report distinguishes movement, overshoot, equivalence, conjunction evidence and unknown training overlap; no causal rescue prose is generated.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V090: Circularity/null/overshoot tests

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** shared-untreated null exposes negative effect covariance but cannot yield a spurious inferential score path; exact old examples and RI overshoot regressions pass.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

## Edge Cases and Success Criteria

Test the named invalid/missing/empty/reordered/boundary cases without weakening the shared scientific contract. Every requirement must have a passing eligible-case test and a relevant explicit failure or boundary case where applicable. Preserve finite/NA distinctions and scientific eligibility reasons. Success requires real behavior, schema-valid outputs, independently rerun gates, reviewed diff and evidence linked to a commit; process exit alone is insufficient.

## Assumptions and Scope Boundary

Implement only the capabilities assigned above and the narrow helpers they require. Read the live constitution and parent contracts. An optional per-study analysis is still a mandatory implemented adapter when assigned here. Do not implement unrelated services, raw-MS analysis, private-data transmission or speculative interfaces. Use synthetic/public-permitted fixtures. Maintainer handles local private regression in the final slice.
