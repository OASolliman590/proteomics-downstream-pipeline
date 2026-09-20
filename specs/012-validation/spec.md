# Feature Specification: Reproducibility, calibration and continuous validation

**Feature Branch:** `spec/012-validation` (logical slice; stay on the current worktree branch unless a branch change is safe and needed).
**Created:** 2026-09-12. **Status:** Corrected draft awaiting Independent Reviewer audit and Maintainer freeze; not verified.
**Input:** Parent roadmap `specs/001-downstream-proteomics/roadmap.md` → **R11**. Deliver reproducibility, calibration and continuous validation for user journey US6.

## User Scenarios & Testing

### US6 — Reproducibility, calibration and continuous validation (Priority P1)

This slice delivers a usable and independently verifiable part of the parent user journey. It consumes the shared canonical contracts and exposes the behavior below through maintained package interfaces. It must not silently change the historical baseline or scientific interpretation.

**Why this priority:** dependent stages cannot reliably interpret their input without this contract being enforced.
**Independent Test:** run the acceptance cases below against actual implementations, including the negative cases; downstream slices may use controlled canonical fixtures rather than requiring an unfinished upstream feature.
**Dependencies:** R10 (`specs/011-reporting/`).

## Requirements

- **FR-101**: The system MUST provide pinned compatible environments. Clean Python/R/Bioconductor restore is tested; version pair is compatible and installed package APIs match adapters.
- **FR-102**: The system MUST provide immutable resource-offline reproduction. A run with cached snapshots completes without network and semantic output hashes match the documented tolerance.
- **FR-103**: The system MUST provide resume and cache invalidation. Code/config/data/resource changes invalidate dependent stages; compatible verified outputs reuse safely and partial artifacts never count as cache.
- **FR-104**: The system MUST provide scientific golden reference matrix. Independent calls cover every backend/hypothesis/design eligibility; tolerances/versions are frozen before candidate comparison.
- **FR-105**: The system MUST provide null and mixture calibration. Prespecified all-null/mixture scenarios calculate actual error/coverage metrics with Monte Carlo uncertainty; failing methods are not labeled validated.
- **FR-106**: The system MUST provide correlated pathway calibration. Known correlated null and spike-in pathways assess method-null-specific error/power; fgsea caveat is retained rather than demanding a false universal guarantee.
- **FR-107**: The system MUST provide cross-platform ci. Windows/Linux core Python+R jobs execute real examples; expensive release calibration is distinct and not silently skipped as passing CI.
- **FR-108**: The system MUST provide performance and resource benchmark. 20k x100 medium case logs time/peak memory/hardware, specialized engine cost separately; targets are assessed on actual work.
- **FR-109**: The system MUST provide report numerical/visual validation. Plot-data values and references match tables; representative rendered reports/figures are visually inspected with explicit evidence.
- **FR-110**: The system MUST provide machine-readable validation ledger. Every requirement/test records command, versions, artifact hash, result and limitations; NOT_RUN/SKIPPED is not PASS.

## Key Entities

Use entities/keys in `specs/001-downstream-proteomics/data-model.md`. All artifacts carry schema_version, run_id/plan_hash where applicable, source lineage and execution status. No alternate implicit identity or inference representation is permitted.

## Acceptance Scenarios

### V101: Pinned compatible environments

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** clean Python/R/Bioconductor restore is tested; version pair is compatible and installed package APIs match adapters.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V102: Immutable resource-offline reproduction

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** a run with cached snapshots completes without network and semantic output hashes match the documented tolerance.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V103: Resume and cache invalidation

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** code/config/data/resource changes invalidate dependent stages; compatible verified outputs reuse safely and partial artifacts never count as cache.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V104: Scientific golden reference matrix

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** independent calls cover every backend/hypothesis/design eligibility; tolerances/versions are frozen before candidate comparison.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V105: Null and mixture calibration

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** prespecified all-null/mixture scenarios calculate actual error/coverage metrics with Monte Carlo uncertainty; failing methods are not labeled validated.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V106: Correlated pathway calibration

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** known correlated null and spike-in pathways assess method-null-specific error/power; fgsea caveat is retained rather than demanding a false universal guarantee.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V107: Cross-platform CI

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** windows/Linux core Python+R jobs execute real examples; expensive release calibration is distinct and not silently skipped as passing CI.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V108: Performance and resource benchmark

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** 20k x100 medium case logs time/peak memory/hardware, specialized engine cost separately; targets are assessed on actual work.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V109: Report numerical/visual validation

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** plot-data values and references match tables; representative rendered reports/figures are visually inspected with explicit evidence.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V110: Machine-readable validation ledger

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** every requirement/test records command, versions, artifact hash, result and limitations; NOT_RUN/SKIPPED is not PASS.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

## Edge Cases and Success Criteria

Test the named invalid/missing/empty/reordered/boundary cases without weakening the shared scientific contract. Every requirement must have a passing eligible-case test and a relevant explicit failure or boundary case where applicable. Preserve finite/NA distinctions and scientific eligibility reasons. Success requires real behavior, schema-valid outputs, independently rerun gates, reviewed diff and evidence linked to a commit; process exit alone is insufficient.

## Assumptions and Scope Boundary

Implement only the capabilities assigned above and the narrow helpers they require. Read the live constitution and parent contracts. An optional per-study analysis is still a mandatory implemented adapter when assigned here. Do not implement unrelated services, raw-MS analysis, private-data transmission or speculative interfaces. Use synthetic/public-permitted fixtures. Maintainer handles local private regression in the final slice.
