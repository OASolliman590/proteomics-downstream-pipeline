# Feature Specification: Preprocessing, missingness and quality control

**Feature Branch:** `spec/004-preprocessing-qc` (logical slice; stay on the current worktree branch unless a branch change is safe and needed).
**Created:** 2026-09-12. **Status:** Corrected draft awaiting Independent Reviewer audit and Maintainer freeze; not verified.
**Input:** Parent roadmap `specs/001-downstream-proteomics/roadmap.md` → **R03**. Deliver preprocessing, missingness and quality control for user journey US1.

## User Scenarios & Testing

### US1 — Preprocessing, missingness and quality control (Priority P1)

This slice delivers a usable and independently verifiable part of the parent user journey. It consumes the shared canonical contracts and exposes the behavior below through maintained package interfaces. It must not silently change the historical baseline or scientific interpretation.

**Why this priority:** dependent stages cannot reliably interpret their input without this contract being enforced.
**Independent Test:** run the acceptance cases below against actual implementations, including the negative cases; downstream slices may use controlled canonical fixtures rather than requiring an unfinished upstream feature.
**Dependencies:** R02 (`specs/003-intake/`).

## Requirements

- **FR-021**: The system MUST provide preserve processed abundance by default. An already normalized log2 fixture remains unchanged; every optional transform has a new artifact and recorded policy.
- **FR-022**: The system MUST provide lfq normalization policies. Median/reference normalization factors match analytic fixtures; quantile is explicit sensitivity; missing reference coverage is not silently ignored.
- **FR-023**: The system MUST provide tmt plex and bridge policies. Two-plex synthetic bridges align loading offsets; missing bridge blocks bridge mode, while identifiable no-bridge plex-covariate inference is separately supported and confounding fails.
- **FR-024**: The system MUST provide contrast-aware coverage masks. Original observed coverage is computed at biological-unit grain; ordinary available-case all-missing groups are nonestimable while eligible proDA follows native dropout rules with explicit uncertainty.
- **FR-025**: The system MUST provide missingness mechanism diagnostics. Observation/group/abundance summaries distinguish observed versus previously imputed values and do not assert MCAR/MNAR from a threshold.
- **FR-026**: The system MUST provide pca and correlation diagnostics. PCA-only fill is isolated; constant/empty/small matrices return valid diagnostic states; variance and feature/sample identities match direct calculations.
- **FR-027**: The system MUST provide prespecified exclusions and watchlists. Outlier flags never auto-delete observations; a changed exclusion reason/policy creates a new plan and sensitivity lineage.
- **FR-028**: The system MUST provide explicit imputation sensitivities. Deterministic minimum, stochastic Gaussian and KNN have correct names/seeds/masks, remain secondary and never alter the primary matrix.
- **FR-029**: The system MUST provide detection-only exploratory endpoint. Independent two-group detection uses exact 2x2 counts and a separate BH family; paired/complex detection requests are marked unsupported.
- **FR-030**: The system MUST provide qc source data and stage report. All plots have exact source tables; QC-only run completes without R modeling or any significant features.

## Key Entities

Use entities/keys in `specs/001-downstream-proteomics/data-model.md`. All artifacts carry schema_version, run_id/plan_hash where applicable, source lineage and execution status. No alternate implicit identity or inference representation is permitted.

## Acceptance Scenarios

### V021: Preserve processed abundance by default

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** an already normalized log2 fixture remains unchanged; every optional transform has a new artifact and recorded policy.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V022: LFQ normalization policies

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** median/reference normalization factors match analytic fixtures; quantile is explicit sensitivity; missing reference coverage is not silently ignored.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V023: TMT plex and bridge policies

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** two-plex synthetic bridges align loading offsets; missing bridge blocks bridge mode, while identifiable no-bridge plex-covariate inference is separately supported and confounding fails.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V024: Contrast-aware coverage masks

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** original observed coverage is computed at biological-unit grain; ordinary available-case all-missing groups are nonestimable while eligible proDA follows native dropout rules with explicit uncertainty.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V025: Missingness mechanism diagnostics

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** observation/group/abundance summaries distinguish observed versus previously imputed values and do not assert MCAR/MNAR from a threshold.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V026: PCA and correlation diagnostics

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** pCA-only fill is isolated; constant/empty/small matrices return valid diagnostic states; variance and feature/sample identities match direct calculations.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V027: Prespecified exclusions and watchlists

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** outlier flags never auto-delete observations; a changed exclusion reason/policy creates a new plan and sensitivity lineage.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V028: Explicit imputation sensitivities

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** deterministic minimum, stochastic Gaussian and KNN have correct names/seeds/masks, remain secondary and never alter the primary matrix.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V029: Detection-only exploratory endpoint

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** independent two-group detection uses exact 2x2 counts and a separate BH family; paired/complex detection requests are marked unsupported.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V030: QC source data and stage report

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** all plots have exact source tables; QC-only run completes without R modeling or any significant features.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

## Edge Cases and Success Criteria

Test the named invalid/missing/empty/reordered/boundary cases without weakening the shared scientific contract. Every requirement must have a passing eligible-case test and a relevant explicit failure or boundary case where applicable. Preserve finite/NA distinctions and scientific eligibility reasons. Success requires real behavior, schema-valid outputs, independently rerun gates, reviewed diff and evidence linked to a commit; process exit alone is insufficient.

## Assumptions and Scope Boundary

Implement only the capabilities assigned above and the narrow helpers they require. Read the live constitution and parent contracts. An optional per-study analysis is still a mandatory implemented adapter when assigned here. Do not implement unrelated services, raw-MS analysis, private-data transmission or speculative interfaces. Use synthetic/public-permitted fixtures. Maintainer handles local private regression in the final slice.
