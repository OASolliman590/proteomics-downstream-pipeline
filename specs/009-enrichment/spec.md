# Feature Specification: Design-compatible pathways and enrichment

**Feature Branch:** `spec/009-enrichment` (logical slice; stay on the current worktree branch unless a branch change is safe and needed).
**Created:** 2026-09-12. **Status:** Corrected draft awaiting Independent Reviewer audit and Maintainer freeze; not verified.
**Input:** Parent roadmap `specs/001-downstream-proteomics/roadmap.md` → **R08**. Deliver design-compatible pathways and enrichment for user journey US3.

## User Scenarios & Testing

### US3 — Design-compatible pathways and enrichment (Priority P1)

This slice delivers a usable and independently verifiable part of the parent user journey. It consumes the shared canonical contracts and exposes the behavior below through maintained package interfaces. It must not silently change the historical baseline or scientific interpretation.

**Why this priority:** dependent stages cannot reliably interpret their input without this contract being enforced.
**Independent Test:** run the acceptance cases below against actual implementations, including the negative cases; downstream slices may use controlled canonical fixtures rather than requiring an unfinished upstream feature.
**Dependencies:** R06 (`specs/007-assay-engines/`) and R07 (`specs/008-resources-mapping/`).

## Requirements

- **FR-071**: The system MUST provide pathway null and method dispatcher. Plan distinguishes competitive, rotation self-contained, preranked and ORA nulls and rejects unsupported method/design combinations.
- **FR-072**: The system MUST provide independent-design camera. Finite gene matrix/design/contrast match official CAMERA with estimated correlation; block/correlation arguments cannot be ignored.
- **FR-073**: The system MUST provide blocked-design roast. Fixed-subject pairing passes its design only; duplicateCorrelation passes block/correlation without duplicate subject encoding; midpFALSE/rotation settings and reference outputs match.
- **FR-074**: The system MUST provide directional and mixed rotation endpoints. Directional and mixed raw P/statistics/families remain distinct; central adjustment uses correct columns rather than package-native FDR mislabeled global.
- **FR-075**: The system MUST provide exploratory fgsea adapter. Declared native ranks, ties, seeds and eps settings are honored; NES/log2err/warnings/leading edges and nonfinite tests are retained.
- **FR-076**: The system MUST provide exact background-aware ora. Hand-enumerated hypergeometric examples match; all eligible sets enter the family and empty foreground gives a valid explanation.
- **FR-077**: The system MUST provide pathway family adjustment. Collections/contrasts/null types/primary status obey the frozen families; original and central q fields remain distinguishable.
- **FR-078**: The system MUST provide alternative-engine pathway sensitivity. proDA/DEqMS ranks and linear gene-model sensitivity are labeled correctly and never sold as equivalent correlation-adjusted primary-engine inference.
- **FR-079**: The system MUST provide leading-edge redundancy and diagnostics. Overlap relationships preserve original set IDs/P values and do not merge independent discoveries or invent causal mechanisms.
- **FR-080**: The system MUST provide enrichment golden and failure suite. Exact ORA, reference CAMERA/ROAST/fgsea, missing snapshots, NA/Inf matrix, all nonfinite tests and no discoveries exercise honest statuses.

## Key Entities

Use entities/keys in `specs/001-downstream-proteomics/data-model.md`. All artifacts carry schema_version, run_id/plan_hash where applicable, source lineage and execution status. No alternate implicit identity or inference representation is permitted.

## Acceptance Scenarios

### V071: Pathway null and method dispatcher

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** plan distinguishes competitive, rotation self-contained, preranked and ORA nulls and rejects unsupported method/design combinations.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V072: Independent-design CAMERA

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** finite gene matrix/design/contrast match official CAMERA with estimated correlation; block/correlation arguments cannot be ignored.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V073: Blocked-design ROAST

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** fixed-subject pairing passes its design only; duplicateCorrelation passes block/correlation without duplicate subject encoding; midpFALSE/rotation settings and reference outputs match.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V074: Directional and mixed rotation endpoints

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** directional and mixed raw P/statistics/families remain distinct; central adjustment uses correct columns rather than package-native FDR mislabeled global.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V075: Exploratory fgsea adapter

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** declared native ranks, ties, seeds and eps settings are honored; NES/log2err/warnings/leading edges and nonfinite tests are retained.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V076: Exact background-aware ORA

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** hand-enumerated hypergeometric examples match; all eligible sets enter the family and empty foreground gives a valid explanation.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V077: Pathway family adjustment

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** collections/contrasts/null types/primary status obey the frozen families; original and central q fields remain distinguishable.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V078: Alternative-engine pathway sensitivity

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** proDA/DEqMS ranks and linear gene-model sensitivity are labeled correctly and never sold as equivalent correlation-adjusted primary-engine inference.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V079: Leading-edge redundancy and diagnostics

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** overlap relationships preserve original set IDs/P values and do not merge independent discoveries or invent causal mechanisms.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V080: Enrichment golden and failure suite

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** exact ORA, reference CAMERA/ROAST/fgsea, missing snapshots, NA/Inf matrix, all nonfinite tests and no discoveries exercise honest statuses.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

## Edge Cases and Success Criteria

Test the named invalid/missing/empty/reordered/boundary cases without weakening the shared scientific contract. Every requirement must have a passing eligible-case test and a relevant explicit failure or boundary case where applicable. Preserve finite/NA distinctions and scientific eligibility reasons. Success requires real behavior, schema-valid outputs, independently rerun gates, reviewed diff and evidence linked to a commit; process exit alone is insufficient.

## Assumptions and Scope Boundary

Implement only the capabilities assigned above and the narrow helpers they require. Read the live constitution and parent contracts. An optional per-study analysis is still a mandatory implemented adapter when assigned here. Do not implement unrelated services, raw-MS analysis, private-data transmission or speculative interfaces. Use synthetic/public-permitted fixtures. Maintainer handles local private regression in the final slice.
