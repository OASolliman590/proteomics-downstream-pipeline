# Feature Specification: Assay-qualified DEqMS and proDA backends

**Feature Branch:** `spec/007-assay-engines` (logical slice; stay on the current worktree branch unless a branch change is safe and needed).
**Created:** 2026-09-12. **Status:** Corrected draft awaiting Independent Reviewer audit and Maintainer freeze; not verified.
**Input:** Parent roadmap `specs/001-downstream-proteomics/roadmap.md` → **R06**. Deliver assay-qualified deqms and proda backends for user journey US2.

## User Scenarios & Testing

### US2 — Assay-qualified DEqMS and proDA backends (Priority P1)

This slice delivers a usable and independently verifiable part of the parent user journey. It consumes the shared canonical contracts and exposes the behavior below through maintained package interfaces. It must not silently change the historical baseline or scientific interpretation.

**Why this priority:** dependent stages cannot reliably interpret their input without this contract being enforced.
**Independent Test:** run the acceptance cases below against actual implementations, including the negative cases; downstream slices may use controlled canonical fixtures rather than requiring an unfinished upstream feature.
**Dependencies:** R05 (`specs/006-limma-inference/`).

## Requirements

- **FR-051**: The system MUST provide count evidence schema. True peptide/PSM count grain and provenance validate; observed-sample count or abundance proxy is rejected as unsupported evidence.
- **FR-052**: The system MUST provide deqms fit pipeline. Installed official fit/eBayes/count/spectraCounteBayes sequence reproduces reference count-adjusted statistics.
- **FR-053**: The system MUST provide deqms statistic alignment. sca P/statistic/variance fields are kept separate from original limma fields; shuffled count IDs align or fail, never silently reorder.
- **FR-054**: The system MUST provide deqms hypotheses and design guard. Unsupported block/precision-weight/threshold-null requests are rejected with explicit contract; eligible sparse exact-contrast use is verified.
- **FR-055**: The system MUST provide lfq dropout provenance gate. proDA accepts qualified unimputed LFQ data and rejects TMT, unknown incompatible imputation and unsupported repeated models.
- **FR-056**: The system MUST provide proda fitting and test_diff. Zero-null/reduced-model reference calls match native output and retain model-specific uncertainty and dropout diagnostics.
- **FR-057**: The system MUST provide alternative-engine uncertainty. Engine-appropriate intervals/statistics are verified against native/reference calculations; no borrowed limma df or fake TREAT output.
- **FR-058**: The system MUST provide method sensitivity comparison. Matched-universe effects/discoveries are displayed as sensitivities, with unmatched universes explained; no hit-count winner changes primary.
- **FR-059**: The system MUST provide dependencies and failure propagation. Missing package/backend error is visible, required analysis fails, optional applicability is distinguished from installation failure.
- **FR-060**: The system MUST provide assay-engine golden/calibration fixtures. Count-dependent variance and LFQ dropout simulations have pinned reference comparisons, unsuitable-input negative cases and recorded diagnostics.

## Key Entities

Use entities/keys in `specs/001-downstream-proteomics/data-model.md`. All artifacts carry schema_version, run_id/plan_hash where applicable, source lineage and execution status. No alternate implicit identity or inference representation is permitted.

## Acceptance Scenarios

### V051: Count evidence schema

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** true peptide/PSM count grain and provenance validate; observed-sample count or abundance proxy is rejected as unsupported evidence.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V052: DEqMS fit pipeline

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** installed official fit/eBayes/count/spectraCounteBayes sequence reproduces reference count-adjusted statistics.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V053: DEqMS statistic alignment

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** sca P/statistic/variance fields are kept separate from original limma fields; shuffled count IDs align or fail, never silently reorder.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V054: DEqMS hypotheses and design guard

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** unsupported block/precision-weight/threshold-null requests are rejected with explicit contract; eligible sparse exact-contrast use is verified.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V055: LFQ dropout provenance gate

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** proDA accepts qualified unimputed LFQ data and rejects TMT, unknown incompatible imputation and unsupported repeated models.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V056: proDA fitting and test_diff

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** zero-null/reduced-model reference calls match native output and retain model-specific uncertainty and dropout diagnostics.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V057: Alternative-engine uncertainty

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** engine-appropriate intervals/statistics are verified against native/reference calculations; no borrowed limma df or fake TREAT output.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V058: Method sensitivity comparison

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** matched-universe effects/discoveries are displayed as sensitivities, with unmatched universes explained; no hit-count winner changes primary.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V059: Dependencies and failure propagation

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** missing package/backend error is visible, required analysis fails, optional applicability is distinguished from installation failure.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V060: Assay-engine golden/calibration fixtures

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** count-dependent variance and LFQ dropout simulations have pinned reference comparisons, unsuitable-input negative cases and recorded diagnostics.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

## Edge Cases and Success Criteria

Test the named invalid/missing/empty/reordered/boundary cases without weakening the shared scientific contract. Every requirement must have a passing eligible-case test and a relevant explicit failure or boundary case where applicable. Preserve finite/NA distinctions and scientific eligibility reasons. Success requires real behavior, schema-valid outputs, independently rerun gates, reviewed diff and evidence linked to a commit; process exit alone is insufficient.

## Assumptions and Scope Boundary

Implement only the capabilities assigned above and the narrow helpers they require. Read the live constitution and parent contracts. An optional per-study analysis is still a mandatory implemented adapter when assigned here. Do not implement unrelated services, raw-MS analysis, private-data transmission or speculative interfaces. Use synthetic/public-permitted fixtures. Maintainer handles local private regression in the final slice.
