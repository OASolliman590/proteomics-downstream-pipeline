# Feature Specification: Unified offline report and streamlined workflow

**Feature Branch:** `spec/011-reporting` (logical slice; stay on the current worktree branch unless a branch change is safe and needed).
**Created:** 2026-09-12. **Status:** Corrected draft awaiting Independent Reviewer audit and Maintainer freeze; not verified.
**Input:** Parent roadmap `specs/001-downstream-proteomics/roadmap.md` → **R10**. Deliver unified offline report and streamlined workflow for user journey US5.

## User Scenarios & Testing

### US5 — Unified offline report and streamlined workflow (Priority P1)

This slice delivers a usable and independently verifiable part of the parent user journey. It consumes the shared canonical contracts and exposes the behavior below through maintained package interfaces. It must not silently change the historical baseline or scientific interpretation.

**Why this priority:** dependent stages cannot reliably interpret their input without this contract being enforced.
**Independent Test:** run the acceptance cases below against actual implementations, including the negative cases; downstream slices may use controlled canonical fixtures rather than requiring an unfinished upstream feature.
**Dependencies:** R03 (`specs/004-preprocessing-qc/`), R06 (`specs/007-assay-engines/`), R08 (`specs/009-enrichment/`) and R09 (`specs/010-treatment-response/`).

## Requirements

- **FR-091**: The system MUST provide typed report data assembler. All sections derive from validated artifact statuses/counts and missing inputs never become healthy zero or PASS.
- **FR-092**: The system MUST provide one-run command integration. Single run executes the valid DAG, stops failed required stages and assembles a useful partial report with failed stage evidence.
- **FR-093**: The system MUST provide offline html report. Report opens without external assets/network, has clear navigation and accessible tables, and no paths expose private content in source release.
- **FR-094**: The system MUST provide qc and inclusion/exclusion sections. Observation hierarchy, n, scale, missingness, normalization and exclusion rationale are visible with source-table links.
- **FR-095**: The system MUST provide complete differential results views. All planned contrasts/hypotheses appear, including nonsignificant/auxiliary contrasts; estimates, CI, q family and estimability are clear.
- **FR-096**: The system MUST provide pathway and response views. Null type, resource universe, leading edges, exploratory status, response caveats and unavailable inference remain explicit.
- **FR-097**: The system MUST provide publication figure and source exports. Vector PDF/SVG and raster PNG as appropriate have exact plotted data, units/n/thresholds and readable labels without unexplained stars.
- **FR-098**: The system MUST provide methods and limitations generation. Methods derive actual config/packages/hypotheses; unknown tissue/batch/provenance restricts language and no canned study claim appears.
- **FR-099**: The system MUST provide null/empty/partial report robustness. No DEPs, no pathways, tiny cohorts, inapplicable methods, failed backend and interrupted stage yield valid honest reports.
- **FR-100**: The system MUST provide quickstart and user workflow. A fresh user can validate, plan, run, inspect, verify and compare examples without manually invoking individual R scripts; commands are actually exercised.

## Key Entities

Use entities/keys in `specs/001-downstream-proteomics/data-model.md`. All artifacts carry schema_version, run_id/plan_hash where applicable, source lineage and execution status. No alternate implicit identity or inference representation is permitted.

## Acceptance Scenarios

### V091: Typed report data assembler

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** all sections derive from validated artifact statuses/counts and missing inputs never become healthy zero or PASS.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V092: One-run command integration

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** single run executes the valid DAG, stops failed required stages and assembles a useful partial report with failed stage evidence.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V093: Offline HTML report

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** report opens without external assets/network, has clear navigation and accessible tables, and no paths expose private content in source release.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V094: QC and inclusion/exclusion sections

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** observation hierarchy, n, scale, missingness, normalization and exclusion rationale are visible with source-table links.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V095: Complete differential results views

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** all planned contrasts/hypotheses appear, including nonsignificant/auxiliary contrasts; estimates, CI, q family and estimability are clear.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V096: Pathway and response views

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** null type, resource universe, leading edges, exploratory status, response caveats and unavailable inference remain explicit.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V097: Publication figure and source exports

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** vector PDF/SVG and raster PNG as appropriate have exact plotted data, units/n/thresholds and readable labels without unexplained stars.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V098: Methods and limitations generation

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** methods derive actual config/packages/hypotheses; unknown tissue/batch/provenance restricts language and no canned study claim appears.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V099: Null/empty/partial report robustness

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** no DEPs, no pathways, tiny cohorts, inapplicable methods, failed backend and interrupted stage yield valid honest reports.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V100: Quickstart and user workflow

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** a fresh user can validate, plan, run, inspect, verify and compare examples without manually invoking individual R scripts; commands are actually exercised.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

## Edge Cases and Success Criteria

Test the named invalid/missing/empty/reordered/boundary cases without weakening the shared scientific contract. Every requirement must have a passing eligible-case test and a relevant explicit failure or boundary case where applicable. Preserve finite/NA distinctions and scientific eligibility reasons. Success requires real behavior, schema-valid outputs, independently rerun gates, reviewed diff and evidence linked to a commit; process exit alone is insufficient.

## Assumptions and Scope Boundary

Implement only the capabilities assigned above and the narrow helpers they require. Read the live constitution and parent contracts. An optional per-study analysis is still a mandatory implemented adapter when assigned here. Do not implement unrelated services, raw-MS analysis, private-data transmission or speculative interfaces. Use synthetic/public-permitted fixtures. Maintainer handles local private regression in the final slice.
