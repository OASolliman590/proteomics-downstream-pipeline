# Feature Specification: Runtime, package and contract foundation

**Feature Branch:** `spec/002-foundation` (logical slice; stay on the current worktree branch unless a branch change is safe and needed).
**Created:** 2026-09-12. **Status:** Corrected draft awaiting Independent Reviewer audit and Maintainer freeze; not verified.
**Input:** Parent roadmap `specs/001-downstream-proteomics/roadmap.md` → **R01**. Deliver runtime, package and contract foundation for user journey US1.

## User Scenarios & Testing

### US1 — Runtime, package and contract foundation (Priority P1)

This slice delivers a usable and independently verifiable part of the parent user journey. It consumes the shared canonical contracts and exposes the behavior below through maintained package interfaces. It must not silently change the historical baseline or scientific interpretation.

**Why this priority:** dependent stages cannot reliably interpret their input without this contract being enforced.
**Independent Test:** run the acceptance cases below against actual implementations, including the negative cases; downstream slices may use controlled canonical fixtures rather than requiring an unfinished upstream feature.
**Dependencies:** Recovered baseline and approved parent contracts only.

## Requirements

- **FR-001**: The system MUST provide installable maintained cli. An isolated install exposes proteomics and python -m entry points; --help/version never require R or data.
- **FR-002**: The system MUST provide strict versioned configuration. Valid examples pass; unknown keys, malformed contrasts and conflicting modes fail with field paths.
- **FR-003**: The system MUST provide typed run and stage states. A simulated failure cannot become COMPLETED; optional inapplicability has an eligibility reason.
- **FR-004**: The system MUST provide safe python-to-r bridge. Arguments with spaces are passed literally, nonzero R exit fails the stage, and stdout/stderr/session metadata survive.
- **FR-005**: The system MUST provide atomic artifact publication. An interrupted writer leaves no apparently completed artifact and a second writer cannot corrupt the run.
- **FR-006**: The system MUST provide environment doctor. Doctor lists actual binaries/package/resource versions; missing R/package is NOT_AVAILABLE and exits nonzero when required.
- **FR-007**: The system MUST provide maintained r package skeleton. A real minimal R package installs and its executable test entry runs a nonmocked I/O roundtrip; it returns no fabricated analysis.
- **FR-008**: The system MUST provide independent test harness. Existing ten tests remain intact; pytest/testthat invoke real future test locations and missing prerequisites are not recorded as passes.
- **FR-009**: The system MUST provide baseline preservation registry. All recovered pipeline/legacy/audit hashes match and tests fail on a deliberate scratch-copy mutation.
- **FR-010**: The system MUST provide one-command developer bootstrap. Fresh documented setup provisions only project-local environments, records solved versions and runs baseline plus foundation gates.

## Key Entities

Use entities/keys in `specs/001-downstream-proteomics/data-model.md`. All artifacts carry schema_version, run_id/plan_hash where applicable, source lineage and execution status. No alternate implicit identity or inference representation is permitted.

## Acceptance Scenarios

### V001: Installable maintained CLI

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** an isolated install exposes proteomics and python -m entry points; --help/version never require R or data.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V002: Strict versioned configuration

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** valid examples pass; unknown keys, malformed contrasts and conflicting modes fail with field paths.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V003: Typed run and stage states

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** a simulated failure cannot become COMPLETED; optional inapplicability has an eligibility reason.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V004: Safe Python-to-R bridge

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** arguments with spaces are passed literally, nonzero R exit fails the stage, and stdout/stderr/session metadata survive.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V005: Atomic artifact publication

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** an interrupted writer leaves no apparently completed artifact and a second writer cannot corrupt the run.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V006: Environment doctor

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** doctor lists actual binaries/package/resource versions; missing R/package is NOT_AVAILABLE and exits nonzero when required.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V007: Maintained R package skeleton

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** a real minimal R package installs and its executable test entry runs a nonmocked I/O roundtrip; it returns no fabricated analysis.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V008: Independent test harness

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** existing ten tests remain intact; pytest/testthat invoke real future test locations and missing prerequisites are not recorded as passes.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V009: Baseline preservation registry

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** all recovered pipeline/legacy/audit hashes match and tests fail on a deliberate scratch-copy mutation.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V010: One-command developer bootstrap

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** fresh documented setup provisions only project-local environments, records solved versions and runs baseline plus foundation gates.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

## Edge Cases and Success Criteria

Test the named invalid/missing/empty/reordered/boundary cases without weakening the shared scientific contract. Every requirement must have a passing eligible-case test and a relevant explicit failure or boundary case where applicable. Preserve finite/NA distinctions and scientific eligibility reasons. Success requires real behavior, schema-valid outputs, independently rerun gates, reviewed diff and evidence linked to a commit; process exit alone is insufficient.

## Assumptions and Scope Boundary

Implement only the capabilities assigned above and the narrow helpers they require. Read the live constitution and parent contracts. An optional per-study analysis is still a mandatory implemented adapter when assigned here. Do not implement unrelated services, raw-MS analysis, private-data transmission or speculative interfaces. Use synthetic/public-permitted fixtures. Maintainer handles local private regression in the final slice.
