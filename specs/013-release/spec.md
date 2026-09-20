# Feature Specification: Archived-study regression, documentation and versioned successor

**Feature Branch:** `spec/013-release` (logical slice; stay on the current worktree branch unless a branch change is safe and needed).
**Created:** 2026-09-12. **Status:** Corrected draft awaiting Independent Reviewer audit and Maintainer freeze; not verified.
**Input:** Parent roadmap `specs/001-downstream-proteomics/roadmap.md` → **R12**. Deliver archived-study regression, documentation and versioned successor for user journey US6.

## User Scenarios & Testing

### US6 — Archived-study regression, documentation and versioned successor (Priority P1)

This slice delivers a usable and independently verifiable part of the parent user journey. It consumes the shared canonical contracts and exposes the behavior below through maintained package interfaces. It must not silently change the historical baseline or scientific interpretation.

**Why this priority:** dependent stages cannot reliably interpret their input without this contract being enforced.
**Independent Test:** run the acceptance cases below against actual implementations, including the negative cases; downstream slices may use controlled canonical fixtures rather than requiring an unfinished upstream feature.
**Dependencies:** R11 (`specs/012-validation/`).

## Requirements

- **FR-111**: The system MUST provide private local legacy adapter regression. Maintainer locally reconstructs archive values/masks/meta and confirms original source hashes without exposing study data to Packet Implementer, Independent Reviewer or any documentation worker.
- **FR-112**: The system MUST provide legacy core numerical comparison. Matched-universe legacy mode compares coefficients/BH/declared settings; missing historical environment is explicit NOT_REPRODUCED, not an invented success.
- **FR-113**: The system MUST provide scientific-change reconciliation. Negr1 auxiliary result, exact permutation discrepancy, different filters/families/pathway nulls and corrected response labels are explained without hardcoding outputs.
- **FR-114**: The system MUST provide full maintained study run. Qualified local study runs new pipeline, reports inherited unknown metadata and returns new independent artifacts without modifying originals.
- **FR-115**: The system MUST provide comprehensive usage and methods docs. Supported assays/engines/designs/unsupported combinations and exact command examples reflect implemented behavior and actual validation.
- **FR-116**: The system MUST provide spec-task-evidence reconciliation. Every requirement maps to a reviewed task/test/evidence/commit; unchecked or unverified work remains visible.
- **FR-117**: The system MUST provide source release hygiene. Artifact whitelist excludes private data, credentials, environments and execution logs; ownership/license is not invented and resource terms are recorded.
- **FR-118**: The system MUST provide versioned successor export. New releases/version folder receives verified source/spec/tests/docs/manifests; original archive package files and private archive hashes are unchanged.
- **FR-119**: The system MUST provide final integrated acceptance run. All runnable core/golden/integration/scientific gates rerun on the final reviewed tree, clean setup is demonstrated and residual external limits are explicit.
- **FR-120**: The system MUST provide review-and-handoff completion. Maintainer reviews final staged/unstaged/untracked diff, records verified commits and delivery paths; no unsupported completion claim or unauthorized release action occurs.

## Key Entities

Use entities/keys in `specs/001-downstream-proteomics/data-model.md`. All artifacts carry schema_version, run_id/plan_hash where applicable, source lineage and execution status. No alternate implicit identity or inference representation is permitted.

## Acceptance Scenarios

### V111: Private local legacy adapter regression

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** Maintainer locally reconstructs archive values/masks/meta and confirms original source hashes without exposing study data to Packet Implementer, Independent Reviewer or any documentation worker.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V112: Legacy core numerical comparison

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** matched-universe legacy mode compares coefficients/BH/declared settings; missing historical environment is explicit NOT_REPRODUCED, not an invented success.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V113: Scientific-change reconciliation

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** negr1 auxiliary result, exact permutation discrepancy, different filters/families/pathway nulls and corrected response labels are explained without hardcoding outputs.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V114: Full maintained study run

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** qualified local study runs new pipeline, reports inherited unknown metadata and returns new independent artifacts without modifying originals.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V115: Comprehensive usage and methods docs

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** supported assays/engines/designs/unsupported combinations and exact command examples reflect implemented behavior and actual validation.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V116: Spec-task-evidence reconciliation

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** every requirement maps to a reviewed task/test/evidence/commit; unchecked or unverified work remains visible.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V117: Source release hygiene

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** artifact whitelist excludes private data, credentials, environments and execution logs; ownership/license is not invented and resource terms are recorded.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V118: Versioned successor export

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** new releases/version folder receives verified source/spec/tests/docs/manifests; original archive package files and private archive hashes are unchanged.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V119: Final integrated acceptance run

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** all runnable core/golden/integration/scientific gates rerun on the final reviewed tree, clean setup is demonstrated and residual external limits are explicit.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V120: Review-and-handoff completion

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** the maintainer reviews the final staged/unstaged/untracked diff, records verified commits and delivery paths; no unsupported completion claim or unauthorized release action occurs.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

## Edge Cases and Success Criteria

Test the named invalid/missing/empty/reordered/boundary cases without weakening the shared scientific contract. Every requirement must have a passing eligible-case test and a relevant explicit failure or boundary case where applicable. Preserve finite/NA distinctions and scientific eligibility reasons. Success requires real behavior, schema-valid outputs, independently rerun gates, reviewed diff and evidence linked to a commit; process exit alone is insufficient.

## Assumptions and Scope Boundary

Implement only the capabilities assigned above and the narrow helpers they require. Read the live constitution and parent contracts. An optional per-study analysis is still a mandatory implemented adapter when assigned here. Do not implement unrelated services, raw-MS analysis, private-data transmission or speculative interfaces. Use synthetic/public-permitted fixtures. Maintainer handles local private regression in the final slice.
