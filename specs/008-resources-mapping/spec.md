# Feature Specification: Versioned annotation, protein groups and gene sets

**Feature Branch:** `spec/008-resources-mapping` (logical slice; stay on the current worktree branch unless a branch change is safe and needed).
**Created:** 2026-09-12. **Status:** Corrected draft awaiting Independent Reviewer audit and Maintainer freeze; not verified.
**Input:** Parent roadmap `specs/001-downstream-proteomics/roadmap.md` → **R07**. Deliver versioned annotation, protein groups and gene sets for user journey US3.

## User Scenarios & Testing

### US3 — Versioned annotation, protein groups and gene sets (Priority P1)

This slice delivers a usable and independently verifiable part of the parent user journey. It consumes the shared canonical contracts and exposes the behavior below through maintained package interfaces. It must not silently change the historical baseline or scientific interpretation.

**Why this priority:** dependent stages cannot reliably interpret their input without this contract being enforced.
**Independent Test:** run the acceptance cases below against actual implementations, including the negative cases; downstream slices may use controlled canonical fixtures rather than requiring an unfinished upstream feature.
**Dependencies:** R05 (`specs/006-limma-inference/`).

## Requirements

- **FR-061**: The system MUST provide local resource snapshot registry. Snapshots identify source/release/species/terms/hash; changed content or missing required files fails before production analysis.
- **FR-062**: The system MUST provide explicit resource preparation command. Fetch/build occurs only in a separate explicit preparation action, records exact versions and never auto-fetches during a frozen run.
- **FR-063**: The system MUST provide species-aware identifier mapping. Stable IDs/species keys align; symbol collisions, retired IDs and missing mappings have explicit records.
- **FR-064**: The system MUST provide orthology evidence preservation. Source/target species and mapping evidence/version are retained; projected sets are labeled projected and ambiguous mappings are reported.
- **FR-065**: The system MUST provide label-independent gene representatives. Coverage/median/ID representative choice is invariant to sample-label permutation and never chooses minimum P or maximum statistic.
- **FR-066**: The system MUST provide ambiguous protein-group policy. Multi-gene groups are excluded from default mapping with reasons; declared aggregation sensitivity retains its own universe and estimand.
- **FR-067**: The system MUST provide finite pathway matrix and model. Gene matrix is finite for required design, models are refit consistently and missing-value loss is quantified rather than silently imputed.
- **FR-068**: The system MUST provide gene-set overlap eligibility. Declared overlap-size filters are applied before tests, full/eligible/mapped membership counts are retained and tiny universe is explicit.
- **FR-069**: The system MUST provide measured gene background for ora. Universe equals eligible tested mapped genes and changes predictably under filtering; whole-genome background is never default.
- **FR-070**: The system MUST provide resource/mapping regression fixtures. Tests cover one-to-many mapping, duplicated symbols, zero overlap, ortholog projection, offline checksum errors and stable representative selection.

## Key Entities

Use entities/keys in `specs/001-downstream-proteomics/data-model.md`. All artifacts carry schema_version, run_id/plan_hash where applicable, source lineage and execution status. No alternate implicit identity or inference representation is permitted.

## Acceptance Scenarios

### V061: Local resource snapshot registry

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** snapshots identify source/release/species/terms/hash; changed content or missing required files fails before production analysis.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V062: Explicit resource preparation command

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** fetch/build occurs only in a separate explicit preparation action, records exact versions and never auto-fetches during a frozen run.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V063: Species-aware identifier mapping

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** stable IDs/species keys align; symbol collisions, retired IDs and missing mappings have explicit records.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V064: Orthology evidence preservation

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** source/target species and mapping evidence/version are retained; projected sets are labeled projected and ambiguous mappings are reported.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V065: Label-independent gene representatives

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** coverage/median/ID representative choice is invariant to sample-label permutation and never chooses minimum P or maximum statistic.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V066: Ambiguous protein-group policy

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** multi-gene groups are excluded from default mapping with reasons; declared aggregation sensitivity retains its own universe and estimand.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V067: Finite pathway matrix and model

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** gene matrix is finite for required design, models are refit consistently and missing-value loss is quantified rather than silently imputed.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V068: Gene-set overlap eligibility

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** declared overlap-size filters are applied before tests, full/eligible/mapped membership counts are retained and tiny universe is explicit.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V069: Measured gene background for ORA

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** universe equals eligible tested mapped genes and changes predictably under filtering; whole-genome background is never default.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V070: Resource/mapping regression fixtures

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** tests cover one-to-many mapping, duplicated symbols, zero overlap, ortholog projection, offline checksum errors and stable representative selection.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

## Edge Cases and Success Criteria

Test the named invalid/missing/empty/reordered/boundary cases without weakening the shared scientific contract. Every requirement must have a passing eligible-case test and a relevant explicit failure or boundary case where applicable. Preserve finite/NA distinctions and scientific eligibility reasons. Success requires real behavior, schema-valid outputs, independently rerun gates, reviewed diff and evidence linked to a commit; process exit alone is insufficient.

## Assumptions and Scope Boundary

Implement only the capabilities assigned above and the narrow helpers they require. Read the live constitution and parent contracts. An optional per-study analysis is still a mandatory implemented adapter when assigned here. Do not implement unrelated services, raw-MS analysis, private-data transmission or speculative interfaces. Use synthetic/public-permitted fixtures. Maintainer handles local private regression in the final slice.
