# Feature Specification: Canonical protein input and biological identity

**Feature Branch:** `spec/003-intake` (logical slice; stay on the current worktree branch unless a branch change is safe and needed).
**Created:** 2026-09-12. **Status:** Corrected draft awaiting Independent Reviewer audit and Maintainer freeze; not verified.
**Input:** Parent roadmap `specs/001-downstream-proteomics/roadmap.md` → **R02**. Deliver canonical protein input and biological identity for user journey US1.

## User Scenarios & Testing

### US1 — Canonical protein input and biological identity (Priority P1)

This slice delivers a usable and independently verifiable part of the parent user journey. It consumes the shared canonical contracts and exposes the behavior below through maintained package interfaces. It must not silently change the historical baseline or scientific interpretation.

**Why this priority:** dependent stages cannot reliably interpret their input without this contract being enforced.
**Independent Test:** run the acceptance cases below against actual implementations, including the negative cases; downstream slices may use controlled canonical fixtures rather than requiring an unfinished upstream feature.
**Dependencies:** R01 (`specs/002-foundation/`).

## Requirements

- **FR-011**: The system MUST provide canonical wide and long intake. Equivalent shuffled wide/long tables yield identical aligned numeric values and masks; duplicate feature-observation pairs fail.
- **FR-012**: The system MUST provide explicit scale and missing encoding. Linear-positive input is transformed once; log2 zero is retained; configured linear zero-as-missing and invalid negative values are distinguished.
- **FR-013**: The system MUST provide biological observation hierarchy. Technical injections never inflate biological n; paired and repeated observations retain subject identity; invalid shared IDs fail.
- **FR-014**: The system MUST provide explicit technical aggregation. Linear-mean and log-mean produce hand-calculated different results; coverage and specimen counts are preserved.
- **FR-015**: The system MUST provide feature annotations and protein groups. Opaque row IDs remain unique while repeated accessions/gene mappings remain explicit; contaminant/decoy flags distinguish unknown from false.
- **FR-016**: The system MUST provide versioned vendor mappings. Synthetic DIA-NN, MaxQuant, FragPipe protein and Spectronaut protein fixtures exercise declared mapping/version; unsupported peptide grain fails.
- **FR-017**: The system MUST provide historical workbook importer. The existing parser is wrapped without editing its source; source-sheet contrast and Method A/B naming are mapped by biological meaning.
- **FR-018**: The system MUST provide immutable canonical bundle. Matrix, masks, annotations, observation metadata and hashes roundtrip without numerical change; reruns refuse collisions.
- **FR-019**: The system MUST provide intake validation report. Reports include grain, scale, n, missingness, duplicate keys, provenance gaps and errors; no hardcoded biological confirmation.
- **FR-020**: The system MUST provide input edge-case fixtures. Exercise ragged files, duplicate IDs, mismatched columns, empty values, Inf, quoted delimiters, unknown encoding and unsupported assay scope.

## Key Entities

Use entities/keys in `specs/001-downstream-proteomics/data-model.md`. All artifacts carry schema_version, run_id/plan_hash where applicable, source lineage and execution status. No alternate implicit identity or inference representation is permitted.

## Acceptance Scenarios

### V011: Canonical wide and long intake

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** equivalent shuffled wide/long tables yield identical aligned numeric values and masks; duplicate feature-observation pairs fail.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V012: Explicit scale and missing encoding

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** linear-positive input is transformed once; log2 zero is retained; configured linear zero-as-missing and invalid negative values are distinguished.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V013: Biological observation hierarchy

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** technical injections never inflate biological n; paired and repeated observations retain subject identity; invalid shared IDs fail.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V014: Explicit technical aggregation

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** linear-mean and log-mean produce hand-calculated different results; coverage and specimen counts are preserved.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V015: Feature annotations and protein groups

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** opaque row IDs remain unique while repeated accessions/gene mappings remain explicit; contaminant/decoy flags distinguish unknown from false.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V016: Versioned vendor mappings

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** synthetic DIA-NN, MaxQuant, FragPipe protein and Spectronaut protein fixtures exercise declared mapping/version; unsupported peptide grain fails.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V017: Historical workbook importer

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** the existing parser is wrapped without editing its source; source-sheet contrast and Method A/B naming are mapped by biological meaning.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V018: Immutable canonical bundle

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** matrix, masks, annotations, observation metadata and hashes roundtrip without numerical change; reruns refuse collisions.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V019: Intake validation report

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** reports include grain, scale, n, missingness, duplicate keys, provenance gaps and errors; no hardcoded biological confirmation.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

### V020: Input edge-case fixtures

**Given** a fixture exercising the declared scientific/input conditions, **when** this capability executes, **then** exercise ragged files, duplicate IDs, mismatched columns, empty values, Inf, quoted delimiters, unknown encoding and unsupported assay scope.

Use an analytic calculation, independently invoked package reference, or observable failure/invariance behavior. A mock of the function under test is not evidence. Record the expected value/state before inspecting candidate output.

## Edge Cases and Success Criteria

Test the named invalid/missing/empty/reordered/boundary cases without weakening the shared scientific contract. Every requirement must have a passing eligible-case test and a relevant explicit failure or boundary case where applicable. Preserve finite/NA distinctions and scientific eligibility reasons. Success requires real behavior, schema-valid outputs, independently rerun gates, reviewed diff and evidence linked to a commit; process exit alone is insufficient.

## Assumptions and Scope Boundary

Implement only the capabilities assigned above and the narrow helpers they require. Read the live constitution and parent contracts. An optional per-study analysis is still a mandatory implemented adapter when assigned here. Do not implement unrelated services, raw-MS analysis, private-data transmission or speculative interfaces. Use synthetic/public-permitted fixtures. Maintainer handles local private regression in the final slice.
