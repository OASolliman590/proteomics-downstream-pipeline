# Tasks: Canonical protein input and biological identity

**Input:** [spec.md](spec.md), [plan.md](plan.md), parent shared contracts.
**Roadmap:** R02. **Dependencies:** R01.

## Phase A: Preconditions

Read the live constitution, frozen packet contract, prerequisite evidence/commits and existing interfaces. Maintainer issues one packet starter only after audit, freeze and explicit authorization; parallel work is limited to a frozen disjoint-file group.

## Phase B: Independently testable implementation

- [ ] T011 [US1] Implement and independently test **Canonical wide and long intake** in `src/proteomics_pipeline/intake/canonical.py`. Acceptance V011: Equivalent shuffled wide/long tables yield identical aligned numeric values and masks; duplicate feature-observation pairs fail. Evidence: command/output/version/hash in `docs/validation/003-intake/V011.json`. Covers FR-011.
- [ ] T012 [US1] Implement and independently test **Explicit scale and missing encoding** in `src/proteomics_pipeline/intake/scale.py`. Acceptance V012: Linear-positive input is transformed once; log2 zero is retained; configured linear zero-as-missing and invalid negative values are distinguished. Evidence: command/output/version/hash in `docs/validation/003-intake/V012.json`. Covers FR-012.
- [ ] T013 [US1] Implement and independently test **Biological observation hierarchy** in `src/proteomics_pipeline/intake/samples.py`. Acceptance V013: Technical injections never inflate biological n; paired and repeated observations retain subject identity; invalid shared IDs fail. Evidence: command/output/version/hash in `docs/validation/003-intake/V013.json`. Covers FR-013.
- [ ] T014 [US1] Implement and independently test **Explicit technical aggregation** in `r/proteomicsCore/R/intake.R; tests/scientific/test_technical_units.py`. Acceptance V014: Linear-mean and log-mean produce hand-calculated different results; coverage and specimen counts are preserved. Evidence: command/output/version/hash in `docs/validation/003-intake/V014.json`. Covers FR-014.
- [ ] T015 [US1] Implement and independently test **Feature annotations and protein groups** in `src/proteomics_pipeline/intake/features.py`. Acceptance V015: Opaque row IDs remain unique while repeated accessions/gene mappings remain explicit; contaminant/decoy flags distinguish unknown from false. Evidence: command/output/version/hash in `docs/validation/003-intake/V015.json`. Covers FR-015.
- [ ] T016 [US1] Implement and independently test **Versioned vendor mappings** in `src/proteomics_pipeline/intake/vendor.py; configs/mappings/`. Acceptance V016: Synthetic DIA-NN, MaxQuant, FragPipe protein and Spectronaut protein fixtures exercise declared mapping/version; unsupported peptide grain fails. Evidence: command/output/version/hash in `docs/validation/003-intake/V016.json`. Covers FR-016.
- [ ] T017 [US1] Implement and independently test **Historical workbook importer** in `src/proteomics_pipeline/intake/legacy.py`. Acceptance V017: The existing parser is wrapped without editing its source; source-sheet contrast and Method A/B naming are mapped by biological meaning. Evidence: command/output/version/hash in `docs/validation/003-intake/V017.json`. Covers FR-017.
- [ ] T018 [US1] Implement and independently test **Immutable canonical bundle** in `src/proteomics_pipeline/intake/bundle.py`. Acceptance V018: Matrix, masks, annotations, observation metadata and hashes roundtrip without numerical change; reruns refuse collisions. Evidence: command/output/version/hash in `docs/validation/003-intake/V018.json`. Covers FR-018.
- [ ] T019 [US1] Implement and independently test **Intake validation report** in `src/proteomics_pipeline/intake/validation.py`. Acceptance V019: Reports include grain, scale, n, missingness, duplicate keys, provenance gaps and errors; no hardcoded biological confirmation. Evidence: command/output/version/hash in `docs/validation/003-intake/V019.json`. Covers FR-019.
- [ ] T020 [US1] Implement and independently test **Input edge-case fixtures** in `tests/contract/test_intake.py; tests/fixtures/synthetic/`. Acceptance V020: Exercise ragged files, duplicate IDs, mismatched columns, empty values, Inf, quoted delimiters, unknown encoding and unsupported assay scope. Evidence: command/output/version/hash in `docs/validation/003-intake/V020.json`. Covers FR-020.

## Phase C: Review checkpoint

All checkboxes start unchecked. Complete a checkbox only when Maintainer has read the diff, rerun the relevant tests and recorded acceptance evidence plus commit. A task can span multiple Packet Implementer correction tasks. If a runtime/resource is unavailable, keep its verification pending and continue independent work; do not write a successful fixture-shaped substitute.

## Dependencies & Execution Order

The order above follows shared contract → implementation → integration → adversarial verification. Preserve data-model and method boundaries. Resolve packet prerequisites before dependent tasks and obey the exact allowlist in the frozen packet index. A thin vertical demonstration is a checkpoint, not the end of the parent roadmap.
