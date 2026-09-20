# Tasks: Runtime, package and contract foundation

**Input:** [spec.md](spec.md), [plan.md](plan.md), parent shared contracts.
**Roadmap:** R01. **Dependencies:** Recovered baseline and approved parent contracts only.

## Phase A: Preconditions

Read the live constitution, frozen packet contract, prerequisite evidence/commits and existing interfaces. Maintainer issues one packet starter only after audit, freeze and explicit authorization; parallel work is limited to a frozen disjoint-file group.

## Phase B: Independently testable implementation

- [ ] T001 [US1] Implement and independently test **Installable maintained CLI** in `src/proteomics_pipeline/cli.py; pyproject.toml`. Acceptance V001: An isolated install exposes proteomics and python -m entry points; --help/version never require R or data. Evidence: command/output/version/hash in `docs/validation/002-foundation/V001.json`. Covers FR-001.
- [ ] T002 [US1] Implement and independently test **Strict versioned configuration** in `src/proteomics_pipeline/config.py; schemas/analysis.schema.json`. Acceptance V002: Valid examples pass; unknown keys, malformed contrasts and conflicting modes fail with field paths. Evidence: command/output/version/hash in `docs/validation/002-foundation/V002.json`. Covers FR-002.
- [ ] T003 [US1] Implement and independently test **Typed run and stage states** in `src/proteomics_pipeline/runtime.py; schemas/stage-result.schema.json`. Acceptance V003: A simulated failure cannot become COMPLETED; optional inapplicability has an eligibility reason. Evidence: command/output/version/hash in `docs/validation/002-foundation/V003.json`. Covers FR-003.
- [ ] T004 [US1] Implement and independently test **Safe Python-to-R bridge** in `scripts/maintained/run_stage.R; src/proteomics_pipeline/runtime.py`. Acceptance V004: Arguments with spaces are passed literally, nonzero R exit fails the stage, and stdout/stderr/session metadata survive. Evidence: command/output/version/hash in `docs/validation/002-foundation/V004.json`. Covers FR-004.
- [ ] T005 [US1] Implement and independently test **Atomic artifact publication** in `src/proteomics_pipeline/provenance.py`. Acceptance V005: An interrupted writer leaves no apparently completed artifact and a second writer cannot corrupt the run. Evidence: command/output/version/hash in `docs/validation/002-foundation/V005.json`. Covers FR-005.
- [ ] T006 [US1] Implement and independently test **Environment doctor** in `src/proteomics_pipeline/doctor.py`. Acceptance V006: Doctor lists actual binaries/package/resource versions; missing R/package is NOT_AVAILABLE and exits nonzero when required. Evidence: command/output/version/hash in `docs/validation/002-foundation/V006.json`. Covers FR-006.
- [ ] T007 [US1] Implement and independently test **Maintained R package skeleton** in `r/proteomicsCore/DESCRIPTION; r/proteomicsCore/NAMESPACE; scripts/maintained/test_r.R`. Acceptance V007: A real minimal R package installs and its executable test entry runs a nonmocked I/O roundtrip; it returns no fabricated analysis. Evidence: command/output/version/hash in `docs/validation/002-foundation/V007.json`. Covers FR-007.
- [ ] T008 [US1] Implement and independently test **Independent test harness** in `tests/contract/; tests/unit/; r/proteomicsCore/tests/testthat/`. Acceptance V008: Existing ten tests remain intact; pytest/testthat invoke real future test locations and missing prerequisites are not recorded as passes. Evidence: command/output/version/hash in `docs/validation/002-foundation/V008.json`. Covers FR-008.
- [ ] T009 [US1] Implement and independently test **Baseline preservation registry** in `docs/validation/baseline-hashes.json; scripts/maintained/check_baseline.py`. Acceptance V009: All recovered pipeline/legacy/audit hashes match and tests fail on a deliberate scratch-copy mutation. Evidence: command/output/version/hash in `docs/validation/002-foundation/V009.json`. Covers FR-009.
- [ ] T010 [US1] Implement and independently test **One-command developer bootstrap** in `scripts/maintained/bootstrap.ps1; scripts/maintained/bootstrap.sh; docs/development.md`. Acceptance V010: Fresh documented setup provisions only project-local environments, records solved versions and runs baseline plus foundation gates. Evidence: command/output/version/hash in `docs/validation/002-foundation/V010.json`. Covers FR-010.

## Phase C: Review checkpoint

All checkboxes start unchecked. Complete a checkbox only when Maintainer has read the diff, rerun the relevant tests and recorded acceptance evidence plus commit. A task can span multiple Packet Implementer correction tasks. If a runtime/resource is unavailable, keep its verification pending and continue independent work; do not write a successful fixture-shaped substitute.

## Dependencies & Execution Order

The order above follows shared contract → implementation → integration → adversarial verification. Preserve data-model and method boundaries. Resolve packet prerequisites before dependent tasks and obey the exact allowlist in the frozen packet index. A thin vertical demonstration is a checkpoint, not the end of the parent roadmap.
