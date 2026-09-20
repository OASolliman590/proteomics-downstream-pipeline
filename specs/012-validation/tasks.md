# Tasks: Reproducibility, calibration and continuous validation

**Input:** [spec.md](spec.md), [plan.md](plan.md), parent shared contracts.
**Roadmap:** R11. **Dependencies:** R10.

## Phase A: Preconditions

Read the live constitution, frozen packet contract, prerequisite evidence/commits and existing interfaces. Maintainer issues one packet starter only after audit, freeze and explicit authorization; parallel work is limited to a frozen disjoint-file group.

## Phase B: Independently testable implementation

- [ ] T101 [US6] Implement and independently test **Pinned compatible environments** in `renv.lock; Python lockfile; containers/`. Acceptance V101: Clean Python/R/Bioconductor restore is tested; version pair is compatible and installed package APIs match adapters. Evidence: command/output/version/hash in `docs/validation/012-validation/V101.json`. Covers FR-101.
- [ ] T102 [US6] Implement and independently test **Immutable resource-offline reproduction** in `src/proteomics_pipeline/provenance.py; tests/integration/test_offline.py`. Acceptance V102: A run with cached snapshots completes without network and semantic output hashes match the documented tolerance. Evidence: command/output/version/hash in `docs/validation/012-validation/V102.json`. Covers FR-102.
- [ ] T103 [US6] Implement and independently test **Resume and cache invalidation** in `src/proteomics_pipeline/runtime.py; tests/integration/test_resume.py`. Acceptance V103: Code/config/data/resource changes invalidate dependent stages; compatible verified outputs reuse safely and partial artifacts never count as cache. Evidence: command/output/version/hash in `docs/validation/012-validation/V103.json`. Covers FR-103.
- [ ] T104 [US6] Implement and independently test **Scientific golden reference matrix** in `tests/scientific/reference/; docs/validation/reference-matrix.json`. Acceptance V104: Independent calls cover every backend/hypothesis/design eligibility; tolerances/versions are frozen before candidate comparison. Evidence: command/output/version/hash in `docs/validation/012-validation/V104.json`. Covers FR-104.
- [ ] T105 [US6] Implement and independently test **Null and mixture calibration** in `scripts/maintained/run_calibration.py; r/proteomicsCore/R/simulate.R`. Acceptance V105: Prespecified all-null/mixture scenarios calculate actual error/coverage metrics with Monte Carlo uncertainty; failing methods are not labeled validated. Evidence: command/output/version/hash in `docs/validation/012-validation/V105.json`. Covers FR-105.
- [ ] T106 [US6] Implement and independently test **Correlated pathway calibration** in `tests/scientific/calibration_pathways.R`. Acceptance V106: Known correlated null and spike-in pathways assess method-null-specific error/power; fgsea caveat is retained rather than demanding a false universal guarantee. Evidence: command/output/version/hash in `docs/validation/012-validation/V106.json`. Covers FR-106.
- [ ] T107 [US6] Implement and independently test **Cross-platform CI** in `.github/workflows/maintained.yml; scripts/maintained/`. Acceptance V107: Windows/Linux core Python+R jobs execute real examples; expensive release calibration is distinct and not silently skipped as passing CI. Evidence: command/output/version/hash in `docs/validation/012-validation/V107.json`. Covers FR-107.
- [ ] T108 [US6] Implement and independently test **Performance and resource benchmark** in `scripts/maintained/benchmark.py; docs/validation/benchmarks/`. Acceptance V108: 20k x100 medium case logs time/peak memory/hardware, specialized engine cost separately; targets are assessed on actual work. Evidence: command/output/version/hash in `docs/validation/012-validation/V108.json`. Covers FR-108.
- [ ] T109 [US6] Implement and independently test **Report numerical/visual validation** in `tests/integration/test_report_artifacts.py; docs/validation/`. Acceptance V109: Plot-data values and references match tables; representative rendered reports/figures are visually inspected with explicit evidence. Evidence: command/output/version/hash in `docs/validation/012-validation/V109.json`. Covers FR-109.
- [ ] T110 [US6] Implement and independently test **Machine-readable validation ledger** in `schemas/validation-evidence.schema.json; docs/validation/acceptance.json`. Acceptance V110: Every requirement/test records command, versions, artifact hash, result and limitations; NOT_RUN/SKIPPED is not PASS. Evidence: command/output/version/hash in `docs/validation/012-validation/V110.json`. Covers FR-110.

## Phase C: Review checkpoint

All checkboxes start unchecked. Complete a checkbox only when Maintainer has read the diff, rerun the relevant tests and recorded acceptance evidence plus commit. A task can span multiple Packet Implementer correction tasks. If a runtime/resource is unavailable, keep its verification pending and continue independent work; do not write a successful fixture-shaped substitute.

## Dependencies & Execution Order

The order above follows shared contract → implementation → integration → adversarial verification. Preserve data-model and method boundaries. Resolve packet prerequisites before dependent tasks and obey the exact allowlist in the frozen packet index. A thin vertical demonstration is a checkpoint, not the end of the parent roadmap.
