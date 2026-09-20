# Tasks: Archived-study regression, documentation and versioned successor

**Input:** [spec.md](spec.md), [plan.md](plan.md), parent shared contracts.
**Roadmap:** R12. **Dependencies:** R11.

## Phase A: Preconditions

Read the live constitution, frozen packet contract, prerequisite evidence/commits and existing interfaces. Maintainer issues one packet starter only after audit, freeze and explicit authorization; parallel work is limited to a frozen disjoint-file group.

## Phase B: Independently testable implementation

- [ ] T111 [US6] Implement and independently test **Private local legacy adapter regression** in `scripts/maintained/legacy_regression.py`. Acceptance V111: Maintainer locally reconstructs archive values/masks/meta and confirms original source hashes without exposing study data to Packet Implementer, Independent Reviewer or any documentation worker. Evidence: command/output/version/hash in `docs/validation/013-release/V111.json`. Covers FR-111.
- [ ] T112 [US6] Implement and independently test **Legacy core numerical comparison** in `tests/regression/; docs/validation/legacy-comparison.md`. Acceptance V112: Matched-universe legacy mode compares coefficients/BH/declared settings; missing historical environment is explicit NOT_REPRODUCED, not an invented success. Evidence: command/output/version/hash in `docs/validation/013-release/V112.json`. Covers FR-112.
- [ ] T113 [US6] Implement and independently test **Scientific-change reconciliation** in `docs/validation/legacy-comparison.md`. Acceptance V113: Negr1 auxiliary result, exact permutation discrepancy, different filters/families/pathway nulls and corrected response labels are explained without hardcoding outputs. Evidence: command/output/version/hash in `docs/validation/013-release/V113.json`. Covers FR-113.
- [ ] T114 [US6] Implement and independently test **Full maintained study run** in `configs/local example template; scripts/maintained/legacy_regression.py`. Acceptance V114: Qualified local study runs new pipeline, reports inherited unknown metadata and returns new independent artifacts without modifying originals. Evidence: command/output/version/hash in `docs/validation/013-release/V114.json`. Covers FR-114.
- [ ] T115 [US6] Implement and independently test **Comprehensive usage and methods docs** in `docs/user-guide/; docs/methods/; README.md`. Acceptance V115: Supported assays/engines/designs/unsupported combinations and exact command examples reflect implemented behavior and actual validation. Evidence: command/output/version/hash in `docs/validation/013-release/V115.json`. Covers FR-115.
- [ ] T116 [US6] Implement and independently test **Spec-task-evidence reconciliation** in `specs/001-downstream-proteomics/traceability.json; progress.md`. Acceptance V116: Every requirement maps to a reviewed task/test/evidence/commit; unchecked or unverified work remains visible. Evidence: command/output/version/hash in `docs/validation/013-release/V116.json`. Covers FR-116.
- [ ] T117 [US6] Implement and independently test **Source release hygiene** in `scripts/maintained/build_release.py; .gitignore`. Acceptance V117: Artifact whitelist excludes private data, credentials, environments and execution logs; ownership/license is not invented and resource terms are recorded. Evidence: command/output/version/hash in `docs/validation/013-release/V117.json`. Covers FR-117.
- [ ] T118 [US6] Implement and independently test **Versioned successor export** in `scripts/maintained/build_release.py; docs/validation/delivery.json`. Acceptance V118: New releases/version folder receives verified source/spec/tests/docs/manifests; original archive package files and private archive hashes are unchanged. Evidence: command/output/version/hash in `docs/validation/013-release/V118.json`. Covers FR-118.
- [ ] T119 [US6] Implement and independently test **Final integrated acceptance run** in `docs/validation/final-acceptance.md`. Acceptance V119: All runnable core/golden/integration/scientific gates rerun on the final reviewed tree, clean setup is demonstrated and residual external limits are explicit. Evidence: command/output/version/hash in `docs/validation/013-release/V119.json`. Covers FR-119.
- [ ] T120 [US6] Implement and independently test **Review-and-handoff completion** in `CHANGELOG.md; specs/001-downstream-proteomics/progress.md`. Acceptance V120: Maintainer reviews final staged/unstaged/untracked diff, records verified commits and delivery paths; no unsupported completion claim or unauthorized release action occurs. Evidence: command/output/version/hash in `docs/validation/013-release/V120.json`. Covers FR-120.

## Phase C: Review checkpoint

All checkboxes start unchecked. Complete a checkbox only when Maintainer has read the diff, rerun the relevant tests and recorded acceptance evidence plus commit. A task can span multiple Packet Implementer correction tasks. If a runtime/resource is unavailable, keep its verification pending and continue independent work; do not write a successful fixture-shaped substitute.

## Dependencies & Execution Order

The order above follows shared contract → implementation → integration → adversarial verification. Preserve data-model and method boundaries. Resolve packet prerequisites before dependent tasks and obey the exact allowlist in the frozen packet index. A thin vertical demonstration is a checkpoint, not the end of the parent roadmap.
