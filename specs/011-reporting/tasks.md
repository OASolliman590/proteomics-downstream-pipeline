# Tasks: Unified offline report and streamlined workflow

**Input:** [spec.md](spec.md), [plan.md](plan.md), parent shared contracts.
**Roadmap:** R10. **Dependencies:** R03, R06, R08, R09.

## Phase A: Preconditions

Read the live constitution, frozen packet contract, prerequisite evidence/commits and existing interfaces. Maintainer issues one packet starter only after audit, freeze and explicit authorization; parallel work is limited to a frozen disjoint-file group.

## Phase B: Independently testable implementation

- [ ] T091 [US5] Implement and independently test **Typed report data assembler** in `src/proteomics_pipeline/reporting/data.py; schemas/report.schema.json`. Acceptance V091: All sections derive from validated artifact statuses/counts and missing inputs never become healthy zero or PASS. Evidence: command/output/version/hash in `docs/validation/011-reporting/V091.json`. Covers FR-091.
- [ ] T092 [US5] Implement and independently test **One-run command integration** in `src/proteomics_pipeline/cli.py; src/proteomics_pipeline/runtime.py`. Acceptance V092: Single run executes the valid DAG, stops failed required stages and assembles a useful partial report with failed stage evidence. Evidence: command/output/version/hash in `docs/validation/011-reporting/V092.json`. Covers FR-092.
- [ ] T093 [US5] Implement and independently test **Offline HTML report** in `src/proteomics_pipeline/reporting/templates/; src/proteomics_pipeline/reporting/render.py`. Acceptance V093: Report opens without external assets/network, has clear navigation and accessible tables, and no paths expose private content in source release. Evidence: command/output/version/hash in `docs/validation/011-reporting/V093.json`. Covers FR-093.
- [ ] T094 [US5] Implement and independently test **QC and inclusion/exclusion sections** in `src/proteomics_pipeline/reporting/qc_data.py`. Acceptance V094: Observation hierarchy, n, scale, missingness, normalization and exclusion rationale are visible with source-table links. Evidence: command/output/version/hash in `docs/validation/011-reporting/V094.json`. Covers FR-094.
- [ ] T095 [US5] Implement and independently test **Complete differential results views** in `src/proteomics_pipeline/reporting/model_data.py`. Acceptance V095: All planned contrasts/hypotheses appear, including nonsignificant/auxiliary contrasts; estimates, CI, q family and estimability are clear. Evidence: command/output/version/hash in `docs/validation/011-reporting/V095.json`. Covers FR-095.
- [ ] T096 [US5] Implement and independently test **Pathway and response views** in `src/proteomics_pipeline/reporting/pathway_data.py; response_data.py`. Acceptance V096: Null type, resource universe, leading edges, exploratory status, response caveats and unavailable inference remain explicit. Evidence: command/output/version/hash in `docs/validation/011-reporting/V096.json`. Covers FR-096.
- [ ] T097 [US5] Implement and independently test **Publication figure and source exports** in `r/proteomicsCore/R/plots.R`. Acceptance V097: Vector PDF/SVG and raster PNG as appropriate have exact plotted data, units/n/thresholds and readable labels without unexplained stars. Evidence: command/output/version/hash in `docs/validation/011-reporting/V097.json`. Covers FR-097.
- [ ] T098 [US5] Implement and independently test **Methods and limitations generation** in `src/proteomics_pipeline/reporting/methods.py`. Acceptance V098: Methods derive actual config/packages/hypotheses; unknown tissue/batch/provenance restricts language and no canned study claim appears. Evidence: command/output/version/hash in `docs/validation/011-reporting/V098.json`. Covers FR-098.
- [ ] T099 [US5] Implement and independently test **Null/empty/partial report robustness** in `tests/integration/test_reports.py`. Acceptance V099: No DEPs, no pathways, tiny cohorts, inapplicable methods, failed backend and interrupted stage yield valid honest reports. Evidence: command/output/version/hash in `docs/validation/011-reporting/V099.json`. Covers FR-099.
- [ ] T100 [US5] Implement and independently test **Quickstart and user workflow** in `docs/user-guide/; configs/examples/`. Acceptance V100: A fresh user can validate, plan, run, inspect, verify and compare examples without manually invoking individual R scripts; commands are actually exercised. Evidence: command/output/version/hash in `docs/validation/011-reporting/V100.json`. Covers FR-100.

## Phase C: Review checkpoint

All checkboxes start unchecked. Complete a checkbox only when Maintainer has read the diff, rerun the relevant tests and recorded acceptance evidence plus commit. A task can span multiple Packet Implementer correction tasks. If a runtime/resource is unavailable, keep its verification pending and continue independent work; do not write a successful fixture-shaped substitute.

## Dependencies & Execution Order

The order above follows shared contract → implementation → integration → adversarial verification. Preserve data-model and method boundaries. Resolve packet prerequisites before dependent tasks and obey the exact allowlist in the frozen packet index. A thin vertical demonstration is a checkpoint, not the end of the parent roadmap.
