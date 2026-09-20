# Tasks: Assay-qualified DEqMS and proDA backends

**Input:** [spec.md](spec.md), [plan.md](plan.md), parent shared contracts.
**Roadmap:** R06. **Dependencies:** R05.

## Phase A: Preconditions

Read the live constitution, frozen packet contract, prerequisite evidence/commits and existing interfaces. Maintainer issues one packet starter only after audit, freeze and explicit authorization; parallel work is limited to a frozen disjoint-file group.

## Phase B: Independently testable implementation

- [ ] T051 [US2] Implement and independently test **Count evidence schema** in `schemas/count-evidence.schema.json; src/proteomics_pipeline/intake/counts.py`. Acceptance V051: True peptide/PSM count grain and provenance validate; observed-sample count or abundance proxy is rejected as unsupported evidence. Evidence: command/output/version/hash in `docs/validation/007-assay-engines/V051.json`. Covers FR-051.
- [ ] T052 [US2] Implement and independently test **DEqMS fit pipeline** in `r/proteomicsCore/R/model_deqms.R`. Acceptance V052: Installed official fit/eBayes/count/spectraCounteBayes sequence reproduces reference count-adjusted statistics. Evidence: command/output/version/hash in `docs/validation/007-assay-engines/V052.json`. Covers FR-052.
- [ ] T053 [US2] Implement and independently test **DEqMS statistic alignment** in `r/proteomicsCore/R/model_deqms.R`. Acceptance V053: sca P/statistic/variance fields are kept separate from original limma fields; shuffled count IDs align or fail, never silently reorder. Evidence: command/output/version/hash in `docs/validation/007-assay-engines/V053.json`. Covers FR-053.
- [ ] T054 [US2] Implement and independently test **DEqMS hypotheses and design guard** in `r/proteomicsCore/R/capabilities.R`. Acceptance V054: Unsupported block/precision-weight/threshold-null requests are rejected with explicit contract; eligible sparse exact-contrast use is verified. Evidence: command/output/version/hash in `docs/validation/007-assay-engines/V054.json`. Covers FR-054.
- [ ] T055 [US2] Implement and independently test **LFQ dropout provenance gate** in `src/proteomics_pipeline/planning.py; r/proteomicsCore/R/model_proda.R`. Acceptance V055: proDA accepts qualified unimputed LFQ data and rejects TMT, unknown incompatible imputation and unsupported repeated models. Evidence: command/output/version/hash in `docs/validation/007-assay-engines/V055.json`. Covers FR-055.
- [ ] T056 [US2] Implement and independently test **proDA fitting and test_diff** in `r/proteomicsCore/R/model_proda.R`. Acceptance V056: Zero-null/reduced-model reference calls match native output and retain model-specific uncertainty and dropout diagnostics. Evidence: command/output/version/hash in `docs/validation/007-assay-engines/V056.json`. Covers FR-056.
- [ ] T057 [US2] Implement and independently test **Alternative-engine uncertainty** in `r/proteomicsCore/R/inference.R`. Acceptance V057: Engine-appropriate intervals/statistics are verified against native/reference calculations; no borrowed limma df or fake TREAT output. Evidence: command/output/version/hash in `docs/validation/007-assay-engines/V057.json`. Covers FR-057.
- [ ] T058 [US2] Implement and independently test **Method sensitivity comparison** in `r/proteomicsCore/R/sensitivity.R; src/proteomics_pipeline/reporting/model_data.py`. Acceptance V058: Matched-universe effects/discoveries are displayed as sensitivities, with unmatched universes explained; no hit-count winner changes primary. Evidence: command/output/version/hash in `docs/validation/007-assay-engines/V058.json`. Covers FR-058.
- [ ] T059 [US2] Implement and independently test **Dependencies and failure propagation** in `src/proteomics_pipeline/doctor.py; r/proteomicsCore/R/capabilities.R`. Acceptance V059: Missing package/backend error is visible, required analysis fails, optional applicability is distinguished from installation failure. Evidence: command/output/version/hash in `docs/validation/007-assay-engines/V059.json`. Covers FR-059.
- [ ] T060 [US2] Implement and independently test **Assay-engine golden/calibration fixtures** in `tests/scientific/test_assay_engines.py; r/proteomicsCore/tests/testthat/test-assay-engines.R`. Acceptance V060: Count-dependent variance and LFQ dropout simulations have pinned reference comparisons, unsuitable-input negative cases and recorded diagnostics. Evidence: command/output/version/hash in `docs/validation/007-assay-engines/V060.json`. Covers FR-060.

## Phase C: Review checkpoint

All checkboxes start unchecked. Complete a checkbox only when Maintainer has read the diff, rerun the relevant tests and recorded acceptance evidence plus commit. A task can span multiple Packet Implementer correction tasks. If a runtime/resource is unavailable, keep its verification pending and continue independent work; do not write a successful fixture-shaped substitute.

## Dependencies & Execution Order

The order above follows shared contract → implementation → integration → adversarial verification. Preserve data-model and method boundaries. Resolve packet prerequisites before dependent tasks and obey the exact allowlist in the frozen packet index. A thin vertical demonstration is a checkpoint, not the end of the parent roadmap.
