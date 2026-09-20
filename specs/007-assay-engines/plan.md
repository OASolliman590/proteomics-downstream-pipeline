# Implementation Plan: Assay-qualified DEqMS and proDA backends

**Date:** 2026-09-12. **Spec:** [spec.md](spec.md). **Roadmap:** R06.

## Summary and Technical Context

Implement this slice in the Python CLI / R package architecture from `specs/001-downstream-proteomics/plan.md`. Python owns configuration, filesystem/state and presentation. R owns statistical calculations. Verify installed package APIs against the research contract before using them; do not invent method arguments or output fields.

## Constitution Check

Pass evidence identity, prespecified scientific plan, uncertainty semantics, reproducibility, baseline preservation, Spec Kit traceability, real implementation and reviewed delegation principles. No exception is approved. The normative methodological constraints are in `../001-downstream-proteomics/contracts/scientific-methods.md` and override a convenient but incompatible implementation.

## Dependencies and Interfaces

Prerequisites: R05 (`specs/006-limma-inference/`). Before dispatch, Maintainer verifies prerequisite commits and updates the packet starter with the actual helper names, schema versions and gate commands now present. Later packets remain outside this work. Canonical data/result/status fields are shared and cannot be redefined locally.

## Implementation Order

1. Read this spec, tasks and parent contracts; inspect actual source and existing tests.
2. Write meaningful boundary/reference tests for the first task, establish the expected failing behavior, and implement the minimal coherent path.
3. Proceed through tasks in dependency order. A helper shared by multiple tasks belongs at the narrowest common seam, not duplicated in CLI and R.
4. Integrate with existing package/CLI interfaces and test observable end-to-end behavior for this slice.
5. Record command/version/hash evidence per acceptance ID. Maintainer independently reruns it and checks altered tests first.

## Planned Source Seams

- `schemas/count-evidence.schema.json; src/proteomics_pipeline/intake/counts.py` — Count evidence schema.
- `r/proteomicsCore/R/model_deqms.R` — DEqMS fit pipeline.
- `r/proteomicsCore/R/model_deqms.R` — DEqMS statistic alignment.
- `r/proteomicsCore/R/capabilities.R` — DEqMS hypotheses and design guard.
- `src/proteomics_pipeline/planning.py; r/proteomicsCore/R/model_proda.R` — LFQ dropout provenance gate.
- `r/proteomicsCore/R/model_proda.R` — proDA fitting and test_diff.
- `r/proteomicsCore/R/inference.R` — Alternative-engine uncertainty.
- `r/proteomicsCore/R/sensitivity.R; src/proteomics_pipeline/reporting/model_data.py` — Method sensitivity comparison.
- `src/proteomics_pipeline/doctor.py; r/proteomicsCore/R/capabilities.R` — Dependencies and failure propagation.
- `tests/scientific/test_assay_engines.py; r/proteomicsCore/tests/testthat/test-assay-engines.R` — Assay-engine golden/calibration fixtures.

## Gate Contract

Current baseline gates: `python -m unittest discover -s tests -v` and `python -m compileall -q pipeline legacy scripts tests`.
After foundation, add actual `python -m pytest tests/unit tests/contract -q` and `Rscript --vanilla scripts/maintained/test_r.R --suite unit` plus this slice's acceptance tests. Scientific adapters additionally require `Rscript --vanilla scripts/maintained/test_r.R --suite scientific` scoped to implemented cases and direct package comparisons. Integration/report slices run their documented canonical example and `proteomics verify`. The first slice creates these target interfaces; they must never be claimed to exist before implementation.

Packet Implementer reports only the packet receipt, including tests and outcomes. Maintainer discovers and runs the actual resulting commands independently. Missing dependencies are an environmental limitation, not a scientific PASS. Do not mark the packet done until its eligible code paths are tested; continue independent work if only an external fixture is unavailable.

## Review, Rollback and Delivery

Maintainer owns final diff review, independent gates, commit and closure; Packet Implementer must not commit or push, and Independent Reviewer never closes. Review any changed baseline tests before trusting green results. Preserve uncommitted implementer work on failure; inspect unstaged, staged and untracked files before cleanup. Changes to scientific meaning return to Maintainer as a spec blocker and require synchronized contract/spec/tests before dependent dispatch. Update progress and traceability only after verification, never from a worker receipt alone.
