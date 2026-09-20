# Implementation Plan: Unified offline report and streamlined workflow

**Date:** 2026-09-12. **Spec:** [spec.md](spec.md). **Roadmap:** R10.

## Summary and Technical Context

Implement this slice in the Python CLI / R package architecture from `specs/001-downstream-proteomics/plan.md`. Python owns configuration, filesystem/state and presentation. R owns statistical calculations. Verify installed package APIs against the research contract before using them; do not invent method arguments or output fields.

## Constitution Check

Pass evidence identity, prespecified scientific plan, uncertainty semantics, reproducibility, baseline preservation, Spec Kit traceability, real implementation and reviewed delegation principles. No exception is approved. The normative methodological constraints are in `../001-downstream-proteomics/contracts/scientific-methods.md` and override a convenient but incompatible implementation.

## Dependencies and Interfaces

Prerequisites: R03 (`specs/004-preprocessing-qc/`), R06 (`specs/007-assay-engines/`), R08 (`specs/009-enrichment/`) and R09 (`specs/010-treatment-response/`). Before dispatch, Maintainer verifies prerequisite commits and updates the packet starter with the actual helper names, schema versions and gate commands now present. Later packets remain outside this work. Canonical data/result/status fields are shared and cannot be redefined locally.

## Implementation Order

1. Read this spec, tasks and parent contracts; inspect actual source and existing tests.
2. Write meaningful boundary/reference tests for the first task, establish the expected failing behavior, and implement the minimal coherent path.
3. Proceed through tasks in dependency order. A helper shared by multiple tasks belongs at the narrowest common seam, not duplicated in CLI and R.
4. Integrate with existing package/CLI interfaces and test observable end-to-end behavior for this slice.
5. Record command/version/hash evidence per acceptance ID. Maintainer independently reruns it and checks altered tests first.

## Planned Source Seams

- `src/proteomics_pipeline/reporting/data.py; schemas/report.schema.json` — Typed report data assembler.
- `src/proteomics_pipeline/cli.py; src/proteomics_pipeline/runtime.py` — One-run command integration.
- `src/proteomics_pipeline/reporting/templates/; src/proteomics_pipeline/reporting/render.py` — Offline HTML report.
- `src/proteomics_pipeline/reporting/qc_data.py` — QC and inclusion/exclusion sections.
- `src/proteomics_pipeline/reporting/model_data.py` — Complete differential results views.
- `src/proteomics_pipeline/reporting/pathway_data.py; response_data.py` — Pathway and response views.
- `r/proteomicsCore/R/plots.R` — Publication figure and source exports.
- `src/proteomics_pipeline/reporting/methods.py` — Methods and limitations generation.
- `tests/integration/test_reports.py` — Null/empty/partial report robustness.
- `docs/user-guide/; configs/examples/` — Quickstart and user workflow.

## Gate Contract

Current baseline gates: `python -m unittest discover -s tests -v` and `python -m compileall -q pipeline legacy scripts tests`.
After foundation, add actual `python -m pytest tests/unit tests/contract -q` and `Rscript --vanilla scripts/maintained/test_r.R --suite unit` plus this slice's acceptance tests. Scientific adapters additionally require `Rscript --vanilla scripts/maintained/test_r.R --suite scientific` scoped to implemented cases and direct package comparisons. Integration/report slices run their documented canonical example and `proteomics verify`. The first slice creates these target interfaces; they must never be claimed to exist before implementation.

Packet Implementer reports only the packet receipt, including tests and outcomes. Maintainer discovers and runs the actual resulting commands independently. Missing dependencies are an environmental limitation, not a scientific PASS. Do not mark the packet done until its eligible code paths are tested; continue independent work if only an external fixture is unavailable.

## Review, Rollback and Delivery

Maintainer owns final diff review, independent gates, commit and closure; Packet Implementer must not commit or push, and Independent Reviewer never closes. Review any changed baseline tests before trusting green results. Preserve uncommitted implementer work on failure; inspect unstaged, staged and untracked files before cleanup. Changes to scientific meaning return to Maintainer as a spec blocker and require synchronized contract/spec/tests before dependent dispatch. Update progress and traceability only after verification, never from a worker receipt alone.
