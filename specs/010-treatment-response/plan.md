# Implementation Plan: Treatment response, equivalence and independent scores

**Date:** 2026-09-12. **Spec:** [spec.md](spec.md). **Roadmap:** R09.

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

- `r/proteomicsCore/R/response.R` — Explicit disease-treatment-residual axes.
- `r/proteomicsCore/R/response.R` — Descriptive reversal without double dipping.
- `r/proteomicsCore/R/response.R` — Bounded response categories.
- `r/proteomicsCore/R/response_uncertainty.R` — Ratio uncertainty contract.
- `r/proteomicsCore/R/equivalence.R` — Model-correct residual equivalence.
- `r/proteomicsCore/R/equivalence.R` — Conjunction claim eligibility.
- `schemas/independent-score.schema.json; r/proteomicsCore/R/scores.R` — Independent fixed-score artifacts.
- `r/proteomicsCore/R/randomization.R` — Exact and Monte Carlo randomization.
- `src/proteomics_pipeline/reporting/response_data.py` — Score/equivalence reporting semantics.
- `tests/scientific/test_response.py; r/proteomicsCore/tests/testthat/test-response.R` — Circularity/null/overshoot tests.

## Gate Contract

Current baseline gates: `python -m unittest discover -s tests -v` and `python -m compileall -q pipeline legacy scripts tests`.
After foundation, add actual `python -m pytest tests/unit tests/contract -q` and `Rscript --vanilla scripts/maintained/test_r.R --suite unit` plus this slice's acceptance tests. Scientific adapters additionally require `Rscript --vanilla scripts/maintained/test_r.R --suite scientific` scoped to implemented cases and direct package comparisons. Integration/report slices run their documented canonical example and `proteomics verify`. The first slice creates these target interfaces; they must never be claimed to exist before implementation.

Packet Implementer reports only the packet receipt, including tests and outcomes. Maintainer discovers and runs the actual resulting commands independently. Missing dependencies are an environmental limitation, not a scientific PASS. Do not mark the packet done until its eligible code paths are tested; continue independent work if only an external fixture is unavailable.

## Review, Rollback and Delivery

Maintainer owns final diff review, independent gates, commit and closure; Packet Implementer must not commit or push, and Independent Reviewer never closes. Review any changed baseline tests before trusting green results. Preserve uncommitted implementer work on failure; inspect unstaged, staged and untracked files before cleanup. Changes to scientific meaning return to Maintainer as a spec blocker and require synchronized contract/spec/tests before dependent dispatch. Update progress and traceability only after verification, never from a worker receipt alone.
