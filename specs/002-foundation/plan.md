# Implementation Plan: Runtime, package and contract foundation

**Date:** 2026-09-12. **Spec:** [spec.md](spec.md). **Roadmap:** R01.

## Summary and Technical Context

Implement this slice in the Python CLI / R package architecture from `specs/001-downstream-proteomics/plan.md`. Python owns configuration, filesystem/state and presentation. R owns statistical calculations. Verify installed package APIs against the research contract before using them; do not invent method arguments or output fields.

## Constitution Check

Pass evidence identity, prespecified scientific plan, uncertainty semantics, reproducibility, baseline preservation, Spec Kit traceability, real implementation and reviewed delegation principles. No exception is approved. The normative methodological constraints are in `../001-downstream-proteomics/contracts/scientific-methods.md` and override a convenient but incompatible implementation.

## Dependencies and Interfaces

Prerequisites: recovered baseline and approved parent contracts only; no packet dependency. Before dispatch, Maintainer verifies the baseline commit and reconciles the unverified prior-worktree candidate, then updates the packet starter with the actual helper names, schema versions and gate commands now present. Later packets remain outside this work. Canonical data/result/status fields are shared and cannot be redefined locally.

## Implementation Order

1. Read this spec, tasks and parent contracts; inspect actual source and existing tests.
2. Write meaningful boundary/reference tests for the first task, establish the expected failing behavior, and implement the minimal coherent path.
3. Proceed through tasks in dependency order. A helper shared by multiple tasks belongs at the narrowest common seam, not duplicated in CLI and R.
4. Integrate with existing package/CLI interfaces and test observable end-to-end behavior for this slice.
5. Record command/version/hash evidence per acceptance ID. Maintainer independently reruns it and checks altered tests first.

## Planned Source Seams

- `src/proteomics_pipeline/cli.py; pyproject.toml` — Installable maintained CLI.
- `src/proteomics_pipeline/config.py; schemas/analysis.schema.json` — Strict versioned configuration.
- `src/proteomics_pipeline/runtime.py; schemas/stage-result.schema.json` — Typed run and stage states.
- `scripts/maintained/run_stage.R; src/proteomics_pipeline/runtime.py` — Safe Python-to-R bridge.
- `src/proteomics_pipeline/provenance.py` — Atomic artifact publication.
- `src/proteomics_pipeline/doctor.py` — Environment doctor.
- `r/proteomicsCore/DESCRIPTION; r/proteomicsCore/NAMESPACE; scripts/maintained/test_r.R` — Maintained R package skeleton.
- `tests/contract/; tests/unit/; r/proteomicsCore/tests/testthat/` — Independent test harness.
- `docs/validation/baseline-hashes.json; scripts/maintained/check_baseline.py` — Baseline preservation registry.
- `scripts/maintained/bootstrap.ps1; scripts/maintained/bootstrap.sh; docs/development.md` — One-command developer bootstrap.

## Gate Contract

Current baseline gates: `python -m unittest discover -s tests -v` and `python -m compileall -q pipeline legacy scripts tests`.
After foundation, add actual `python -m pytest tests/unit tests/contract -q` and `Rscript --vanilla scripts/maintained/test_r.R --suite unit` plus this slice's acceptance tests. Scientific adapters additionally require `Rscript --vanilla scripts/maintained/test_r.R --suite scientific` scoped to implemented cases and direct package comparisons. Integration/report slices run their documented canonical example and `proteomics verify`. The first slice creates these target interfaces; they must never be claimed to exist before implementation.

Packet Implementer reports only the packet receipt, including tests and outcomes. Maintainer discovers and runs the actual resulting commands independently. Missing dependencies are an environmental limitation, not a scientific PASS. Do not mark the packet done until its eligible code paths are tested; continue independent work if only an external fixture is unavailable.

## Review, Rollback and Delivery

Maintainer owns final diff review, independent gates, commit and closure; Packet Implementer must not commit or push, and Independent Reviewer never closes. Review any changed baseline tests before trusting green results. Preserve uncommitted implementer work on failure; inspect unstaged, staged and untracked files before cleanup. Changes to scientific meaning return to Maintainer as a spec blocker and require synchronized contract/spec/tests before dependent dispatch. Update progress and traceability only after verification, never from a worker receipt alone.
