# ADR 0007 — Claude route substitution and Maintainer integration amendments

Status: accepted by the operator (Omar, Maintainer) on 2026-10-01.

## Context

AGENTS.md names GPT-5.6 Sol (coordination/verification), GPT-6 Astra (independent audit) and GPT-5.6 Luna (packet implementation) and forbids *silent* route substitution. On 2026-10-01 the operator authorized continuing past R01 (explicit GO for R02 onward, in the repository's own order) and chose **Claude (Opus)** for all three roles. Several frozen interfaces needed small shared changes once real handlers existed.

## Decision

1. **Route.** The Claude (Opus) route performed specification coordination, implementation and verification for R02, R03, R04, R05, R10a, R13, R06, R07, R09, R08, R10b and R11 (R12 tooling only; its private gates are Maintainer-only and NOT_RUN). This substitution is recorded here, in each receipt and in traceability (`verification.route`). The Sol/Astra/Luna routes did not do this work. No separate independent (Astra-equivalent) audit was performed: every receipt records `independent_audit: NOT_RUN`. Acceptance below is Maintainer-route self-verification with re-executed gates, not an independent review.
2. **A-2026-10-01-02 — integration seam.** `src/proteomics_pipeline/workflow.py` (new, Maintainer-owned) implements full validate/plan/run and plan-consistency verify; `cli.py` delegates to it. Planning stages (intake, preprocessing, design) run before a plan exists and carry `plan_hash = null`; `runtime.validate_run_status` accepts that only for those capabilities. The report template is packaged (`pyproject.toml` package-data).
3. **A-2026-10-01-03 — R01 regression tests.** Four R01 tests asserted that limma/plan/run were *not yet implemented*. They now simulate absence explicitly (monkeypatched absent module, an unregistered R capability, a still-unimplemented Phase 2 capability) and keep their invariants: absent handlers are unavailable, unknown R capabilities are NOT_RUN, required unimplemented capabilities exit 3. Doctor's package list may grow with implemented capabilities. R01's historical acceptance at `88e7fb6773ab9403e8b57037a017204357ba75cd` is unchanged.
4. **A-2026-10-01-04 — tooling.** `scripts/check_spec_kit.py` checks the amended kit; `scripts/maintained/install_r_dependencies.R` installs the R packages explicitly into a project-local library (never global, never during analysis).
5. **A-2026-10-01-05 — Phase 2 integration.** Report-section providers for response, pathways and assay engines; the plan freezes DEqMS count-evidence hashes; the R05 score-test guard applies only to Phase 1 plans; the workflow runs assay engines, resources, response and pathways; `resources prepare`; the R01 doctor test simulates an absent resources module.
6. **A-2026-10-01-06 — full report and comparison.** `report --run` and `compare --left --right --output` are activated as contracted; the workflow appends an optional `report_full` stage after the unchanged R10a stub; the R10a stage-order assertion allows that trailing stage.
7. **A-2026-10-01-07 — resume and cache.** `resume --run` is activated; stage reuse requires an identical fingerprint of normalized request, code manifest and environment, and any mismatch re-executes the stage and its dependants.
8. **A-2026-10-01-08 — R04 regression test.** The V038 negative asserted that `assay_engines` was unimplemented; after R06 it simulates an absent `assay_service` module (invariant unchanged). Found by the full gate re-run on 2026-10-02.
9. **A-2026-10-01-09 — R01 CLI regression test.** The required-but-unimplemented example was a DEqMS model; after R06 the child process blocks `assay_service` explicitly (invariant unchanged). The re-run also exposed an untyped error for a missing count-evidence file, fixed inside R06.

## Consequences

Commits, pushes, PRs, remotes, licensing and GitHub/Notion changes remain unauthorized; all work is an uncommitted working tree on branch `claude/full-pipeline`, so every `verified_commit` stays null and receipts identify working-source manifests instead.
