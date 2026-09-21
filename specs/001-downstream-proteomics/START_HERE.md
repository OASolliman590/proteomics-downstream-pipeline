# Start here — Spec Kit v1.2.0-frozen

This checkout's baseline is recovered study code at `3aefdd95c46a1c56dae89e79ddf441b6cc978148`. The maintained successor is specified, not implemented or validated. Independent audit and Maintainer freeze are outstanding; no implementation packet is authorized.

Read [PHASES.md](PHASES.md), then [packet-index.md](packet-index.md). The only first implementation brief is [IMPLEMENTATION_BRIEF_R01.md](IMPLEMENTATION_BRIEF_R01.md), used only after freeze. The [constitution](../../.specify/memory/constitution.md), [parent spec](spec.md), [scientific methods](contracts/scientific-methods.md), [schema](contracts/analysis.schema.json), [semantic validation](contracts/semantic-validation.md), [data model](data-model.md) and [CLI/artifacts](contracts/cli-and-artifacts.md) are the shared contracts.

The [roadmap](roadmap.md), [plan](plan.md), [validation strategy](validation-strategy.md), [traceability](traceability.json), [progress](progress.md) and [quality checklist](checklists/spec-quality.md) separate required software from actual evidence. [Research](research.md) explains changed contract decisions with primary references; older scientific-methods-research.md is background, not a competing normative contract.

Current runnable surface remains scripts/run_pipeline.py, the six historical R stages, pipeline/20_external_limma_intake.py and tests/test_runner.py. The recovery did not rerun R; environment.yml is not a lockfile; resource snapshots and a repository license remain unresolved. No private study files are needed or permitted for this candidate review.

[Hardening report](../../docs/SPEC_HARDENING_REPORT.md) records the changes and verification limits. It is outside the immutable docs/audit/ tree. R01–R05 plus thin R10a is a future Phase 1 milestone, never proof that full R10/R11/R12 or v1.0 shipped.
