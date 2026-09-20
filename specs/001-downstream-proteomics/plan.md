# Implementation Plan: Maintained downstream proteomics

**Date:** 2026-09-12. **Spec:** [spec.md](spec.md). **Status:** proposed implementation architecture.

## Summary

Build a thin Python CLI around a tested R analysis package. Preserve the old scripts as immutable reference. One schema-driven plan dispatches reusable stages with explicit scientific applicability and file contracts. The execution layer owns hashing, status, caching, report inputs and reproducibility; R owns statistical estimation, diagnostics, multiplicity and enrichment. Avoid parallel implementations of the same statistics in Python and R. Independent Python calculations are test oracles only.

## Technical Context

**Languages:** Python 3.12 and an R/Bioconductor pair resolved and pinned together at implementation. Preserve a separate R 4.4.2 legacy reproduction environment when feasible. Do not force legacy packages into the maintained environment.
**Python dependencies:** packaging/CLI, schema validation and YAML parsing; pandas/NumPy or Arrow only when used by real data contracts; Jinja2 for offline report rendering. Prefer Typer or argparse consistently, not both. Resolve exact versions in the first slice and record the choice.
**R dependencies:** limma, edgeR only when directly required, DEqMS, proDA, fgsea, BiocParallel, statmod, matrixStats, jsonlite, readr, ggplot2 and testthat, plus a minimal declared table/plot support set. Resource preparation may use msigdbr; production analysis consumes snapshots. Add packages only after checking their installed API and role.
**Storage:** local TSV/CSV/JSON, optional Parquet for large numeric artifacts, RDS for internal fitted objects. TSV is the interoperable scientific export. No database.
**Environment:** Python lock plus R renv lock; pinned container build with digest available for Linux. Native Windows and Linux mandatory for core paths; macOS documented and tested if a runner is available. Core Python checks alone are not a platform pass for R.
**Testing:** pytest, existing unittest regression, testthat, golden reference calls, independent analytic calculations, Monte Carlo calibration, report content and accessibility checks, clean installation and resource-offline tests.
**Constraints:** private archive remains on SSD; no private uploads; no runtime installation/downloads; no uncontrolled shell interpolation or arbitrary R eval from config; immutable execution directories.

## Constitution Check

| Principle | Architectural enforcement |
|---|---|
| Evidence first | Canonical schema, biological-unit hierarchy, input-scale state machine, unknown provenance fields |
| Executable scientific plan | Read-only validation then frozen semantic plan before fitting; method applicability matrix |
| Uncertainty survives | Typed hypotheses, backend-native inference, exact contrast handling, explicit families and nonestimability |
| Reproducibility verified | Content-addressed stage keys, per-process session capture, locked resources, independent gates |
| Baseline preserved | Historical directories excluded from implementation edits and checked by hash |
| Spec Kit traceability | Numbered slices, 120 linked tasks, ADR and evidence workflow |
| Coherent interfaces | One Python package plus one R package; file contract only across language boundary |
| Reviewed delegation | Frozen one-packet Packet Implementer tasks with review/test/commit checkpoints and closure by Maintainer; Independent Reviewer audits only when the user requests it |

No constitution exception is currently approved. If method research changes an API assumption, update the relevant contract/ADR and tests before implementing that adapter.

## Maintained Source Structure

```text
src/proteomics_pipeline/
  cli.py                  command dispatch, useful errors and exit status
  config.py               schema + semantic validation
  intake/                 canonical matrix, vendor mappings and sample hierarchy
  planning.py             design/method eligibility and frozen plan
  runtime.py              subprocess, atomic stages, logs and resource limits
  provenance.py           content hashes, manifests, snapshots and semantic comparison
  reporting/              renderer and report-data schema (no statistics)
r/proteomicsCore/
  DESCRIPTION NAMESPACE
  R/io.R R/validation.R R/preprocess.R R/qc.R R/design.R
  R/model_limma.R R/model_deqms.R R/model_proda.R
  R/inference.R R/multiplicity.R R/mapping.R R/pathways.R
  R/response.R R/independent_scores.R R/diagnostics.R
  tests/testthat/
scripts/run_stage.R       maintained bridge in new scripts/maintained/ path if name conflict
schemas/                 versioned canonical schemas copied from approved contracts
configs/                 runnable synthetic examples and annotated templates
resources/               small redistributable fixture sets + metadata, not private datasets
tests/unit/ tests/contract/ tests/integration/ tests/scientific/ tests/regression/
docs/user-guide/ docs/methods/ docs/adr/ docs/validation/
pipeline/ legacy/ docs/audit/  immutable historical reference
```

The existing `scripts/run_stage.R` is part of the recovery helper and may be superseded through a new maintained path rather than inadvertently changing its archival behavior. Python module paths in slice plans are intended seams; a small, documented rename is permitted before dependent tasks begin, with traceability updates.

## Execution DAG and State

`ingest → validate → plan → preprocess/QC → design → fit → inference → mapping/pathways → response → report → verify`.

QC may produce a report before modeling. Pre-fit QC cannot silently exclude observations. All exclusions require frozen policy and recorded rationale. A derived exploratory reanalysis has a new plan, never mutates the primary fit. Pathway/response stages may run independently once their model inputs are available; initial implementation runs deterministically/sequentially. Add controlled parallelism only after reproducibility tests.

States: PENDING, RUNNING, COMPLETED, COMPLETED_WITH_WARNINGS, INAPPLICABLE, FAILED, CANCELLED. A required FAILED/NOT_RUN stage prevents overall completion. INAPPLICABLE requires an explicit eligibility reason and contract permitting omission for that study. Absence of significance is COMPLETED with zero discoveries, not INAPPLICABLE. Record public semantic run status separately from OS return code.

Each stage writes into a temporary directory beneath its run, validates output schema/integrity and atomically promotes it. Cache key includes code hash, schema version, normalized config subset, source hashes, sample/feature ordered IDs, packages/container digest, RNG settings and resource snapshot hashes. A changed dependency invalidates all downstream affected stages. Reuse is permitted only for completed verified artifacts, never a partial directory. Lock concurrent writers to a run and test interrupted-lock recovery; do not use PID existence alone across hosts.

## Baseline Migration

1. Register archive and script hashes. Retain baseline ten tests unchanged.
2. Implement new canonical intake and a local `import-legacy` converter driven by stable feature/sample identity, not D1/T1 prefixes.
3. Add new model/inference independently; compare matched feature universes and configurations to historical outputs.
4. Label scientific corrections separately: global declared multiplicity, robust design handling, independent score testing, exhaustive enumeration, bounded response classes and correlation-aware pathway methods.
5. Regenerate a new report with all principal and auxiliary contrasts. Never overwrite the historical report or use its text as executable truth.
6. Package a successor release on SSD. Maintainer must produce a manifest and transfer only reviewed source/results authorized for local storage, not alter the historical package manifest to pretend it was the original delivery.

## Environments and Command Gates

Existing baseline gates run now: `python -m unittest discover -s tests -v`; `python -m compileall -q pipeline legacy scripts tests`. They validate the recovered baseline only, not the future R01 implementation.

The foundation slice must make these future commands real and document exact environment activation:

```text
python -m pytest tests/unit tests/contract -q
Rscript --vanilla scripts/maintained/test_r.R --suite unit
python -m proteomics_pipeline.cli doctor --json
python -m proteomics_pipeline.cli validate --config configs/examples/synthetic-independent.yml
python -m proteomics_pipeline.cli plan --config configs/examples/synthetic-independent.yml --output runs/example-plan.json
python -m proteomics_pipeline.cli run --config configs/examples/synthetic-independent.yml --output runs/example-01
python -m proteomics_pipeline.cli verify --run runs/example-01
python -m pytest tests/integration tests/scientific -q
Rscript --vanilla scripts/maintained/test_r.R --suite scientific
python scripts/maintained/run_calibration.py --profile release --output runs/calibration
python scripts/maintained/legacy_regression.py --archive-root <local-private-original> --output runs/legacy-comparison
```

These are target contracts, not assertions that commands exist today. Before each Packet Implementer packet starter Maintainer discovers the gate commands actually implemented and lists newly required gates separately. No task may pass by omitting a gate because its command was never implemented.

## Complexity and Delivery Strategy

Twelve independently reviewable Spec Kit packets share these contracts. The corrected graph and exact file ownership are in `packet-index.md`; only R06, R07 and R09 may run in parallel after R05, because their write sets are disjoint and they have no mutual dependencies. Maintainer verifies changes, records evidence and commits; no unattended or receipt-only batch is considered landed. The final gate checks that the whole command path and scientific claim vocabulary remain coherent, not just individual modules.
