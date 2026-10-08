# Proteomics downstream pipeline

This repository contains a **recovered study baseline** and a **maintained successor**, governed by a frozen Spec Kit (v1.2.0 core, with the operator-authorized amendments of ADR 0006 and ADR 0009). It is **not a released platform** and not a v1.0 release. Private study files are not included and must not be uploaded or committed.

## What exists today

| Scope | Status |
|---|---|
| Recovered study scripts + immutable audit (`pipeline/`, `legacy/`, `docs/audit/`) | Inspectable baseline; unchanged (hash-checked). |
| R01 foundation (CLI, configuration, runtime, Python↔R I/O) | Accepted 2026-09-21. |
| R02–R11 (intake, QC, design, limma, DEqMS/proDA, resources/mapping, pathways, response, thin and full reports, locks/reproduction/calibration/validation ledger) and R13 PERMANOVA | Merged into `main` (PR #1). Acceptance cases pass on synthetic fixtures. An independent audit (2026-10-02) and re-review (2026-10-03) returned ACCEPT WITH FIXES; the fixes are done. V107 cross-platform CI passed on Ubuntu and Windows (run 37087074374). |
| R14a–R14f post-differential analysis (spec 015, Phase 4) | On branch `claude/post-de` (PR #2): V131–V167 pass on synthetic fixtures. The post-DE work was independently reviewed on 2026-10-05 (ACCEPT WITH FIXES, fixed). Modules adapt to the data and never hold up the DE results ([ADR 0010](docs/adr/0010-adaptive-post-de-policy.md)); every adaptation is recorded. The adaptive redesign itself has not yet been independently reviewed. Post-DE calibration evidence: `docs/validation/015-post-de-analysis/R14f/calibration`. |
| R12 release tooling and documentation | V115–V117 pass (methods/reason-code docs, spec-task-evidence reconciliation, release hygiene). V111–V114 and V118–V120 are NOT_RUN: Maintainer-only private gates, operator decisions (license, publication) and final review. See the [Maintainer runbook](docs/validation/013-release/MAINTAINER_RUNBOOK.md). |

Every status comes from executed evidence in the [acceptance ledger](docs/validation/acceptance.json); `verified_commit` stays null until CI runs on a pushed commit. Methods and capabilities: [docs/methods](docs/methods/README.md). Release rules: [docs/user-guide/release.md](docs/user-guide/release.md). Changes: [CHANGELOG](CHANGELOG.md).

The maintained target is declared protein-abundance intake/QC → frozen design → qualified limma (and eligible adapters) → design-valid enrichment → descriptive response and restricted independent inference → honest offline reporting. Phase 4 adds opt-in post-differential analysis (sets, sensitivity, phenotype association, leakage-safe biomarker evaluation, co-abundance and interaction networks). It excludes raw-MS search/quantification, PTM localization, single-cell, web UI/databases, network deployment, clinical/diagnostic claims and causal drug-mechanism claims.

## Quick start (synthetic examples only)

```bash
# Python (project-local) and R packages (project-local library)
~/.local/bin/uv venv --python 3.13 .venv && ~/.local/bin/uv pip install --python .venv/bin/python -e '.[test]'
export R_LIBS_USER="$PWD/.r-lib" LANG=en_US.UTF-8
Rscript --vanilla scripts/maintained/install_r_dependencies.R
R CMD INSTALL --library="$R_LIBS_USER" r/proteomicsCore

# Validate, plan and run the Phase 1 example, then verify integrity
.venv/bin/proteomics validate --config configs/examples/example-independent.json --json
.venv/bin/proteomics run --config configs/examples/example-independent.json --output runs/example --json
.venv/bin/proteomics verify --run runs/example --json
open runs/example/report/index.html
```

Tests: `.venv/bin/python -m pytest tests -q -rs` and `Rscript --vanilla scripts/maintained/test_r.R --suite foundation` (runs every testthat file). Full workflow, Phase 2 demo, resume, locks and validation tooling: [docs/user-guide/usage.md](docs/user-guide/usage.md). Packet receipts and evidence are under [docs/validation](docs/validation/phase-gates/README.md).

The historical runner remains documented in [docs/PIPELINE.md](docs/PIPELINE.md) (`python scripts/run_pipeline.py --help`, `python -m unittest discover -s tests -v`).

## Specification

Start with [START_HERE](specs/001-downstream-proteomics/START_HERE.md), [PHASES](specs/001-downstream-proteomics/PHASES.md), the [packet index](specs/001-downstream-proteomics/packet-index.md), the [scientific methods contract](specs/001-downstream-proteomics/contracts/scientific-methods.md) (v1.3.0) and [decisions and open questions](specs/001-downstream-proteomics/decisions.md). No phase is complete merely because its specification exists; acceptance is recorded only from executed evidence in [traceability](specs/001-downstream-proteomics/traceability.json). No license has been selected.
