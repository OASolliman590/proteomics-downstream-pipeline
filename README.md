# Proteomics downstream pipeline

This repository contains a **recovered study baseline** and a **maintained successor in progress**, governed by a frozen v1.2.0 Spec Kit (with operator-authorized 2026-10-01 amendments). It is **not yet a released or independently audited platform**. Private study files are not included and must not be uploaded or committed.

## What exists today

| Scope | Status |
|---|---|
| Recovered study scripts + immutable audit (`pipeline/`, `legacy/`, `docs/audit/`) | Inspectable baseline; unchanged. |
| R01 foundation (CLI, configuration, runtime, Python↔R I/O) | Accepted 2026-09-21. |
| R02 intake, R03 QC, R04 design, R05 limma, R10a thin report (Phase 1, `v0.1-limma-core` scope) | Implemented; all acceptance cases pass on synthetic fixtures (2026-10-01), **uncommitted** and self-verified by the operator-chosen Claude route; independent audit not yet done. |
| R13 PERMANOVA/PERMDISP (operator-authorized amendment, Phase 2) | Implemented; acceptance cases pass on synthetic fixtures (same caveats). |
| R06 DEqMS/proDA, R07 resources/mapping, R09 response, R08 pathways (Phase 2); R10b full report/compare; R11 locks, offline reproduction, resume, calibration, benchmark, validation ledger (Phase 3) | Implemented; acceptance cases pass on synthetic fixtures (same caveats). R11 V107 cross-platform CI passes on Ubuntu and Windows (CI run 37087074374). |
| R12 private regression and release | NOT_RUN (Maintainer-only). See [progress](specs/001-downstream-proteomics/progress.md). |

The maintained target is declared protein-abundance intake/QC → frozen design → qualified limma (and eligible adapters) → design-valid enrichment → descriptive response and restricted independent inference → honest offline reporting. It excludes raw-MS search/quantification, PTM localization, single-cell, classifier training, web UI/databases, network deployment and causal drug-mechanism claims.

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

Start with [START_HERE](specs/001-downstream-proteomics/START_HERE.md), [PHASES](specs/001-downstream-proteomics/PHASES.md), the [packet index](specs/001-downstream-proteomics/packet-index.md), the [scientific methods contract](specs/001-downstream-proteomics/contracts/scientific-methods.md) (v1.2.0) and [decisions and open questions](specs/001-downstream-proteomics/decisions.md). No phase is complete merely because its specification exists; acceptance is recorded only from executed evidence in [traceability](specs/001-downstream-proteomics/traceability.json). No license has been selected.
