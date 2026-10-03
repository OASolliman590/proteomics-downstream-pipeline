# Phase-gate regression evidence (2026-10-02, Claude route; re-run after the independent-audit fixes, ADR 0008)

Re-executed after R02–R11 and R13 were implemented (and after amendments A-2026-10-01-08/09 fixed two regression tests that assumed DEqMS was still unimplemented), on macOS x86_64, Python 3.13.15, R 4.6.1, project-local `.venv/` and `.r-lib/` (restored from `renv.lock`), `LANG=en_US.UTF-8`. Logs are in `evidence/` with local paths redacted. Every packet gate was also re-executed the same day; see each packet's `evidence/gate-results.json` and the ledger `docs/validation/acceptance.json`.

| Gate | Command | Result |
|---|---|---|
| Historical recovered tests | `.venv/bin/python -m unittest discover -s tests -v` | 10/10 OK |
| R01 regression (with recorded test amendments A-03, A-09; config change A-10) | `.venv/bin/python -m pytest tests/contract/foundation tests/unit/foundation -q -rs` | 70 passed |
| All R testthat files | `Rscript --vanilla scripts/maintained/test_r.R --suite foundation` | 319 expectations, 0 failures (assay-engines 13, design 28, foundation 46, limma 20, mapping 11, pathways 19, permanova 87, preprocessing 44, qc 16, response 35) |
| Protected baseline registry | `.venv/bin/python scripts/maintained/check_baseline.py` | valid, 0 mismatches (pipeline/, legacy/, docs/audit/ unchanged) |
| Spec kit structure/traceability/links | `.venv/bin/python scripts/check_spec_kit.py` | PASS (13 slices, 130 requirements) |
| Whole Python suite | `.venv/bin/python -m pytest tests -q -rs` | 266 passed, 0 skipped |

Accepted by the operator-authorized Claude route (self-verified; fd8dfa9 independently audited ACCEPT WITH FIXES, fixes committed as fa4d948 and re-reviewed ACCEPT WITH FIXES): Phase 1 R01–R05 + R10a; Phase 2 R13, R06, R07, R09, R08; Phase 3 R10b and R11 (V107 by CI run 37087074374). Ledger: 120 PASS, 10 NOT_RUN (V111–V120 R12 private regression/release, Maintainer-only). fd8dfa9 was committed by Omar and independently audited (ACCEPT WITH FIXES); the fixes are committed as fa4d948 and were re-reviewed (ACCEPT WITH FIXES, non-blocking follow-ups). Nothing is tagged or released.
