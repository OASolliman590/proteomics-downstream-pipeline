# R12 Maintainer-route receipt — release tooling, documentation and reconciliation

Packet R12 (FR-111–FR-120 / T111–T120 / V111–V120) was worked on 2026-10-04 by the **Claude (Opus) route**, operator-authorized for the Sol/Astra/Luna roles ([ADR 0007](../../adr/0007-claude-route-and-integration-amendments.md)). Independent audit: **NOT_RUN**. The packet is **partial by design**: the private gates and operator decisions are not an implementer's to run or make.

## Identities

- Contracts: kit 1.2.0 core with ADR 0006/0009 amendments; `specs/013-release/spec.md` (SM12, SM14, SM19, SM25).
- Working-source manifest (R12 allowlist plus the amended `reproduction.py`, 15 files): `6e0260c0c30cbfb8e3ef2c6e52451279862b9d91a5d67f015dd3069d95b32879` — `evidence/working-source-manifest.txt`.
- Shared change: A-2026-10-01-21 (the ledger keeps each R12 row's recorded reason; the traceability identity note records the ADR 0009 identities). Decision D-53.
- `verified_commit`: **null** (local commit on `claude/post-de`, not pushed). Environment: macOS x86_64, Python 3.13.15, git, `LANG=en_US.UTF-8`.

## Gates executed

| Command | Exit | Result | Log |
|---|---|---|---|
| `.venv/bin/python -m pytest tests/regression -q -rs -p no:cacheprovider` | 0 | 38 passed | `evidence/python-tests-1.log` |
| `.venv/bin/python scripts/maintained/build_release.py check` | 0 | no blocking finding on this tree; license and ownership unresolved | `evidence/release-hygiene-check.log` |
| `.venv/bin/python scripts/check_spec_kit.py` | 0 | PASS | `evidence/spec-kit-check.log` |
| `.venv/bin/python -m pytest tests -q -p no:cacheprovider -rs` (whole tree incl. R14a–R14f, calibration evidence and R12) | 0 | 394 passed, 0 skipped | `evidence/full-python-suite.log` |
| all R testthat files (`testthat::test_dir(..., load_package = "installed")`) | 0 | 677 expectations, 0 failed, 0 skipped | `evidence/r-all-testthat.log` |
| `Rscript --vanilla scripts/maintained/test_r.R --suite foundation` | 0 | every testthat file, no failure | `evidence/r-test-suite-foundation.log` |
| R11 lock qualification (`check_r_lock('renv.lock')`, `check_python_lock('requirements.lock')`) | 0 | qualified; no Python lock problem | `evidence/lock-qualification.log` |

The first full-suite run found one failure: the R11 encoding lint flagged `ZipFile.open` in `legacy_service.py` as a text-mode `open()`. Archive members are now hashed from bytes read in memory (nothing is extracted), and every gate was re-run on the corrected tree (the results above).

## Acceptance cases

| Case | Status | Evidence / reason |
|---|---|---|
| V111 | NOT_RUN | Maintainer-only private regression; run locally by Omar ([runbook](MAINTAINER_RUNBOOK.md) §1). The tooling (`legacy_regression.py verify-archive`) is tested on a synthetic archive: hashes verified against an inventory without extraction; a tampered member fails; a missing archive is NOT_RUN (`E_LEGACY_UNAVAILABLE`). |
| V112 | NOT_RUN | Maintainer-only private regression; run locally by Omar. It needs historical tables and a provisioned historical R/resource environment (runbook §3). The comparison tooling is tested on hand-made tables: frozen tolerances, preserved differences, NOT_REPRODUCED rows, refusal of copied values (`E_LEGACY_COPIED`) and of unverified runs. |
| V113 | NOT_RUN | Needs Maintainer-approved local comparison summaries and a disclosure review. The implementer route did not write study-derived numbers into the public repository (runbook §4). |
| V114 | NOT_RUN | Maintainer-only private regression; run locally by Omar (runbook §2). |
| V115 | PASS | `tests/regression/test_docs.py`: the [methods page](../../methods/README.md) lists exactly the implemented capabilities (runtime discovery) and every schema enum value. The [reason-code catalogue](../../methods/reason-codes.md) equals the codes found in the sources in both directions and is regenerated, not hand-edited. Negative: no license claim in public docs, no license file, "not a v1.0 release" stated, and a fake capability row is detected. The documented commands are executed by V100. |
| V116 | PASS | `tests/regression/test_traceability_reconciliation.py`: identities parsed independently from every slice's spec/tasks equal the 167 traceability rows: exactly 120 core identities plus the ADR 0006 and 0009 amendments. Checkboxes agree with statuses, every PASS has existing non-pending evidence, every unresolved row has an explicit reason, and the ledger agrees in status and reason. Eight negative mutations are detected (pending evidence, rename without redirect, duplicate, missing spec identity, missing reason, checkbox drift, redirect without ADR, generic ledger reason). Fail-before for A-2026-10-01-21: `evidence/ledger-reason-fail-before.log` (10 R12 reasons overwritten). |
| V117 | PASS | `tests/regression/test_release.py`: a scratch git tree with seeded dummy secrets (private key, GitHub, AWS, password), private paths (Unix and Windows), a private workbook, `.env`, a private log and a fake MIT license. Each one blocks the release and the export, while untracked environments, runs and private evidence are never inventoried (independent inventory oracle). A changed protected baseline blocks the release; CI runner paths are allowed. This repository has no blocking finding. |
| V118 | NOT_RUN | operator decision required (AGENTS.md): no release authorization and no license decision exist. The tooling is tested: a synthetic private-successor export gives byte-identical files and an independently verified `MANIFEST.sha256` with the source unchanged. Collision, missing or mismatched authorization, a public release with an unresolved license, a destination inside the tree and a bad version are each refused. |
| V119 | NOT_RUN | Needs the final reviewed tree, CI on its pushed commit and V111/V114. The local re-run of every gate on this tree is recorded in [final-acceptance.md](../final-acceptance.md) and cannot support v1.0 completion. |
| V120 | NOT_RUN | Needs the Maintainer's review of all changes, evidence and authorization records; a worker receipt cannot close the roadmap. |

## Not done by this route (on purpose)

No license was chosen, nothing was published or pushed, and no private data was requested or read. `docs/validation/legacy-comparison.md` (V113) and `docs/validation/delivery.json` (V120) are left for the Maintainer.

## Acceptance

All executed gates exited 0. V115–V117 are accepted by the operator-authorized Claude route (self-verified; independent audit NOT_RUN). V111–V114 and V118–V120 remain NOT_RUN for the reasons above.
