# Final acceptance record — local gate run of 2026-10-04

This records the last local re-run of every executable gate on branch `claude/post-de`, on the tree that contains R14a–R14f, the calibration extension and the R12 tooling (the gates ran on the working tree committed as bb3bb33; this record and the status documents were committed after it). It is **not** the V119 final integrated acceptance run: that run needs the final reviewed tree, CI on its pushed commit and the Maintainer-only gates V111 and V114. This is **not a v1.0 release**. No license has been selected.

Environment: macOS x86_64 (8 cores), Python 3.13.15, R 4.6.1, project-local `.r-lib/`, `LANG=en_US.UTF-8`. `verified_commit` stays null.

| Gate | Command | Exit | Result |
|---|---|---|---|
| Python suite | `.venv/bin/python -m pytest tests -q -p no:cacheprovider -rs` | 0 | 394 passed, 0 skipped (`013-release/evidence/full-python-suite.log`) |
| R suite | `Rscript --vanilla scripts/maintained/test_r.R --suite foundation` and `testthat::test_dir(...)` | 0 | 677 expectations, 0 failed, 0 skipped (`013-release/evidence/r-all-testthat.log`) |
| Spec kit | `.venv/bin/python scripts/check_spec_kit.py` | 0 | PASS |
| Release hygiene | `.venv/bin/python scripts/maintained/build_release.py check` | 0 | no blocking finding; license and ownership unresolved |
| R11 locks | `check_r_lock('renv.lock')`, `check_python_lock('requirements.lock')` | 0 | qualified |
| Post-DE calibration | `scripts/maintained/run_post_de_calibration.py --cores 8` (commit 3df8f9c) | 0 | all four gates PASS (`015-post-de-analysis/R14f/calibration`) |

Acceptance ledger after this run: 160 PASS, 7 NOT_RUN, 0 FAIL, 0 SKIPPED, 0 INAPPLICABLE (167 identities). The NOT_RUN cases are V111–V114 and V118–V120. Each has the reason recorded in [traceability](../../specs/001-downstream-proteomics/traceability.json): Maintainer-only private regression, Maintainer-approved summaries and disclosure review, operator decision required (AGENTS.md), final reviewed tree with CI, or Maintainer review.

External limits, reported separately and not counted as done:

- Independent audit of R14a–R14f and R12: NOT_RUN.
- Cross-platform CI for the post-DE commits: not run (nothing pushed). The last CI pass (run 37087074374) was on fa4d948.
- Maintainer-only gates V111/V114 and the decisions behind V113, V118 and V120: see the [Maintainer runbook](013-release/MAINTAINER_RUNBOOK.md).
