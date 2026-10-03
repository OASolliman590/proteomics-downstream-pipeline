# ADR 0008: Fixes for the independent audit of fd8dfa9

Status: accepted by the Maintainer route (Claude, operator-authorized; see ADR 0007). Date: 2026-10-02. `verified_commit`: null (uncommitted).

## Context

Omar committed the R02–R11/R13 work as `fd8dfa9` on `claude/full-pipeline`. An independent audit of that commit returned **ACCEPT WITH FIXES**, with two major and six minor findings. This is the first independent review of this work. The earlier receipts all recorded `independent_audit: NOT_RUN`.

## Decision

Each finding was fixed inside its owning packet. Every fix has a test that **fails on fd8dfa9 and passes after the fix**. Evidence: the new tests were run against an extracted fd8dfa9 tree with its own R package build (`docs/validation/audit-2026-10-02/evidence/before.log`) and against the fixed tree (`after.log`).

| Finding | Packet(s) | Fix |
|---|---|---|
| MAJOR 1: the primary limma model could be fitted on an imputed sensitivity matrix and still labelled primary | R01 (config, amendment A-2026-10-01-10), R05, R08, R09, R13 | Config rejects a sensitivity that targets the primary model (`E_SENSITIVITY_PRIMARY_MODEL`). `matrix_for_model` and an R guard refuse a primary fit on any non-observed matrix. Settings record the matrix and imputation actually used. Downstream builders call `guard_downstream` (`E_PRIMARY_NOT_OBSERVED`). |
| MAJOR 2: the subject-blocked PERMANOVA group test was degenerate (P = 1) | R13 | `pm_group_scheme` chooses one of three schemes. Whole subjects are permuted when group is constant within subject (balanced subjects required). Permutation is within subject when group varies within every subject. Mixed designs are refused. PERMDISP and pairwise tests use the same rule, and the scheme is recorded. Plan-time checks give the same refusals. |
| m1: calibration re-implemented the fit | R11 (and R03/R04 vectorization) | `calibrate_once` runs coverage_tables → featurewise_estimability → fit_limma_model → adjust_family. A higher-power mixture scenario was added. Calibration was re-run (results in the R11 receipt). |
| m2: weak oracles | R13 | Exact-label fixtures. The pipeline's own covariate and subject permutation designs are checked with `shuffleSet`. The either-label assertions were removed. |
| m3: silent pairwise and empty-group drops | R13 | `refusals.tsv`, warnings, `adjustment_family_size`/`adjustment_family_planned`, and a report table. |
| m4: an empty interaction cell refused the whole optional stage | R13 | Only the interaction term is refused. A required analysis is still rejected at plan time. |
| m5: contract and code disagreed | R13 | `not_in_primary_matrix` was removed from the contract. The R verifier checks the DEP table hash explicitly; it was already verified on read by `.pc_find_input`. |
| m6: unvalidated covariate names in the formula | R13 | Internal syntactic names; the declared name is reported. |

Further amendments: R03 `tests/integration/test_qc.py` V028 now targets a declared `role=sensitivity` model. The R13 V125 tests expect a term-level refusal for optional runs. The R13 V124/V126/V130 oracles were tightened.

## Addendum: first CI run (V107)

GitHub Actions run 36982402784 on fd8dfa9 failed: Ubuntu 2 tests (V059 assumed a local `.r-lib`) and Windows 16 tests. The Windows causes were argument loss with `Rscript -e`, cp1252 decoding of UTF-8 files and R output, the console-script lookup, and backslash paths split by POSIX `shlex`. Amendment A-2026-10-01-11 fixes these portably (D-40). V107 is recorded as FAIL until a CI run on the pushed fixes passes.

## Consequences

The affected packet gates, the full R testthat suite and the full Python suite were re-run (`docs/validation/phase-gates`). Acceptance stays self-verified apart from this audit, and nothing is committed. Decisions D-31 to D-38 record the details.
