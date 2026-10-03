# R09 Maintainer-route receipt — Treatment response, equivalence and independent scores

Packet R09 (FR-081–FR-090 / T081–T090 / V081–V090) was implemented and verified on 2026-10-01 (all gates re-executed 2026-10-02 after the Phase 2/3 amendments) by the **Claude (Opus) route**, which the operator authorized to act for the Sol/Astra/Luna roles ([ADR 0007](../../adr/0007-claude-route-and-integration-amendments.md)). This is not work by the routes named in AGENTS.md.

## Identities

- Contract versions: kit 1.2.0 with the 2026-10-01 amendments (ADR 0006/0007); scientific methods contract v1.2.0.
- Working-source manifest (packet allowlist, 14 files): `58ca1c832f7edde14f330a6637bb9ae96ab2092f8aeab7ea0641067292610eac` — lines in `evidence/working-source-manifest.txt` (`sha256  path`, sorted, LF).
- `verified_commit`: **null** — the packet was committed by Omar as `fd8dfa9`; the audit fixes and this re-verification are uncommitted on top of it (branch `claude/full-pipeline`).
- Environment: macOS-26.7-x86_64-i386-64bit-Mach-O; Python 3.13.15; R 4.6.1 with limma 3.68.5, statmod 1.5.2, impute 1.86.0, vegan 2.7-6, permute 0.9-10, DEqMS, proDA and fgsea at the versions pinned in `renv.lock`, in a project-local library (`.r-lib/`); `LANG=en_US.UTF-8`.

## Gates executed

| Command | Exit | Result | Log |
|---|---|---|---|
| `.venv/bin/python -m pytest tests/scientific/test_response.py -q -rs -p no:cacheprovider` | 0 | 7 passed in 170.22s (0:02:50) | `python-tests-1.log` |
| `Rscript --vanilla -e res <- testthat::test_file('r/proteomicsCore/tests/testthat/test-response.R', reporter='summary', stop_on_failure=TRUE); df <- as` | 0 | expectations=35 failed=0 skipped=0 errors=0 | `r-test-response.log` |

Phase-level regression (R01 suite, historical unittest, all R tests, baseline registry, spec checker, whole Python suite) is in [`docs/validation/phase-gates`](../phase-gates/README.md).

## Acceptance cases

| Case | Status | Tests | Positive oracle | Negative case |
|---|---|---|---|---|
| V081 | PASS | test_response.py::test_v081_* | r = d + t; unscaled Var(d) = .5, Cov(d,t) = −.25 per feature; shared-control warning | incoherent axis E_AXIS_INCOHERENT |
| V082 | PASS | test_response.py::test_v082_* (overlap) | overlapping training subjects → DescriptiveScore with no P/q/statistic columns | no IndependentScoreTest or rescue label |
| V083 | PASS | test-response.R (V083) | literal class boundaries incl. 1.2 vs 1.200001; RI = 3 is overshoot | no full_reversal class, no rounding |
| V084 | PASS | test-response.R (V084) | Fieller bounded set equals the quadratic-root oracle; unbounded/disjoint kept | no endpoint division or truncation |
| V085 | PASS | test-response.R (V085); test_response.py::test_v089_* | TOST p_lower/p_upper/max equal direct pt; 90% interval | nonsignificance alone is never equivalence |
| V086 | PASS | test-response.R (V086); test_response.py::test_v086_* | joint P = max of components; direction hash retained | overlapping direction cohort E_RESCUE_DIRECTION_NOT_INDEPENDENT |
| V087 | PASS | test-response.R; test_response.py::test_v082_* (independent) | fixed transform applied exactly; disjoint subjects → inferential | overlap/unknown provenance → descriptive |
| V088 | PASS | test-response.R (V088); test_response.py | exact k/70 = 2/70 and 6/70; end-to-end k equals itertools enumeration; MC (k+1)/(B+1) | 3/71 impossible by construction |
| V089 | PASS | test_response.py::test_v089_* | near_restoration that fails TOST stays not equivalent; tables separate | no equivalence wording in descriptive rows |
| V090 | PASS | test-response.R; test_response.py | shared-control negative covariance; exact denominators | in-sample scores never get P values |


## Independent audit 2026-10-02 (ADR 0008)

The response request builder now refuses primary limma results not fitted on observed data (`guard_downstream`, MAJOR 1). Before/after evidence: `docs/validation/audit-2026-10-02/`.

## Notes

- Decisions taken on ambiguous points are D-01–D-38 in [decisions.md](../../../specs/001-downstream-proteomics/decisions.md).
- Independent audit: an external audit of fd8dfa9 (2026-10-02) returned ACCEPT WITH FIXES; its findings are fixed (ADR 0008). Private-study regression: **NOT_RUN** (Maintainer-only, out of scope). No private data or private-derived numbers were used.

## Acceptance

All listed gates exited 0 and every acceptance case above, including its negative case, passed on re-execution. R09 is recorded as accepted by the operator-authorized Claude route (self-verified); the independent audit of fd8dfa9 returned ACCEPT WITH FIXES, and the fixes are self-verified and not yet re-audited.
