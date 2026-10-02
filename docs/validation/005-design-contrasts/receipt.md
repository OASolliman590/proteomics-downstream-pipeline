# R04 Maintainer-route receipt — Design validation, exact contrasts and blocking

Packet R04 (FR-031–FR-040 / T031–T040 / V031–V040) was implemented and verified on 2026-10-01 (all gates re-executed 2026-10-02 after the Phase 2/3 amendments) by the **Claude (Opus) route**, which the operator authorized to act for the Sol/Astra/Luna roles ([ADR 0007](../../adr/0007-claude-route-and-integration-amendments.md)). This is not work by the routes named in AGENTS.md.

## Identities

- Contract versions: kit 1.2.0 with the 2026-10-01 amendments (ADR 0006/0007); scientific methods contract v1.2.0.
- Working-source manifest (packet allowlist, 9 files): `d1ab57e690e2e9c38a28705e688fd39c3f947e99ee59e3e41c05844c92bc4696` — lines in `evidence/working-source-manifest.txt` (`sha256  path`, sorted, LF).
- `verified_commit`: **null** — nothing is committed (branch `claude/full-pipeline`, uncommitted working tree).
- Environment: macOS-26.7-x86_64-i386-64bit-Mach-O; Python 3.13.15; R 4.6.1 with limma 3.68.5, statmod 1.5.2, impute 1.86.0, vegan 2.7-6, permute 0.9-10, DEqMS, proDA and fgsea at the versions pinned in `renv.lock`, in a project-local library (`.r-lib/`); `LANG=en_US.UTF-8`.

## Gates executed

| Command | Exit | Result | Log |
|---|---|---|---|
| `.venv/bin/python -m pytest tests/scientific/test_design_contract.py -q -rs -p no:cacheprovider` | 0 | 28 passed in 81.39s (0:01:21) | `python-tests-1.log` |
| `Rscript --vanilla -e res <- testthat::test_file('r/proteomicsCore/tests/testthat/test-design.R', reporter='summary', stop_on_failure=TRUE); df <- as.d` | 0 | expectations=27 failed=0 skipped=0 errors=0 | `r-test-design.log` |

Phase-level regression (R01 suite, historical unittest, all R tests, baseline registry, spec checker, whole Python suite) is in [`docs/validation/phase-gates`](../phase-gates/README.md).

## Acceptance cases

| Case | Status | Tests | Positive oracle | Negative case |
|---|---|---|---|---|
| V031 | PASS | test_design_contract.py::test_v031_* | coefficient map and numeric rows equal the explicit indicator/centred/product oracle; reversible percent-encoding | system()/source()/path traversal/undeclared terms E_DESIGN_TERM before R is invoked |
| V032 | PASS | test_design_contract.py::test_v032_* | reordering keeps aligned coefficients and reference meaning | missing reference E_DESIGN_LEVEL; nonfinite covariate E_DESIGN_COVARIATE_NONFINITE |
| V033 | PASS | test_design_contract.py::test_v033_* | rank/alias equal independent qr/svd; crossed design accepted | batch = treatment E_DESIGN_CONFOUNDED with aliases; no plan, no fit |
| V034 | PASS | test_design_contract.py::test_v034_* | fixed-subject coefficients; consensus correlation equals a direct duplicateCorrelation call; visits do not inflate n | missing subject IDs E_SUBJECT_MISSING; fixed+random encoding E_BLOCKING_CONFLICT |
| V035 | PASS | test_design_contract.py::test_v035_*; test-design.R | d=2, t=−1, r=1, r=d+t; interaction from coefficient weights | unknown coefficient E_CONTRAST_COEFFICIENT; reversed weights E_CONTRAST_DIRECTION |
| V036 | PASS | test_design_contract.py::test_v036_*; test-design.R | featurewise n/df/estimability and named reasons; excluded rows retained | no zero effect/P=1 for unestimable rows |
| V037 | PASS | test-design.R (V037 ×2) | contrast-as-coefficient SE/effect/sigma equal direct GLS with weights, a missing row and a blocked variant | contrasts.fit shortcut differs (>1e-6) on the sparse feature |
| V038 | PASS | test_design_contract.py::test_v038_*; test-design.R | DEqMS/proDA eligibility errors from the frozen table, independent of installation | absent adapter is NOT_RUN / E_CAPABILITY_NOT_IMPLEMENTED, never limma under its name |
| V039 | PASS | test_design_contract.py::test_v039_* | plan with masks/designs/contrasts/families before any fit; identical inputs give identical hash; weight/exclusion/reason/resource changes change it | fit with wrong plan hash or changed input E_PLAN_CHANGED; missing plan E_PLAN_REQUIRED |
| V040 | PASS | test_design_contract.py::test_v040_* | unequal n accepted; repeated-unit imbalance accepted; reorder-invariant estimability | singleton E_INSUFFICIENT_REPLICATION (QC-only still works); absent level E_DESIGN_GROUP_EMPTY |


## Notes

- Decisions taken on ambiguous points are D-01–D-30 in [decisions.md](../../../specs/001-downstream-proteomics/decisions.md).
- Independent audit (Astra-equivalent): **NOT_RUN**. Private-study regression: **NOT_RUN** (Maintainer-only, out of scope). No private data or private-derived numbers were used.

## Acceptance

All listed gates exited 0 and every acceptance case above, including its negative case, passed on re-execution. R04 is recorded as accepted by the operator-authorized Claude route (self-verified, uncommitted); independent review remains open.
