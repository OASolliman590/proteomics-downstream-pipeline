# R03 Maintainer-route receipt — Preprocessing, missingness and quality control

Packet R03 (FR-021–FR-030 / T021–T030 / V021–V030) was implemented and verified on 2026-10-01 (all gates re-executed 2026-10-02 after the Phase 2/3 amendments) by the **Claude (Opus) route**, which the operator authorized to act for the Sol/Astra/Luna roles ([ADR 0007](../../adr/0007-claude-route-and-integration-amendments.md)). This is not work by the routes named in AGENTS.md.

## Identities

- Contract versions: kit 1.2.0 with the 2026-10-01 amendments (ADR 0006/0007); scientific methods contract v1.2.0.
- Working-source manifest (packet allowlist, 16 files): `99d7adaee57c892e47e5508d77ef84a55db3d0af6a7e4d4f1e08574d348d9905` — lines in `evidence/working-source-manifest.txt` (`sha256  path`, sorted, LF).
- `verified_commit`: **null** — the packet was committed by Omar as `fd8dfa9`; the audit fixes and this re-verification are uncommitted on top of it (branch `claude/full-pipeline`).
- Environment: macOS-26.7-x86_64-i386-64bit-Mach-O; Python 3.13.15; R 4.6.1 with limma 3.68.5, statmod 1.5.2, impute 1.86.0, vegan 2.7-6, permute 0.9-10, DEqMS, proDA and fgsea at the versions pinned in `renv.lock`, in a project-local library (`.r-lib/`); `LANG=en_US.UTF-8`.

## Gates executed

| Command | Exit | Result | Log |
|---|---|---|---|
| `.venv/bin/python -m pytest tests/integration/test_qc.py -q -rs -p no:cacheprovider` | 0 | 8 passed in 6.09s | `python-tests-1.log` |
| `Rscript --vanilla -e res <- testthat::test_file('r/proteomicsCore/tests/testthat/test-preprocessing.R', reporter='summary', stop_on_failure=TRUE); df ` | 0 | expectations=44 failed=0 skipped=0 errors=0 | `r-test-preprocessing.log` |
| `Rscript --vanilla -e res <- testthat::test_file('r/proteomicsCore/tests/testthat/test-qc.R', reporter='summary', stop_on_failure=TRUE); df <- as.data.` | 0 | expectations=16 failed=0 skipped=0 errors=0 | `r-test-qc.log` |

Phase-level regression (R01 suite, historical unittest, all R tests, baseline registry, spec checker, whole Python suite) is in [`docs/validation/phase-gates`](../phase-gates/README.md).

## Acceptance cases

| Case | Status | Tests | Positive oracle | Negative case |
|---|---|---|---|---|
| V021 | PASS | test_qc.py::test_v021_* | preserve leaves values/masks equal; requested median is a distinct artifact with lineage | hidden normalization on a non-log2 scale E_NORMALIZATION_SCALE; display PCA input is a different artifact |
| V022 | PASS | test-preprocessing.R (V022 ×3); test_qc.py::test_v021_v022_* | median/reference factors match hand arithmetic; NA preserved; quantile is a complete-case sensitivity via limma | missing reference coverage E_REFERENCE_COVERAGE |
| V023 | PASS | test-preprocessing.R (V023) | within-plex loading factors and bridge adjustment match hand calculation; bridge-missing feature ineligible | missing bridge E_TMT_BRIDGE_REQUIRED; plex = treatment E_TMT_PLEX_CONFOUNDED |
| V024 | PASS | test-preprocessing.R (V024) | 4/4, 2/4 eligible; 1/4 excluded; 0/4 nonestimable; counts from the original mask; vectorized coverage_tables identical to the former loop (audit oracle test) | native_dropout deferred (not_run); unknown mask E_ORIGINAL_MASK_REQUIRED |
| V025 | PASS | test-preprocessing.R (V025) | observed / prior-imputed / numeric availability separated; no mechanism asserted | unknown mask stays unknown |
| V026 | PASS | test-qc.R (V026 ×4); test_design_contract.py::test_v026_negative_* | eigenvalues/scores match prcomp up to sign; one-NA display fill; pairwise n | constant/empty inputs INAPPLICABLE without variance; display matrix fed to limma E_DISPLAY_MATRIX_REJECTED |
| V027 | PASS | test_qc.py::test_v027_*; rerun on the real AnalysisPlan in test_design_contract.py::test_v039_* | extreme sample flagged not dropped; declared exclusion honoured; envelope hash changes with reason | blank reason E_EXCLUSION_REASON_REQUIRED |
| V028 | PASS | test-preprocessing.R (V028); test_qc.py::test_v028_* | MinDet = feature minimum; seeded Gaussian repeats; knn equals a direct impute.knn call; primary unchanged | MinDet with stochastic parameters E_SENSITIVITY_PARAMETERS; zero SD E_IMPUTATION_DISTRIBUTION |
| V029 | PASS | test-preprocessing.R (V029); post-R05 rerun test_limma.py::test_v029_* | Fisher P = 34/70 and 2/70 by enumeration; detection-only BH; separate family | paired design E_DETECTION_DESIGN_UNSUPPORTED (no P column); detection never enters abundance families |
| V030 | PASS | test_qc.py::test_v030_* | QC-only scope without model fit; every QC value has a source table | constant matrix gives INAPPLICABLE PCA, no fabricated panel |


## Independent audit 2026-10-02 (ADR 0008)

The V028 test now targets a declared role=sensitivity model, because sensitivities may no longer target the primary model. `coverage_tables` is vectorized; an oracle test keeps the former loop verbatim and checks identical output. Before/after evidence: `docs/validation/audit-2026-10-02/`.

## Notes

- Decisions taken on ambiguous points are D-01–D-38 in [decisions.md](../../../specs/001-downstream-proteomics/decisions.md).
- Independent audit: an external audit of fd8dfa9 (2026-10-02) returned ACCEPT WITH FIXES; its findings are fixed (ADR 0008). Private-study regression: **NOT_RUN** (Maintainer-only, out of scope). No private data or private-derived numbers were used.

## Acceptance

All listed gates exited 0 and every acceptance case above, including its negative case, passed on re-execution. R03 is recorded as accepted by the operator-authorized Claude route (self-verified); the independent audit of fd8dfa9 returned ACCEPT WITH FIXES, and the fixes are self-verified and not yet re-audited.
