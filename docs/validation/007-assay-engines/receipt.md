# R06 Maintainer-route receipt — Assay-qualified DEqMS and proDA backends

Packet R06 (FR-051–FR-060 / T051–T060 / V051–V060) was implemented and verified on 2026-10-01 (all gates re-executed 2026-10-02 after the Phase 2/3 amendments) by the **Claude (Opus) route**, which the operator authorized to act for the Sol/Astra/Luna roles ([ADR 0007](../../adr/0007-claude-route-and-integration-amendments.md)). This is not work by the routes named in AGENTS.md.

## Identities

- Contract versions: kit 1.2.0 with the 2026-10-01 amendments (ADR 0006/0007); scientific methods contract v1.2.0.
- Working-source manifest (packet allowlist, 11 files): `d7832f09939eb2471b538e0c751c9a0ee276c18f90fb68a6faca79e32b32893d` — lines in `evidence/working-source-manifest.txt` (`sha256  path`, sorted, LF).
- `verified_commit`: **null** — the packet was committed by Omar as `fd8dfa9`; the audit fixes and this re-verification are uncommitted on top of it (branch `claude/full-pipeline`).
- Environment: macOS-26.7-x86_64-i386-64bit-Mach-O; Python 3.13.15; R 4.6.1 with limma 3.68.5, statmod 1.5.2, impute 1.86.0, vegan 2.7-6, permute 0.9-10, DEqMS, proDA and fgsea at the versions pinned in `renv.lock`, in a project-local library (`.r-lib/`); `LANG=en_US.UTF-8`.

## Gates executed

| Command | Exit | Result | Log |
|---|---|---|---|
| `.venv/bin/python -m pytest tests/scientific/test_assay_engines.py -q -rs -p no:cacheprovider` | 0 | 10 passed in 255.29s (0:04:15) | `python-tests-1.log` |
| `Rscript --vanilla -e res <- testthat::test_file('r/proteomicsCore/tests/testthat/test-assay-engines.R', reporter='summary', stop_on_failure=TRUE); df ` | 0 | expectations=13 failed=0 skipped=0 errors=0 | `r-test-assay-engines.log` |

Phase-level regression (R01 suite, historical unittest, all R tests, baseline registry, spec checker, whole Python suite) is in [`docs/validation/phase-gates`](../phase-gates/README.md).

## Acceptance cases

| Case | Status | Tests | Positive oracle | Negative case |
|---|---|---|---|---|
| V051 | PASS | test-assay-engines.R (V051); test_assay_engines.py::test_v051_* | declared peptide counts accepted; zeros preserved; explicit justified offset recorded | proxies E_DEQMS_COUNT_EVIDENCE; zero without policy E_DEQMS_COUNT_NONPOSITIVE (no hidden +1) |
| V052 | PASS | test_assay_engines.py::test_v052_v053_* | sca.P.Value/sca.t equal a direct official lmFit→eBayes→spectraCounteBayes→outputResult call | ordinary limma P kept only as native_limma_* fields |
| V053 | PASS | test_assay_engines.py::test_v052_v053_* | shuffled count rows give identical results; sca and limma fields distinct; central q from sca P | duplicate count IDs fail (optional model FAILED, run PARTIAL) |
| V054 | PASS | test_assay_engines.py::test_v054_* | sparse exact route equals a direct reparameterized DEqMS reference | blocked/weighted/TREAT DEqMS rejected before execution (E_DEQMS_DESIGN_UNSUPPORTED) |
| V055 | PASS | test-assay-engines.R (V055); test_assay_engines.py::test_v055_* | only unimputed LFQ independent fixed designs eligible | TMT, unknown imputation, repeated blocks: named proDA errors, decided before any data stage; never limma fallback |
| V056 | PASS | test_assay_engines.py::test_v055_v056_v057_* | proDA/test_diff P, SE, df, n_obs match a seeded direct reference; all-missing-group feature keeps n_obs and dropout_prior_dependent | no invented zeros, no TREAT endpoint |
| V057 | PASS | test_assay_engines.py::test_v055_v056_v057_* | native SE/df exported with method-specific CI text | df not borrowed from limma |
| V058 | PASS | test_assay_engines.py::test_v058_* | stable-ID comparison discloses unmatched universes; primary limma table unchanged by alternative engines | alternative never promoted to primary |
| V059 | PASS | test_assay_engines.py::test_v059_* (4 cases) | missing package: optional → PARTIAL/3, required → FAILED/3 (E_ENGINE_NOT_AVAILABLE, kept in plan); child failure: optional → PARTIAL/4, required → FAILED/4 | no INAPPLICABLE relabel, no fallback table |
| V060 | PASS | test_assay_engines.py (all) | seeded fixtures with direct references; proDA convergence and seed recorded | genuine DEqMS failure kept as FAILED, not replaced with expected numbers; no calibration claim |


## Notes

- Decisions taken on ambiguous points are D-01–D-38 in [decisions.md](../../../specs/001-downstream-proteomics/decisions.md).
- Independent audit: an external audit of fd8dfa9 (2026-10-02) returned ACCEPT WITH FIXES; its findings are fixed (ADR 0008). Private-study regression: **NOT_RUN** (Maintainer-only, out of scope). No private data or private-derived numbers were used.

## Acceptance

All listed gates exited 0 and every acceptance case above, including its negative case, passed on re-execution. R06 is recorded as accepted by the operator-authorized Claude route (self-verified); the independent audit of fd8dfa9 returned ACCEPT WITH FIXES, and the fixes are self-verified and not yet re-audited.
