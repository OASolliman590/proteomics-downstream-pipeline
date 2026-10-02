# R05 Maintainer-route receipt — Core limma inference and multiplicity

Packet R05 (FR-041–FR-050 / T041–T050 / V041–V050) was implemented and verified on 2026-10-01 (all gates re-executed 2026-10-02 after the Phase 2/3 amendments) by the **Claude (Opus) route**, which the operator authorized to act for the Sol/Astra/Luna roles ([ADR 0007](../../adr/0007-claude-route-and-integration-amendments.md)). This is not work by the routes named in AGENTS.md.

## Identities

- Contract versions: kit 1.2.0 with the 2026-10-01 amendments (ADR 0006/0007); scientific methods contract v1.2.0.
- Working-source manifest (packet allowlist, 9 files): `d0e31163d0549f1e0c702755f884f37f1b9a871eca580efbe62b2ddd14013ae7` — lines in `evidence/working-source-manifest.txt` (`sha256  path`, sorted, LF).
- `verified_commit`: **null** — nothing is committed (branch `claude/full-pipeline`, uncommitted working tree).
- Environment: macOS-26.7-x86_64-i386-64bit-Mach-O; Python 3.13.15; R 4.6.1 with limma 3.68.5, statmod 1.5.2, impute 1.86.0, vegan 2.7-6, permute 0.9-10, DEqMS, proDA and fgsea at the versions pinned in `renv.lock`, in a project-local library (`.r-lib/`); `LANG=en_US.UTF-8`.

## Gates executed

| Command | Exit | Result | Log |
|---|---|---|---|
| `.venv/bin/python -m pytest tests/scientific/test_limma.py -q -rs -p no:cacheprovider` | 0 | 9 passed in 48.74s | `python-tests-1.log` |
| `Rscript --vanilla -e res <- testthat::test_file('r/proteomicsCore/tests/testthat/test-limma.R', reporter='summary', stop_on_failure=TRUE); df <- as.da` | 0 | expectations=20 failed=0 skipped=0 errors=0 | `r-test-limma.log` |

Phase-level regression (R01 suite, historical unittest, all R tests, baseline registry, spec checker, whole Python suite) is in [`docs/validation/phase-gates`](../phase-gates/README.md).

## Acceptance cases

| Case | Status | Tests | Positive oracle | Negative case |
|---|---|---|---|---|
| V041 | PASS | test_limma.py::test_v041_* | effects/SE/t/P/df match a separate lmFit/eBayes reference on the declared universe; original NA retained | a filling adapter would differ (demonstrated) |
| V042 | PASS | test_limma.py::test_v042_* (unblocked and duplicate_correlation) | exported effect/SE/P match direct reparameterized limma; stdev.unscaled equals independent GLS; weights/block provenance recorded | settings record weights and block; no shortcut path |
| V043 | PASS | test_limma.py::test_v043_*; test-limma.R | prior/posterior variances match direct eBayes(trend, robust); settings saved; robust meaning stated | robust=TRUE without statmod fails preflight |
| V044 | PASS | test_limma.py::test_v043_v044_v047_* | TREAT P equals limma::treat(lfc=.5); separate type and family; display filter does not change P | TREAT table is complete, not a filtered zero-null table |
| V045 | PASS | test_limma.py::test_v041_v045_*; test-limma.R | CI = effect ± qt(.975, df_total)·SE | normal-quantile interval differs |
| V046 | PASS | test-limma.R (V046); test_limma.py::test_v048_* | pooled BH q = [.04,.04,2/3,.8]; secondary .001; BY harmonic factor; counts | per-contrast BH differs from the oracle |
| V047 | PASS | test_limma.py::test_v043_v044_v047_* | omnibus F/P equal limma eBayes F on the coefficient subset; interaction from weights | no interaction-significance flag from unequal significance |
| V048 | PASS | test_limma.py::test_v048_*; test-limma.R | all 24 planned rows; coverage exclusion with reason; numerical failure typed; family incomplete | missing rows, NA→0, orphan q rejected by table verification |
| V049 | PASS | test_limma.py::test_v049_* | whole-subject omission after technical aggregation; counts 2/2; effect change equals a direct refit | never single-injection or single-visit omission |
| V050 | PASS | test_limma.py::test_v050_* | zero discoveries COMPLETED/0; Cov(d,t) = −0.25, Var(r) = 0.5 | no score P artifact or column |


## Notes

- Decisions taken on ambiguous points are D-01–D-30 in [decisions.md](../../../specs/001-downstream-proteomics/decisions.md).
- Independent audit (Astra-equivalent): **NOT_RUN**. Private-study regression: **NOT_RUN** (Maintainer-only, out of scope). No private data or private-derived numbers were used.

## Acceptance

All listed gates exited 0 and every acceptance case above, including its negative case, passed on re-execution. R05 is recorded as accepted by the operator-authorized Claude route (self-verified, uncommitted); independent review remains open.
