# R08 Maintainer-route receipt — Design-compatible pathways and enrichment

Packet R08 (FR-071–FR-080 / T071–T080 / V071–V080) was implemented and verified on 2026-10-01 (all gates re-executed 2026-10-02 after the Phase 2/3 amendments) by the **Claude (Opus) route**, which the operator authorized to act for the Sol/Astra/Luna roles ([ADR 0007](../../adr/0007-claude-route-and-integration-amendments.md)). This is not work by the routes named in AGENTS.md.

## Identities

- Contract versions: kit 1.2.0 with the 2026-10-01 amendments (ADR 0006/0007); scientific methods contract v1.2.0.
- Working-source manifest (packet allowlist, 11 files): `b5e3cdf146e19bb8e2f8f0a2527b1b81d5c81732487a102e84b4727410b54d15` — lines in `evidence/working-source-manifest.txt` (`sha256  path`, sorted, LF).
- `verified_commit`: **null** — nothing is committed (branch `claude/full-pipeline`, uncommitted working tree).
- Environment: macOS-26.7-x86_64-i386-64bit-Mach-O; Python 3.13.15; R 4.6.1 with limma 3.68.5, statmod 1.5.2, impute 1.86.0, vegan 2.7-6, permute 0.9-10, DEqMS, proDA and fgsea at the versions pinned in `renv.lock`, in a project-local library (`.r-lib/`); `LANG=en_US.UTF-8`.

## Gates executed

| Command | Exit | Result | Log |
|---|---|---|---|
| `.venv/bin/python -m pytest tests/scientific/test_pathways.py -q -rs -p no:cacheprovider` | 0 | 4 passed in 40.32s | `python-tests-1.log` |
| `Rscript --vanilla -e res <- testthat::test_file('r/proteomicsCore/tests/testthat/test-pathways.R', reporter='summary', stop_on_failure=TRUE); df <- as` | 0 | expectations=19 failed=0 skipped=0 errors=0 | `r-test-pathways.log` |

Phase-level regression (R01 suite, historical unittest, all R tests, baseline registry, spec checker, whole Python suite) is in [`docs/validation/phase-gates`](../phase-gates/README.md).

## Acceptance cases

| Case | Status | Tests | Positive oracle | Negative case |
|---|---|---|---|---|
| V071 | PASS | test_pathways.py::test_v071_v073_* | independent → CAMERA; fixed-subject/duplicate-correlation → ROAST | CAMERA on a paired design E_CAMERA_BLOCKED_DESIGN |
| V072 | PASS | test_pathways.py::test_v072_* | CAMERA P and estimated correlation equal direct limma::camera(inter.gene.cor = NA) | tiny set ineligible before testing |
| V073 | PASS | test_pathways.py::test_v071_v073_* | mroast with block/consensus correlation equals the direct seeded call (nrot 199, smoke label) | production nrot 199 E_ROTATION_RESOLUTION |
| V074 | PASS | test_pathways.py::test_v072_* | directional and mixed rows/families separate; central q recomputed per family | native per-call FDR kept as native only |
| V075 | PASS | test_pathways.py::test_v072_*; test-pathways.R | fgsea P/NES equal a direct seeded call on representative native t ranks; exploratory label | nonfinite/duplicate ranks E_FGSEA_RANK_INVALID |
| V076 | PASS | test-pathways.R (V076) | p = 1/45 and 1; BH q = 2/45 and 1; empty foreground p = 1 | dropping the zero-hit set changes q1 |
| V077 | PASS | test-pathways.R (V077); test_pathways.py | families never pool different nulls; BH/BY per family | — |
| V078 | PASS | test_pathways.py::test_v078_* | non-limma primary needs a named limma linear sensitivity; rows labelled sensitivity | E_PATHWAY_ENGINE_MISMATCH |
| V079 | PASS | test-pathways.R (V079) | leading-edge Jaccard/intersection displayed; P values unchanged | no merged sets |
| V080 | PASS | test-pathways.R (V080); test_pathways.py::test_v077_v080_*; R07 V062 | empty foreground keeps all eligible sets with p = 1; missing snapshot fails | all-nonfinite → E_PATHWAY_NO_FINITE_TESTS |


## Notes

- Decisions taken on ambiguous points are D-01–D-30 in [decisions.md](../../../specs/001-downstream-proteomics/decisions.md).
- Independent audit (Astra-equivalent): **NOT_RUN**. Private-study regression: **NOT_RUN** (Maintainer-only, out of scope). No private data or private-derived numbers were used.

## Acceptance

All listed gates exited 0 and every acceptance case above, including its negative case, passed on re-execution. R08 is recorded as accepted by the operator-authorized Claude route (self-verified, uncommitted); independent review remains open.
