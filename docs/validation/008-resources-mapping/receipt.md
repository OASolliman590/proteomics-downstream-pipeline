# R07 Maintainer-route receipt — Versioned annotation, protein groups and gene sets

Packet R07 (FR-061–FR-070 / T061–T070 / V061–V070) was implemented and verified on 2026-10-01 (all gates re-executed 2026-10-02 after the Phase 2/3 amendments) by the **Claude (Opus) route**, which the operator authorized to act for the Sol/Astra/Luna roles ([ADR 0007](../../adr/0007-claude-route-and-integration-amendments.md)). This is not work by the routes named in AGENTS.md.

## Identities

- Contract versions: kit 1.2.0 with the 2026-10-01 amendments (ADR 0006/0007); scientific methods contract v1.2.0.
- Working-source manifest (packet allowlist, 15 files): `ebeb8e7a008c80856a30ef94c32c45c83994153183f0ade504999e00c327d899` — lines in `evidence/working-source-manifest.txt` (`sha256  path`, sorted, LF).
- `verified_commit`: **null** — nothing is committed (branch `claude/full-pipeline`, uncommitted working tree).
- Environment: macOS-26.7-x86_64-i386-64bit-Mach-O; Python 3.13.15; R 4.6.1 with limma 3.68.5, statmod 1.5.2, impute 1.86.0, vegan 2.7-6, permute 0.9-10, DEqMS, proDA and fgsea at the versions pinned in `renv.lock`, in a project-local library (`.r-lib/`); `LANG=en_US.UTF-8`.

## Gates executed

| Command | Exit | Result | Log |
|---|---|---|---|
| `.venv/bin/python -m pytest tests/contract/test_resources.py -q -rs -p no:cacheprovider` | 0 | 11 passed in 34.98s | `python-tests-1.log` |
| `Rscript --vanilla -e res <- testthat::test_file('r/proteomicsCore/tests/testthat/test-mapping.R', reporter='summary', stop_on_failure=TRUE); df <- as.` | 0 | expectations=11 failed=0 skipped=0 errors=0 | `r-test-mapping.log` |

Phase-level regression (R01 suite, historical unittest, all R tests, baseline registry, spec checker, whole Python suite) is in [`docs/validation/phase-gates`](../phase-gates/README.md).

## Acceptance cases

| Case | Status | Tests | Positive oracle | Negative case |
|---|---|---|---|---|
| V061 | PASS | test_resources.py::test_v061_* | snapshot manifest and inner files verified against independent SHA-256 before analysis | byte mutation, missing file, all-zero placeholder, inner-file change all fail |
| V062 | PASS | test_resources.py::test_v062_* | 'proteomics resources prepare' builds a local snapshot; analysis runs with network blocked | no snapshot → E_RESOURCE_MISSING (no download); remote source refused |
| V063 | PASS | test_resources.py::test_v063_* | only the declared rat namespace maps; retired/unmapped kept with reasons | cross-species resource E_RESOURCE_TAXONOMY |
| V064 | PASS | test-mapping.R (V064); test_resources.py | projected collection labelled ortholog_projected with evidence; ambiguous orthologs excluded | projection without evidence refused |
| V065 | PASS | test-mapping.R (V065); test_resources.py | coverage → median → feature_id tie stages; invariant to column/label permutation | no statistic-based winner |
| V066 | PASS | test-mapping.R (V066); test_resources.py | multi-gene group excluded with reason; same-gene group retained | no duplication into both genes |
| V067 | PASS | test_resources.py::test_v063_to_v069_* | finite gene matrix with loss table; gene model t/P equal direct limma on that matrix | no imputation, no discarded-representative ranks |
| V068 | PASS | test-mapping.R (V068); test_resources.py | intersections 0/1/3/6 with bounds 2–5 → only the size-3 set eligible, before testing | — |
| V069 | PASS | test-mapping.R (V069) | ORA universe = six eligible measured mapped genes | resource genome never added |
| V070 | PASS | test_resources.py::test_v070_* | repeated runs byte-identical; group-label swap leaves representatives unchanged | corrupted snapshot fails E_RESOURCE_HASH (no refresh) |


## Notes

- Decisions taken on ambiguous points are D-01–D-30 in [decisions.md](../../../specs/001-downstream-proteomics/decisions.md).
- Independent audit (Astra-equivalent): **NOT_RUN**. Private-study regression: **NOT_RUN** (Maintainer-only, out of scope). No private data or private-derived numbers were used.

## Acceptance

All listed gates exited 0 and every acceptance case above, including its negative case, passed on re-execution. R07 is recorded as accepted by the operator-authorized Claude route (self-verified, uncommitted); independent review remains open.
