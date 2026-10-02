# R13 Maintainer-route receipt — Multivariate PERMANOVA and dispersion testing (operator-authorized amendment)

Packet R13 (FR-121–FR-130 / T121–T130 / V121–V130) was implemented and verified on 2026-10-01 (all gates re-executed 2026-10-02 after the Phase 2/3 amendments) by the **Claude (Opus) route**, which the operator authorized to act for the Sol/Astra/Luna roles ([ADR 0007](../../adr/0007-claude-route-and-integration-amendments.md)). This is not work by the routes named in AGENTS.md.

## Identities

- Contract versions: kit 1.2.0 with the 2026-10-01 amendments (ADR 0006/0007); scientific methods contract v1.2.0.
- Working-source manifest (packet allowlist, 8 files): `1ac3a27683ee55820426c6997f07410c856bec5a396005f5c4c4606a5e3fd574` — lines in `evidence/working-source-manifest.txt` (`sha256  path`, sorted, LF).
- `verified_commit`: **null** — nothing is committed (branch `claude/full-pipeline`, uncommitted working tree).
- Environment: macOS-26.7-x86_64-i386-64bit-Mach-O; Python 3.13.15; R 4.6.1 with limma 3.68.5, statmod 1.5.2, impute 1.86.0, vegan 2.7-6, permute 0.9-10, DEqMS, proDA and fgsea at the versions pinned in `renv.lock`, in a project-local library (`.r-lib/`); `LANG=en_US.UTF-8`.

## Gates executed

| Command | Exit | Result | Log |
|---|---|---|---|
| `.venv/bin/python -m pytest tests/scientific/test_permanova.py -q -rs -p no:cacheprovider` | 0 | 8 passed in 71.27s (0:01:11) | `python-tests-1.log` |
| `Rscript --vanilla -e res <- testthat::test_file('r/proteomicsCore/tests/testthat/test-permanova.R', reporter='summary', stop_on_failure=TRUE); df <- a` | 0 | expectations=54 failed=0 skipped=0 errors=0 | `r-test-permanova.log` |

Phase-level regression (R01 suite, historical unittest, all R tests, baseline registry, spec checker, whole Python suite) is in [`docs/validation/phase-gates`](../phase-gates/README.md).

## Acceptance cases

| Case | Status | Tests | Positive oracle | Negative case |
|---|---|---|---|---|
| V121 | PASS | test-permanova.R (V121); test_permanova.py::test_v121_* | z-scored Euclidean equals dist(scale(t(X))) and PCA-space distance; Manhattan equals vegdist; DEP-derived set equals the q≤.05 rows | unknown panel member E_PERMANOVA_FEATURE_SET; <2 usable features refused |
| V122 | PASS | test-permanova.R (V122); test_permanova.py::test_v122_* | Df/SS/R²/pseudo-F equal hand Gower sums of squares; P equals a seeded adonis2 | production <999 permutations E_PERMANOVA_RESOLUTION; P never 0 |
| V123 | PASS | test-permanova.R (V123) | global + 3 pairwise rows; pairwise P/R² equal direct subset adonis2; Holm and BH within the set | global row never adjusted |
| V124 | PASS | test-permanova.R (V124) | group term from the unrestricted margin model; covariate from within-group permutations; shuffleSet keeps groups | no group row from the blocked model |
| V125 | PASS | test-permanova.R (V125); test_permanova.py::test_v125_* | filled design interaction equals direct adonis2 margin model | empty cell E_PERMANOVA_INTERACTION_NONESTIMABLE (plan rejection if required, INAPPLICABLE if optional) |
| V126 | PASS | test-permanova.R (V126); test_permanova.py::test_v130_* | PERMDISP equals permutest; dispersion-only labelled dispersion_difference_location_not_established | location_shift never used when PERMDISP is significant |
| V127 | PASS | test-permanova.R (V127) | per-feature R² equal lm anova; mean equals multivariate R² (≤1e-8); not_applicable for Manhattan | violated identity E_PERMANOVA_IDENTITY |
| V128 | PASS | test-permanova.R (V128); test_permanova.py::test_v121_v128_* | seeded random null reproduces; best-possible = top-k mean; optimism flag | independent panel not flagged; same-data panels always carry context |
| V129 | PASS | test_permanova.py::test_v129_* | tables, PNG/PDF figures with source tables, report section with rule and flags, local links only | failed PERMANOVA shown FAILED with no values; run PARTIAL |
| V130 | PASS | test_permanova.py::test_v130_* | null P>.05; location P at the floor; dispersion-only PERMDISP P<.05; identical seeds give identical tables | — |


## Notes

- Decisions taken on ambiguous points are D-01–D-30 in [decisions.md](../../../specs/001-downstream-proteomics/decisions.md).
- Independent audit (Astra-equivalent): **NOT_RUN**. Private-study regression: **NOT_RUN** (Maintainer-only, out of scope). No private data or private-derived numbers were used.

## Acceptance

All listed gates exited 0 and every acceptance case above, including its negative case, passed on re-execution. R13 is recorded as accepted by the operator-authorized Claude route (self-verified, uncommitted); independent review remains open.
