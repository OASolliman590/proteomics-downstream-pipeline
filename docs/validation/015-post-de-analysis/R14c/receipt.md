# R14c Maintainer-route receipt — Post-DE protein–phenotype association (operator-authorized amendment, ADR 0009)

Packet R14c (FR-144–FR-148 / T144–T148 / V144–V148) was implemented and verified on 2026-10-04 by the **Claude (Opus) route**, operator-authorized for the Sol/Astra/Luna roles ([ADR 0007](../../../adr/0007-claude-route-and-integration-amendments.md)). Independent audit: **NOT_RUN**.

## Identities

- Contracts: kit 1.3.0 (ADR 0009); scientific methods v1.3.0 (SM33, with SM03/SM25); `specs/015-post-de-analysis/contracts/post-de.md`.
- Working-source manifest (packet allowlist with the shared design factory, 7 files): `7996980e435ecee08d108f5f5035d85ecdb115522bcae5a77dbf51c53a5b189e` — `evidence/working-source-manifest.txt`.
- Shared change: A-2026-10-01-16 (`design_service.replan(group_term=False)`). Decision D-47.
- `verified_commit`: **null** (local commit on `claude/post-de`, not pushed). Environment: macOS x86_64, Python 3.13.15, R 4.6.1, limma 3.68.5, `.r-lib/`, `LANG=en_US.UTF-8`.

## Gates executed

| Command | Exit | Result | Log |
|---|---|---|---|
| `.venv/bin/python -m pytest tests/scientific/test_post_de_association.py -q -rs -p no:cacheprovider` | 0 | 9 passed | `evidence/python-tests-1.log` |
| `Rscript --vanilla -e "testthat::test_file('r/proteomicsCore/tests/testthat/test-post-de-association.R', ...)"` | 0 | 14 expectations, 0 failed, 0 skipped | `evidence/r-test-post-de-association.log` |
| `.venv/bin/python -m pytest tests -q -p no:cacheprovider -rs` | 0 | 306 passed, 0 skipped | `evidence/full-python-suite.log` |
| all R testthat files (`testthat::test_dir(..., load_package = "installed")`) | 0 | 539 expectations, 0 failed, 0 skipped | `evidence/r-all-testthat.log` |
| `.venv/bin/python scripts/check_spec_kit.py` | 0 | PASS | `evidence/spec-kit-check.log` |

## Acceptance cases

| Case | Status | Positive oracle | Negative case |
|---|---|---|---|
| V144 | PASS | phenotype coefficient, moderated t, P and family q equal a direct `lmFit(~ group + age + score)` + `eBayes(trend, robust)` + BH; one family `phenotype__score`; the 10 planted features are found | phenotype rows never enter the group-contrast family (`dea/` holds only primary-model rows) |
| V145 | PASS | n = 8 units: Spearman (via `cor.test`) and partial correlations (residuals of `lm(~ age)`) equal the oracle within 1e-12; permutation k and P = (k+1)/(B+1) equal a seeded `sample.int()` loop over units; BH within the family | permutation scheme is `biological_units`; R-level tests show subject-blocked permutations keep subjects whole or permute within subjects only |
| V146 | PASS | phenotype missing for 3 of 20 units: results equal a complete-case fit on 17 units; n = 17 and 3 missing are reported; the request holds exactly 17 phenotype values | `missing: mean_impute` fails `E_PHENOTYPE_IMPUTATION` |
| V147 | PASS | (a) phenotype recorded only in group B: the pooled analysis is refused, the within-B analysis runs and equals a direct fit; (b) Simpson fixture: pooled slope positive and within-group slopes negative (checked with `lm()`), both reported, features flagged and `W_PHENOTYPE_SIMPSON` warned | a required pooled association on the one-group phenotype fails `E_PHENOTYPE_ALIASED` before any fit |
| V148 | PASS | every heatmap cell equals the association table statistic; row orders are a permutation; the report shows family, n units and adjustment set | the heatmap has its source table (`heatmap_source.tsv`) |

## Acceptance

All gates exited 0 and every acceptance case above, including its negative case, passed. R14c is accepted by the operator-authorized Claude route (self-verified; independent audit NOT_RUN). Synthetic fixtures only.
