# R14b Maintainer-route receipt — Post-DE robustness and sensitivity (operator-authorized amendment, ADR 0009)

Packet R14b (FR-137–FR-143 / T137–T143 / V137–V143) was implemented and verified on 2026-10-04 by the **Claude (Opus) route**, operator-authorized for the Sol/Astra/Luna roles ([ADR 0007](../../../adr/0007-claude-route-and-integration-amendments.md)). Independent audit: **NOT_RUN**.

## Identities

- Contracts: kit 1.3.0 (ADR 0009); scientific methods v1.3.0 (SM32, with SM03/SM04/SM07/SM13); `specs/015-post-de-analysis/contracts/post-de.md`.
- Working-source manifest (packet allowlist with the shared design factory, 7 files): `7bb09990e3ff76fa6037bf281c0508b85edef5de54c0322fdb8a80863c0bdcbf` — `evidence/working-source-manifest.txt`.
- Shared changes: A-2026-10-01-14 (re-plan API, used here) and A-2026-10-01-15 (required-module refusal hook in the workflow). Decision D-46.
- `verified_commit`: **null** (local commit on `claude/post-de`, not pushed). Environment: macOS x86_64, Python 3.13.15, R 4.6.1, limma 3.68.5, `.r-lib/`, `LANG=en_US.UTF-8`.

## Gates executed

| Command | Exit | Result | Log |
|---|---|---|---|
| `.venv/bin/python -m pytest tests/scientific/test_post_de_sensitivity.py -q -rs -p no:cacheprovider` | 0 | 10 passed | `evidence/python-tests-1.log` |
| `Rscript --vanilla -e "testthat::test_file('r/proteomicsCore/tests/testthat/test-post-de-sensitivity.R', ...)"` | 0 | 21 expectations, 0 failed, 0 skipped | `evidence/r-test-post-de-sensitivity.log` |
| `.venv/bin/python -m pytest tests -q -p no:cacheprovider -rs` | 0 | 297 passed, 0 skipped | `evidence/full-python-suite.log` |
| all R testthat files (`testthat::test_dir(..., load_package = "installed")`) | 0 | 524 expectations, 0 failed, 0 skipped | `evidence/r-all-testthat.log` |
| `.venv/bin/python scripts/check_spec_kit.py` | 0 | PASS | `evidence/spec-kit-check.log` |

## Acceptance cases

| Case | Status | Positive oracle | Negative case |
|---|---|---|---|
| V137 | PASS | imbalance cross-tabulation equals counts from the fixture; balanced/partial designs: additive and interaction estimable; confounded design (one group entirely male): R `qr()` rank shows additive estimable and interaction rank deficient; classes match, empty cell `C/F` named | a required interaction in the confounded design fails `E_SENSITIVITY_NONESTIMABLE` at plan time, before any fit (no `dea/`) |
| V138 | PASS | effects, moderated t, P and family q of the sex-adjusted model equal a direct `lmFit(~ group + sex)` + `eBayes(trend, robust)` + BH within the family; the model's families are `<family>__adj_sex` | the primary `dea/zero_null.tsv` is byte-identical to the hash recorded by the limma stage and holds no sensitivity rows |
| V139 | PASS | the male subgroup model equals a direct fit on the subset; coverage is re-applied on the subset | a group with one male is refused `E_SUBGROUP_COVERAGE` naming `A (1)`; no model directory or comparison row exists |
| V140 | PASS | retained/lost/gained sets, Jaccard, effect correlation and attenuation slope equal Python set arithmetic and R `cor()`/`lm(y ~ 0 + x)`; the non-interaction warning is present; the estimable interaction is its own family | no "significant in subgroup only" result is labelled an interaction effect (descriptive labels only) |
| V141 | PASS | 31-vs-11 design, 20 draws with seed 4242: discovery counts, drawn units and per-feature selection frequencies equal a direct seeded loop (`sample()` + limma + BH) | every draw lists whole biological units (B1 with two injections appears once); seed, draws and unit level are recorded; label `descriptive` |
| V142 | PASS | per-unit maximum absolute effect change and lost memberships equal leave-one-unit-out refits written in the test; the planted high-leverage unit B1 ranks first and its omitted source observations are both injections | one row per biological unit, never per injection |
| V143 | PASS | robustness fractions, sensitivity-model fractions and influence flags equal fractions computed in the test from the V138–V142 outputs; label `descriptive` | no adjusted-P column; primary q values unchanged |

## Acceptance

All gates exited 0 and every acceptance case above, including its negative case, passed. R14b is accepted by the operator-authorized Claude route (self-verified; independent audit NOT_RUN). Synthetic fixtures only.
