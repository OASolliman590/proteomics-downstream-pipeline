# R14d Maintainer-route receipt — Post-DE biomarker discrimination evaluation (operator-authorized amendment, ADR 0009)

Packet R14d (FR-149–FR-160 / T149–T160 / V149–V160) was implemented and verified on 2026-10-04 by the **Claude (Opus) route**, operator-authorized for the Sol/Astra/Luna roles ([ADR 0007](../../../adr/0007-claude-route-and-integration-amendments.md)). Independent audit: **NOT_RUN**.

## Identities

- Contracts: kit 1.3.0 (ADR 0009); scientific methods v1.3.0 (SM34–SM38 with SM23/SM25); `specs/015-post-de-analysis/contracts/post-de.md`.
- Working-source manifest (packet allowlist with the shared design factory, 7 files): `a5fa8c0c13d01a035f0477d832fca391502f79efa5e4f7e8b4b62b8d9849f9d8` — `evidence/working-source-manifest.txt`.
- Shared change: A-2026-10-01-17 (`classifier_params`; glmnet/e1071 API probes in the R11 lock check). Decision D-48.
- `verified_commit`: **null** (local commit on `claude/post-de`, not pushed). Environment: macOS x86_64 (8 cores), Python 3.13.15, R 4.6.1, glmnet 5.1, pROC 1.19.1, e1071 1.7-17, randomForest 4.7-1.2, limma 3.68.5, `.r-lib/`, `LANG=en_US.UTF-8`.

## Gates executed

| Command | Exit | Result | Log |
|---|---|---|---|
| `.venv/bin/python -m pytest tests/scientific/test_post_de_biomarker.py -q -rs -p no:cacheprovider` | 0 | 22 passed | `evidence/python-tests-1.log` |
| `Rscript --vanilla -e "testthat::test_file('r/proteomicsCore/tests/testthat/test-post-de-biomarker.R', ...)"` | 0 | 108 expectations, 0 failed, 0 skipped | `evidence/r-test-post-de-biomarker.log` |
| `.venv/bin/python -m pytest tests -q -p no:cacheprovider -rs` | 0 | 328 passed, 0 skipped | `evidence/full-python-suite.log` |
| all R testthat files (`testthat::test_dir(..., load_package = "installed")`) | 0 | 648 expectations, 0 failed, 0 skipped | `evidence/r-all-testthat.log` |
| `.venv/bin/python scripts/check_spec_kit.py` | 0 | PASS | `evidence/spec-kit-check.log` |

## Acceptance cases

| Case | Status | Positive oracle | Negative case |
|---|---|---|---|
| V149 | PASS | every R14d table carries exactly one label from the frozen vocabulary; no output contains "diagnostic", "clinical utility" or "validated" outside the label tokens and the mandated "not externally validated"; the report section carries its label | a file with a forbidden claim and a table without a label fail `E_BIOMARKER_CLAIM` |
| V150 | PASS | AUC, DeLong CI and stratified-bootstrap CI (seeded) equal `pROC::roc`/`ci.auc` for planted up, down and noise features; the prespecified-direction AUC of a down-shifted feature is below 0.5 | auto-direction AUCs are always labelled `in_sample` with the direction-optimism flag |
| V151 | PASS | 4 vs 30 units with 5-fold CV: the planner refuses (INAPPLICABLE) and records both `E_BIOMARKER_SMALL_N` and `E_BIOMARKER_FOLDS` | no AUC is emitted; a required module rejects the plan (exit 2) |
| V152 | PASS | 2000 pure-noise features, 20 per class: a deliberately leaky reference written in the test (selection and scaling on all data, then CV) has mean AUC > 0.75; the pipeline's per-repeat mean ± 1.96 SD includes 0.5 and the pooled AUC is within 0.2 of 0.5; the transform audit lists every transform per outer fold, fitted on training observations only | a transform fitted on held-out samples fails `E_BIOMARKER_LEAKAGE`; `imputation: global_median` is refused |
| V153 | PASS | out-of-fold predictions, per-repeat AUCs and their mean equal an independent nested loop written in the test (same documented fold dealing, seeds and selector, direct glmnet); the repeated variant (two observations per subject) never splits a subject and also equals the loop | `group_by_subject: false` in a repeated design and a fold that shares a subject fail `E_BIOMARKER_GROUP_LEAKAGE` |
| V154 | PASS | selection frequencies equal the per-fold selections of the independent loop; planted features are selected in ≥ 90 % of folds and noise features mostly < 20 % | the frequencies are not a single full-data selection |
| V155 | PASS | reversing the class level order leaves the pooled AUC and every score unchanged; R-level fits equal direct glmnet and e1071 (decision values) fits oriented by the model coding | an orientation flip (1 − AUC) is excluded; the SVM uses decision values; random forest is `NOT_RUN` (`E_ENGINE_NOT_AVAILABLE`) when the package is absent |
| V156 | PASS | B = 19 whole-procedure permutations: P = (k+1)/(B+1); signal at the floor 1/20; noise P > 0.1; subject-level permutations keep whole subjects; the compute guard refuses (`E_BIOMARKER_COMPUTE_GUARD`) instead of lowering B | `scope: final_classifier_only` fails `E_BIOMARKER_PERMUTATION_SCOPE` |
| V157 | PASS | a panel selected on the full noise data is `fixed_panel_cv`, flagged selection-optimistic (AUC > 0.7) and reported beside a nested reselection estimate near 0.5; an independent panel is `fixed_panel_cv` without the flag; panel hash recorded | a same-data panel is never reported without the flag and the nested comparison |
| V158 | PASS | Brier score, 10 equal-width calibration bins, sensitivity/specificity and their unit-level bootstrap CIs equal computations in the test from the out-of-fold predictions; thresholds come from the inner out-of-fold scores of each training fold (audit) | `threshold_rule: youden_test` fails `E_BIOMARKER_THRESHOLD_LEAKAGE` |
| V159 | PASS | disjoint cohort: the locked model (`locked_model.json`, SHA-256 recorded) applied by hand in the test gives the reported AUC and sensitivity; label `independently_validated`; cohort hashes recorded | a cohort sharing 3 subjects fails `E_SCORE_SELECTION_OVERLAP` (required: exit 2); `retune: true` fails `E_VALIDATION_RETUNED` |
| V160 | PASS | ROC, CV-repeat, permutation-null, calibration and selection-stability figures each have a source table; every plotted ROC point is re-derived from the out-of-fold predictions; the report shows claim labels, n per class, CV scheme, seeds and B | a figure without a source table would fail this test |

## Calibration

Null calibration of nested CV and the permutation P (R11 extension) is recorded with R14f (`docs/validation/015-post-de-analysis/R14f`).

## Acceptance

All gates exited 0 and every acceptance case above, including its negative case, passed. R14d is accepted by the operator-authorized Claude route (self-verified; independent audit NOT_RUN). Synthetic fixtures only; no diagnostic or clinical-use claim is made.
