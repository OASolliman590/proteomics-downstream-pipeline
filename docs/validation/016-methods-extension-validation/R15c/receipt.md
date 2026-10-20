# R15c receipt — Classifier calibration (FR-179 to FR-184; SM45; ADR 0011)

Packet R15c (FR-179 to FR-184 / T179 to T184 / V179 to V184) was implemented and run on 2026-10-13 by the **Claude route (Haiku builder, Opus review)**. The gate results are recorded **as produced**: two gates FAIL, and they were not changed. Independent audit: **NOT_RUN**.

## Order of the evidence (git)

1. `1c79fad` — the pre-registered gate file, its seeds and its replicate counts (`R15c/gates-preregistration.json`, SHA-256 `da7577a6cf5da61bd8b0a0d0ca0b9e3c47819cdf76be8451d6f0e06a34336a26`). Before this commit only one single-replicate timing pilot ran, in the scratch area; its outputs were not used for any gate.
2. `5a16524` — the calibration driver, the runner and the lasso and elastic-net classifiers (amendment A-2026-10-01-47), before any calibration run.
3. `44b7a71` — the calibration run outputs (`R15c/calibration/`).
4. `72281e4` — the acceptance tests (V179 to V184).

## Identities

- Contracts: scientific methods v1.4.0 (SM45); `specs/016-methods-extension-validation/contracts/methods-extension.md` (calibration gates; counts and seeds in the gate file).
- The runner refuses any gate file whose SHA-256 differs from the pinned value (`E_CALIBRATION_GATE_CHANGED`, tested).
- Working-source manifest (18 files: the driver, the runner, the tests, the R14d change, the schema, the README, the gate file and the calibration outputs): `evidence/working-source-manifest.txt`.
- `verified_commit`: **null** (local commits on `claude/methods-016`, not pushed).
- Environment: macOS (darwin 25.6.0), Python 3.13.15, R 4.6.1, glmnet 5.1, e1071 1.7.17, randomForest 4.7.1.2 (all present in `.r-lib`; nothing was installed), `LANG=en_US.UTF-8`, 12 worker processes.

## Design (as pre-registered)

- Families: `penalized_logistic` (ridge, glmnet alpha 0), `lasso` (alpha 1), `elastic_net` (alpha 0.5, declared; lambda tuned by inner 3-fold CV on the inner out-of-fold AUC), `svm_linear` and `svm_polynomial` (e1071, cost grid {0.1, 1, 10} tuned on inner AUC), `random_forest` (randomForest, 500 trees; the package is present, so the family ran).
- Nested procedure: the R14d production procedure (`pd_bm_nested`): repeated stratified 2 x 5-fold outer, 3-fold inner, top-5 AUC selection, Youden threshold on the training fold. Null: 15 vs 15 units, 200 pure-noise features. Planted signal: five features +1.5 SD in the positive class.
- Counts: nested AUC null 200 per family; permutation P 40 datasets x B = 19 per family; planted AUC 100 per family; planted permutation 40 x 19; leaky reference 60 (a panel chosen on all data before cross-validation); LOOCV (12 vs 12 units, inner 3-fold) 100 null and 50 planted, for ridge and lasso; subject-blocked (10 vs 10 subjects x 2 observations) 100 null AUC (grouped and ungrouped) and 30 x 19 permutations, for ridge, lasso and elastic net.
- Seeds: null bases 600000 (AUC), 700000 (permutation), 800000 (leaky); planted 900000 (AUC) and 950000 (permutation); blocked 1000000; LOOCV 1100000 (null) and 1150000 (planted); family offset 1000 x index. The ridge null seeds equal the R14f seeds.
- Compute: 1875 s on 12 cores (budget 45 min). No replicate count was reduced.

## Gate results (from `calibration/summary.json`)

| Family | Scenario | n | Statistic | Gate |
|---|---|---|---|---|
| `penalized_logistic` | nested AUC null | 200 | mean 0.4919, SE 0.0093, upper95 0.5102 | PASS |
| `penalized_logistic` | permutation P null (B = 19) | 40 | k(P<=0.05) 2, binomial P 0.601, chi-square P 0.657 | PASS |
| `penalized_logistic` | leaky reference | 60 | mean AUC 0.8918 | PASS |
| `penalized_logistic` | planted signal AUC | 100 | mean 0.9331 | PASS |
| `penalized_logistic` | planted signal P | 40 | share P<=0.05 0.950 | PASS |
| `penalized_logistic` | LOOCV null (12 vs 12) | 100 | mean 0.4332, SE 0.0242, upper95 0.4805 | **FAIL** |
| `penalized_logistic` | subject-blocked nested AUC null | 100 | mean 0.4770 | PASS |
| `penalized_logistic` | subject-blocked permutation P | 30 | binomial P 0.785, chi-square P 0.912 | PASS |
| `penalized_logistic` | ungrouped reference minus grouped | 100 | difference 0.1975 | PASS |
| `lasso` | nested AUC null | 200 | mean 0.4918, SE 0.0096, upper95 0.5106 | PASS |
| `lasso` | permutation P null (B = 19) | 40 | k(P<=0.05) 2, binomial P 0.601, chi-square P 0.590 | PASS |
| `lasso` | leaky reference | 60 | mean AUC 0.7891 | PASS |
| `lasso` | planted signal AUC | 100 | mean 0.9327 | PASS |
| `lasso` | planted signal P | 40 | share P<=0.05 0.900 | PASS |
| `lasso` | LOOCV null (12 vs 12) | 100 | mean 0.4233, SE 0.0216, upper95 0.4656 | **FAIL** |
| `lasso` | subject-blocked nested AUC null | 100 | mean 0.5005 | PASS |
| `lasso` | subject-blocked permutation P | 30 | binomial P 0.785, chi-square P 0.436 | PASS |
| `lasso` | ungrouped reference minus grouped | 100 | difference 0.1911 | PASS |
| `elastic_net` | nested AUC null | 200 | mean 0.4975, SE 0.0095, upper95 0.5160 | PASS |
| `elastic_net` | permutation P null (B = 19) | 40 | k(P<=0.05) 0, binomial P 1.000, chi-square P 0.522 | PASS |
| `elastic_net` | leaky reference | 60 | mean AUC 0.8676 | PASS |
| `elastic_net` | planted signal AUC | 100 | mean 0.9439 | PASS |
| `elastic_net` | planted signal P | 40 | share P<=0.05 0.975 | PASS |
| `elastic_net` | subject-blocked nested AUC null | 100 | mean 0.4920 | PASS |
| `elastic_net` | subject-blocked permutation P | 30 | binomial P 0.447, chi-square P 0.912 | PASS |
| `elastic_net` | ungrouped reference minus grouped | 100 | difference 0.1954 | PASS |
| `svm_linear` | nested AUC null | 200 | mean 0.5036, SE 0.0093, upper95 0.5217 | PASS |
| `svm_linear` | permutation P null (B = 19) | 40 | k(P<=0.05) 0, binomial P 1.000, chi-square P 0.395 | PASS |
| `svm_linear` | leaky reference | 60 | mean AUC 0.8541 | PASS |
| `svm_linear` | planted signal AUC | 100 | mean 0.9227 | PASS |
| `svm_linear` | planted signal P | 40 | share P<=0.05 1.000 | PASS |
| `svm_polynomial` | nested AUC null | 200 | mean 0.4914, SE 0.0071, upper95 0.5053 | PASS |
| `svm_polynomial` | permutation P null (B = 19) | 40 | k(P<=0.05) 2, binomial P 0.601, chi-square P 0.987 | PASS |
| `svm_polynomial` | leaky reference | 60 | mean AUC 0.8265 | PASS |
| `svm_polynomial` | planted signal AUC | 100 | mean 0.8858 | PASS |
| `svm_polynomial` | planted signal P | 40 | share P<=0.05 0.925 | PASS |
| `random_forest` | nested AUC null | 200 | mean 0.4909, SE 0.0088, upper95 0.5081 | PASS |
| `random_forest` | permutation P null (B = 19) | 40 | k(P<=0.05) 0, binomial P 1.000, chi-square P 0.0059 | **FAIL** |
| `random_forest` | leaky reference | 60 | mean AUC 0.8521 | PASS |
| `random_forest` | planted signal AUC | 100 | mean 0.9086 | PASS |
| `random_forest` | planted signal P | 40 | share P<=0.05 0.925 | PASS |

## Gates that FAIL (recorded, not changed)

- **`loocv_null` for `penalized_logistic` (mean 0.433) and `lasso` (mean 0.423)**. The gate is |mean - 0.5| <= 0.03 and mean + 1.96 SE <= 0.53. Both are below 0.5 by more than the pre-registered tolerance. This is the known pessimism of pooled leave-one-unit-out AUC on pure noise: a held-out unit's score is anti-correlated with the labels of the remaining units. The gate was not relaxed. V181 is FAIL.
- **`permutation_p_null` for `random_forest`**. The binomial P is 1.000 (no permutation P at or below 0.05 in 40 datasets), but the chi-square goodness-of-fit P over the 20 attainable values is 0.0059 (< 0.01), so the permutation P values are not uniform. V180 is FAIL; the SVM families pass their gates.

## Acceptance cases

| Case | Status | Evidence |
|---|---|---|
| V179 | PASS | Lasso, elastic-net and ridge scores equal a direct glmnet fit at a fixed lambda and panel within 1e-10 (negative: lasso scores differ from ridge by more than 1e-6). The null and planted gates PASS for the three families. |
| V180 | FAIL | SVM decision values equal e1071 within 1e-10 (negative: the score has more than two distinct values, not the predicted classes). Random forest runs. **Its permutation gate FAILS** (chi-square P 0.0059). |
| V181 | FAIL | The pooled LOOCV AUC equals a direct leave-one-unit-out loop within 1e-10 (lasso, fixed panel and lambda). The leaky reference is visibly inflated (negative). **The LOOCV null gate FAILS** for ridge and lasso. |
| V182 | PASS | No subject is split across folds in the subject-blocked design (negative: observation-level folds split subjects). The blocked null gates and the ungrouped reference PASS for ridge, lasso and elastic net. |
| V183 | PASS | The gate file is pinned (SHA-256 above); a changed copy is refused with `E_CALIBRATION_GATE_CHANGED`. |
| V184 | PASS | Seeds, counts, versions and timing are recorded (`calibration_raw.json`); one replicate reproduces exactly; a reduced count without a recorded adaptation is refused (`E_CALIBRATION_REPLICATES`, tested). |

## Deviations and doubts (for the reviewer)

1. **Planted-signal gates are proposed.** The gate file says they require Maintainer confirmation before the first run. Under the coordinator's instruction they were carried into this run; their confirmation is pending. Their PASS results are not evidence until confirmed.
2. **Classifier code (A-2026-10-01-47).** Lasso and elastic net were added to the R14d biomarker module through the glmnet path; the ridge branch is unchanged, and the R14d tests still pass.
3. **Runner addition after the run.** `verify_replicate_counts` (the `E_CALIBRATION_REPLICATES` rule) was added to `run_methods_calibration.py` after the run. It is a pure addition: the run does not call it, and it changes no gate or result.
4. **Metadata path sanitized.** `calibration_raw.json` recorded the absolute gate-file path; that one field was rewritten to the repository-relative path. No number changed.
5. **Exit codes.** The runner returns 1 when a gate fails. The driver itself exited 0, and `gate-results.json` records the driver's exit code (0); the gate verdicts come from `summary.json`.
6. **Scope.** LOOCV is calibrated for ridge and lasso only; the subject-blocked form for ridge, lasso and elastic net only; SVM radial is not calibrated; the blocked permutations use 30 datasets per family. These are the pre-registered counts.
7. **Synthetic only.** Noise and planted signal; no study data (SM45). A calibration pass is not a validation of a real study.

## Full suites

See the table at the bottom of this file.
