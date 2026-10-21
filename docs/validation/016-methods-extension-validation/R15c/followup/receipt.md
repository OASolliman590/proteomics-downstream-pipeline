# R15c follow-up receipt — leave-pair-out AUC and random-forest permutation diagnostic (D-77, D-78)

On 2026-10-21 the Maintainer decided how to follow up the two R15c gates that failed on 2026-10-13 (V180, V181). D-77: the loocv scheme reports the leave-pair-out AUC. D-78: run a larger, pre-registered random-forest permutation diagnostic. Both runs were done by the Claude route (Opus). The gate results are recorded **as produced**. The 2026-10-13 R15c summary and its FAIL verdicts are unchanged. Whether V180 and V181 are accepted on this evidence is the Maintainer's decision. Independent audit: **NOT_RUN**.

## Order of the evidence (git, branch `claude/r15c-followups`)

1. `85486bc` — decisions D-77 to D-80 and the follow-up gate file (`followup/gates-preregistration.json`, SHA-256 `e98bcb2d30b37605430ef905f2c34f08542371e44d108559b06b9f4ec80d50c4`), committed before the leave-pair-out code and before either run.
2. Random-forest diagnostic, started at `85486bc` from a frozen copy of `proteomicsCore` installed from that commit. There is no R change since the R15c run, and the copy is isolated from later reinstalls.
3. `53ac158` — the leave-pair-out code (A-2026-10-01-49), its tests and both drivers. The leave-pair-out calibration then ran from a frozen copy of `proteomicsCore` installed from that commit.
4. This commit — the run outputs, `calibration/summary.json` from `scripts/maintained/evaluate_calibration_followup.py` (which pins the gate file), and this receipt.

## Design (as pre-registered)

- **Leave-pair-out:** 12 vs 12 units and 200 pure-noise features, as in the R15c LOOCV scenario. Planted signal: five features +1.5 SD in the positive class. Every positive × negative unit pair (144) is held out together, and the whole inner procedure (top-5 AUC selection, inner 3-fold tuning) is refitted on the other units. The AUC is the mean pair comparison, with ties counting one half. Families: ridge and lasso. Counts: 100 null and 50 planted per family. Seeds: 1200000 (null) and 1250000 (planted), plus 1000 × the family index.
- **Random-forest diagnostic:** the R15c null design (15 vs 15 units, 2 × 5-fold repeated stratified outer, inner 3-fold, top-5 AUC, 500 trees) with the production whole-procedure permutation rule. 200 datasets with B = 99. Seeds: 1300000 + 5000 (family index 5) + dataset. For every dataset the run records the number of tied permutation AUCs (within 1e-12) and the mid-P.

## Results (`calibration/summary.json`)

| Scenario | Family | n | Statistic | Gate |
|---|---|---|---|---|
| `lpo_null` | ridge | 100 | mean 0.4913, SE 0.0196, upper95 0.5297 | PASS |
| `lpo_null` | lasso | 100 | mean 0.4953, SE 0.0218, upper95 0.5379 | **FAIL** (upper95 > 0.53) |
| `lpo_planted` | ridge | 50 | mean 0.9385 | PASS (proposed rule) |
| `lpo_planted` | lasso | 50 | mean 0.9267 | PASS (proposed rule) |
| `rf_permutation_diagnostic`, type I | random forest | 200 | 7 P ≤ 0.05 (rate 0.035), binomial P 0.876 | PASS |
| `rf_permutation_diagnostic`, uniformity (B + 1 values) | random forest | 200 | chi-square 88.0, df 99, P 0.778 | PASS |
| `rf_permutation_diagnostic`, uniformity (10 bins) | random forest | 200 | counts 19, 11, 18, 22, 21, 19, 23, 17, 24, 26; chi-square 8.1, df 9, P 0.524 | PASS |
| ties (descriptive) | random forest | 200 | 24 datasets with a tie; mean 0.135 and at most 2 tied permutations; mean P 0.544, mean mid-P 0.538 | — |

## Reading

- **V181.** The pooled leave-one-out AUC was biased low under the null (0.433 ridge, 0.423 lasso; 2026-10-13). The leave-pair-out AUC removes that bias (0.491 ridge, 0.495 lasso). Ridge meets the pre-registered gate. Lasso's mean is within 0.005 of 0.5, but its Monte Carlo upper bound (0.538) exceeds 0.53. Its replicate SD is about 0.22, so with 100 replicates the upper bound sits at about 0.538 even for an unbiased mean. The gate is recorded as FAIL; it was not changed.
- **V180.** With 200 datasets and B = 99, the random-forest permutation P meets every pre-registered rule. Tied permutation AUCs are rare, so ties do not explain the 2026-10-13 pattern. The earlier FAIL (40 datasets, B = 19, 2 expected per attainable value, chi-square P 0.0059) is consistent with chance in an under-powered check.

## Reproducibility

The leave-pair-out driver first ran to completion and wrote both tables, then failed while writing `lpo_raw.json` (an openssl hash object is not JSON-serialisable). The driver was fixed and the run repeated from the same frozen library. Both tables were byte-identical to the first run (`evidence/lpo-run.log`). The random-forest driver failed at the same step after writing its table. Its `rf_diagnostic_raw.json` was written afterwards from the recorded versions and states that the elapsed time was not recorded.

## Commands

| Command | Exit | Log |
|---|---|---|
| `Rscript --vanilla tests/scientific/calibration_rf_diagnostic.R <gates> <out> 12` | 1 (after writing the table; provenance step, see above) | `evidence/rf-diagnostic-run.log` |
| `Rscript --vanilla tests/scientific/calibration_lpo.R <gates> <out> 4` (second run) | 0 | `evidence/lpo-run.log` |
| `python scripts/maintained/evaluate_calibration_followup.py` | 0 | `calibration/summary.json` |
| `python -m pytest tests/scientific/methods_extension/test_methods_extension_calibration.py -k v181` before the code | 1 (2 failed) | `evidence/v181-fail-before.log` |
| the same after the code | 0 (3 passed) | `evidence/v181-pass-after.log` |

## Open for the Maintainer

**Resolved 2026-10-21:** item 1 — V180 accepted (D-81); item 2 — 400 fresh lasso replicates run (D-82, `R15c/followup-lasso400/`), still FAIL; item 3 remains open.

1. Accept V180 on the diagnostic evidence, or keep the 2026-10-13 FAIL.
2. V181 for lasso: accept on the unbiased mean, pre-register a larger leave-pair-out run (for example 400 replicates, which narrows the upper bound to about 0.52 if the mean holds), or keep FAIL.
3. Confirm or revise the proposed planted-signal rules (carried unchanged from R15c).
