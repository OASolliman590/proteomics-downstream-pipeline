# R15c follow-up 2 receipt — 400 lasso leave-pair-out null replicates (D-82)

On 2026-10-21 the Maintainer asked for 400 lasso replicates of the leave-pair-out null (D-82). The first follow-up had an unbiased mean but failed the gate's upper bound with 100 replicates. The run was done by the Claude route (Opus). The result is recorded **as produced**, and V181 follows it. Independent audit: **NOT_RUN**.

## Order of the evidence (git, branch `claude/r15c-followups`)

1. `595572c` — decisions D-81 to D-83 and this run's gate file (`gates-preregistration.json`, SHA-256 `e8c3c0dac5e03f60d930fd2078e6bb852ddb7331e17aa3142a79319976b25901`), committed before the run.
2. The run, on 12 cores, from the frozen copy of `proteomicsCore` installed from `53ac158` (the leave-pair-out code). The biomarker code is unchanged since that commit. It took 1085 s.
3. This commit — the outputs, `calibration/summary.json` (from `scripts/maintained/evaluate_calibration_followup.py --gates`), the logs and this receipt.

## Design (as pre-registered)

The design is the leave-pair-out null of `R15c/followup`: 12 vs 12 units, 200 noise features, every positive × negative pair held out, inner 3-fold tuning, and top-5 AUC selection. It ran for lasso only, with 400 replicates and seeds 1401001 to 1401400, which are disjoint from every earlier calibration. The gate is unchanged: |mean − 0.5| ≤ 0.03 and mean + 1.96 SE ≤ 0.53. The 100 earlier replicates are not pooled with these.

## Result

| Run | n | mean | SE | upper95 | median | SD | Gate |
|---|---|---|---|---|---|---|---|
| this run | 400 | 0.5204 | 0.0101 | 0.5402 | 0.528 | 0.203 | **FAIL** |
| `R15c/followup` (for reference) | 100 | 0.4953 | 0.0218 | 0.5379 | 0.469 | 0.218 | FAIL |
| `R15c/followup` ridge (for reference) | 100 | 0.4913 | 0.0196 | 0.5297 | 0.493 | 0.196 | PASS |

## Reading

The lasso leave-pair-out AUC under the null is not biased low, as the pooled leave-one-out AUC was (0.423 on 2026-10-13). With 400 replicates, however, its mean is 0.520, about two Monte Carlo SE above 0.5, and the upper bound (0.540) exceeds 0.53. The gate fails, and V181 stays FAIL. A slight optimism of about 0.02 for lasso in this small design (12 vs 12 units, 200 features) cannot be excluded on this evidence. No mechanism has been identified, and the distribution shows no degenerate pattern: only 5 of 400 replicates equal exactly 0.5. Ridge meets the gate.

## Provenance note

At commit `595572c` the driver wrote a fixed gate path into `lpo_raw.json`. The path was corrected after the run, and the record says so. The recorded `gates_sha256` was computed from the file the run read, and it equals this pre-registration. The driver now records the path it is given.

## Open for the Maintainer

- V181 for lasso stays FAIL. Possible next steps, each to be pre-registered:
  - investigate lasso's inner tuning in very small designs;
  - restrict the loocv scheme to ridge (and the other families that pass);
  - report lasso leave-pair-out AUCs with a documented calibration caveat.
