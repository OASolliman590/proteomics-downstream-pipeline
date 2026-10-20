# R15a receipt — Welch t-test engine (ADR 0011, optional engine under `methods_extension.welch`)

Packet R15a (FR-168 to FR-173 / T168 to T173 / V168 to V173) was implemented on 2026-10-13 by the **Claude route (Haiku builder, Opus review)**, as one stage of the operator-directed packet sequence. The packet is ready for review; it is not accepted. Independent audit: **NOT_RUN** (the reviewer checks this stage after it).

## Identities

- Contracts: kit 1.2.0 with the ADR 0011 identity range; scientific methods v1.4.0 (SM42, SM43); `specs/016-methods-extension-validation/contracts/methods-extension.md` (`methods_extension.welch`, `E_WELCH_*`, `W_WELCH_ADAPTED`).
- Working-source manifest (17 files: the R15a allowlist and the shared files changed under amendment A-2026-10-01-44): aggregate SHA-256 `d1b4e8ed5f78f5a6328946e6efc5180f0a120c35d7a9ba38daccc6776c01a4a9` — `evidence/working-source-manifest.txt`.
- `verified_commit`: **null** (local commit on `claude/methods-016`, not pushed).
- Environment: macOS (darwin 25.6.0), Python 3.13.15, R 4.6.1, limma 3.68.5, testthat 3.3.2, `.r-lib/` project library, `LANG=en_US.UTF-8`.

## What was built

- `r/proteomicsCore/R/model_welch.R`: per-feature two-sample Welch (Welch–Satterthwaite df, effect = mean(numerator) − mean(denominator), CI at the declared level) and the paired variant (within-subject differences, n − 1 df) on the observed log2 cells of the primary matrix (SM04, no new imputation). The stage writes the shared endpoint envelope with `engine = welch` and `result_type = WelchResult`, applies BH inside the declared family through the existing `adjust_family`, verifies the table with the existing `verify_differential_table`, and writes `welch_results.tsv`, `families.tsv`, `paper_compatible.tsv` (raw-P rows only, labelled `not_multiplicity_controlled`, BH q in the same table) and `eligibility.json` (decisions, adaptations, claim labels).
- `src/proteomics_pipeline/welch_service.py`: builds the stage request from the frozen plan (inputs hash-checked by R). It orients each pairwise contrast from its `group.<level>` weights, takes the Welch family from the plan's declared multiplicity families (E_FAMILY_MISSING otherwise), and records covariate adaptations.
- `src/proteomics_pipeline/reporting/welch_data.py`: reads only the stage's own outputs (claim labels, counts, families, adaptations).
- Shared, recorded under A-2026-10-01-44: the analysis schema (`welch` in `primary_engine` and `models[].engine`; optional `methods_extension`; Phase 5), run-status and report-data phase enums (phase 5), the spec copies of the analysis and run-status schemas (kept byte-identical), `config.py` (E_WELCH_NOT_DECLARED; declared role checked), `design_service.py` (welch eligibility; covariate adaptations written into the plan decision), `workflow.py` (welch capability in `methods_extension/welch/`, after the limma stage), `runtime.py` (capability map), `dispatch.R` and `capabilities.R` (handler and engine table), ADR 0002 (recorded amendment), `docs/methods/README.md` and the regenerated `docs/methods/reason-codes.md`.

## Decisions to review

1. **Phase 5.** The spec says Phase 1–4 configurations cannot enable `methods_extension`, but the runtime phase enum stopped at 4. Phase 5 is added to the enum, and a root rule refuses `methods_extension.welch.enabled` at phases 1–4.
2. **Secondary role mapping.** `methods_extension.welch.role` is `secondary` for a Welch model that is not primary. `models[].role` admits only `primary` and `sensitivity`, so such a model uses `role: sensitivity`. The config check enforces this mapping.
3. **Single primary.** A second primary model is refused by the schema's single-item rule (`E_CONFIG_SCHEMA`, "Too many items") before the semantic check runs.
4. **Family.** Welch endpoints must be members of a declared multiplicity family of the plan. A family that mixes Welch with another engine is refused (`E_WELCH_FAMILY_MIXED`).
5. **Paired variant.** `paired = auto` on a subject-blocked design uses the paired form when each subject has one observation per group. When the subjects do not pair, the stage runs the unadjusted two-sample form and records `block:<subject>` as an adaptation (`welch_block_not_pairable`). `paired = false` on a blocked design, or `paired = true` on an unblocked one, is refused (`E_WELCH_PAIRED_UNIT`).
6. **Non-pairwise contrasts.** A declared family member that cannot be oriented as a pairwise contrast is kept as a planned endpoint with state `not_run` and reason `E_WELCH_CONTRAST_NOT_PAIRWISE`, so the family stays complete.
7. **Units.** An unblocked Welch model refuses duplicated biological units (`E_WELCH_BIOLOGICAL_UNIT`, SM03).

## Gates executed

| Command | Exit | Result | Log |
|---|---|---|---|
| `python -m pytest tests/scientific/test_welch_engine.py -q -rs -p no:cacheprovider` | 0 | 13 passed, 0 skipped | `evidence/python-welch-tests.log` |
| testthat `test-model-welch.R` (`stop_on_failure = TRUE`) | 0 | 1222 expectations, 0 failed, 0 skipped | `evidence/r-model-welch-testthat.log` |
| `python -m pytest tests/regression tests/unit tests/scientific/reference/test_text_encoding.py -q -p no:cacheprovider` | 0 | 91 passed (R15a and R15b share this run) | `evidence/gates-regression-unit-encoding.log` |
| `python scripts/check_spec_kit.py` | 0 | PASS | `evidence/spec-kit-check.log` |

Full-suite results are recorded in the section "Full suites" below.

## Acceptance cases

| Case | Status | Positive oracle | Negative case |
|---|---|---|---|
| V168 | PASS | Every feature: t, Welch df, P, effect and CI equal `stats::t.test(var.equal = FALSE)` within 1e-10 (R level on 200 features; pipeline level on 200 features). | Pooled df (`var.equal = TRUE`, df = n − 2 = 13) differs; the effect sign follows the declared numerator (planted +2 in B is positive). |
| V169 | PASS | (a) a config with no `methods_extension` and the same config with the block disabled give byte-identical limma endpoint rows (run_id and plan_hash excluded); (b) a Welch secondary leaves the limma rows unchanged; (c) the Welch primary run contains no `dea/` directory and no limma row. | (d) two primary models are refused by the schema's single-item rule; (e) welch requested without the declaration fails `E_WELCH_NOT_DECLARED` at validation, so no Welch or relabelled limma result exists. |
| V170 | PASS | BH q equals `p.adjust(p, "BH")` within 1e-12 over the 200-feature family; the paper-compatible selection equals the raw-P set exactly; BH q sits in the same table. | Every paper-compatible row carries `not_multiplicity_controlled`; the q column differs from the raw P (it is not the raw P relabelled). |
| V171 | PASS | Paired t: t, df = 7, P, mean difference and CI equal `stats::t.test(paired = TRUE)` within 1e-10 for 8 subjects; `variant = paired`, `n_pairs = 8`. | `paired = false` on the declared paired design fails `E_WELCH_PAIRED_UNIT`; no `welch_results.tsv` is written. Unit level: the two-sample form on a declared pair is refused. |
| V172 | PASS | The covariate-adjusted Welch is not run: the statistic and effect equal the unadjusted `t.test` within 1e-10; the record `{welch, covariate:sex, adjusted, unadjusted, welch_not_covariate_adjusting}` appears in `eligibility.json`, in the plan's `engine_eligibility`, and as `W_WELCH_ADAPTED` in `warnings.json`. | The sex-adjusted group effect differs from the Welch effect on this unbalanced fixture (negative assertion on every feature). |
| V173 | PASS | Welch and limma effects are identical within 1e-10 (limma from the stage's own rows; limma P checked against `lmFit` + `eBayes(trend, robust)` in R); each engine's P equals its own oracle within 1e-10. | The maximum absolute difference between the Welch and limma P values exceeds 1e-6 (identical P values would indicate pooled variances). Each row keeps its own family (`welch_zero_null`, `primary-family`). |

End-to-end: a synthetic three-group study (control, disease, treated; 60 features; planted up and down effects) with `welch` as the primary engine runs to `COMPLETED` (plan, design, welch, report stages), `proteomics verify --run` returns `valid: true`, the Welch P values equal the `t.test` oracle within 1e-10, and `welch_report_data` reports the claim labels `welch` and `not_multiplicity_controlled`.

## Full suites

Final run on the reviewed tree (after the follow-up group F1–F4):

| Suite | Command | Exit | Result | Log |
|---|---|---|---|---|
| Full Python suite | `python -m pytest tests -q -rs -p no:cacheprovider` | 0 | 596 passed, 0 skipped (36 min) | `evidence/full-python-suite.log` |
| Full R suite | `testthat::test_dir("r/proteomicsCore/tests/testthat", load_package = "installed", stop_on_failure = TRUE)` | 0 | 2249 expectations, 0 failed, 0 skipped, 20 files | `evidence/r-all-testthat.log` |

The first full Python run of this stage (before F1 to F4) was stopped at about 24 percent: the post-DE generality test (`test_v167_*`) runs one full pipeline per cell and takes minutes per cell, so the run had not hung; it was stopped to start the final run. No result from that run is recorded.

## Follow-up group F1 (review of stage 1): Welch primary feeds the whole run

Recorded under A-2026-10-01-46. A welch primary now writes the `dea/` contract: `dea/zero_null.tsv` (the Welch rows with the shared envelope), `dea/families.tsv`, `dea/model_status.json` (engine welch, observed cells, no imputation), `welch_results.tsv`, `paper_compatible.tsv` and `eligibility.json`. The statistic label is `Welch t` (`paired t` for the paired variant). The report's dea section and the methods writers pick the primary differential stage (`dea/stage-result.json`), and describe a Welch model without limma moderation. The run catalogue draws the volcano, MA, dot-plot and Venn figures from the Welch `dea/` tables; `proteomics verify` returns exit 0. Post-DE association and eligibility consume the Welch families (the post-DE phase precheck now accepts Phase 5 as well as Phase 4). A welch primary cannot share a run with limma models (`E_PRIMARY_MODEL_CONFLICT`), because both would write `dea/`. The end-to-end tests are `test_f1_welch_primary_feeds_dea_report_figures_and_verify` and `test_f1_post_de_association_consumes_the_welch_dea_families`. The DEP heatmap is not asserted (it needs a discovery count that the synthetic study does not give).

## Gaps and doubts (for the reviewer)

- **Heatmap.** The DEP heatmap is not asserted in the Welch end-to-end test; the figure appears only when the family has discoveries in the synthetic study.
- **V169 (a) baseline.** The "identical to the R05 baseline" check compares the same configuration with and without a disabled `methods_extension` block, not an archived R05 artifact.
- **New reason codes beyond the spec's table:** `E_WELCH_BIOLOGICAL_UNIT`, `E_WELCH_CONTRAST_NOT_PAIRWISE`, `E_WELCH_FAMILY_MIXED`, `E_WELCH_WEIGHTS_UNSUPPORTED`, `E_WELCH_MATRIX_UNSUPPORTED`; the reason `insufficient_units_per_group` (not a code). All are documented in `specs/016-methods-extension-validation/contracts/methods-extension.md` and in the regenerated reason-code catalogue. The Maintainer should confirm them.
- **Role mapping (F4).** The spec's `secondary` is declared as `models[].role: sensitivity`, because the schema admits only `primary` and `sensitivity`. The mapping does not trigger `E_SENSITIVITY_PRIMARY_MODEL` (tested), and a sensitivity matrix on a welch model is refused (`E_WELCH_MATRIX_UNSUPPORTED`, tested).
- **Coverage and weights.** Welch accepts only `available_case` coverage and no precision weights (`E_COVERAGE_POLICY_UNSUPPORTED`, `E_WELCH_WEIGHTS_UNSUPPORTED`).
- **Phase label.** The stub report's phase label text ("Phase 1 thin offline report") is existing text and is not changed for phase 5 runs.
- **Maintainer ratification.** Amendments A-2026-10-01-44 and A-2026-10-01-46, the ADR 0002 amendment and the phase-5 enum are drafted by the implementing route and pending Maintainer confirmation.

## Acceptance

All packet gates exited 0 and every acceptance case above, including its negative case, passed. R15a is ready for review; it is not accepted until the reviewer and the Maintainer confirm the amendments.
