# R14f Maintainer-route receipt — Post-DE eligibility/dependency report and generality matrix (operator-authorized amendment, ADR 0009)

Packet R14f (FR-166–FR-167 / T166–T167 / V166–V167) was implemented and verified on 2026-10-04 by the **Claude (Opus) route**, operator-authorized for the Sol/Astra/Luna roles ([ADR 0007](../../../adr/0007-claude-route-and-integration-amendments.md)). Independent audit: **NOT_RUN**.

## Identities

- Contracts: kit 1.3.0 (ADR 0009); scientific methods v1.3.0 (SM41 with SM25); `specs/015-post-de-analysis/contracts/post-de.md`.
- Working-source manifest (packet allowlist with the shared design factory, 5 files): `6fd59ad6f46502ab3a1a25936742e74b4a5540387730812ab2db8b4cd0ba5bd9` — `evidence/working-source-manifest.txt`.
- Shared plumbing used as already amended: A-2026-10-01-14 (the `post_de_eligibility` capability, its workflow stage and report-section provider). Decision D-51. The generality matrix found an R14e defect, fixed separately before this packet (commit 5327df4, A-2026-10-01-19).
- `verified_commit`: **null** (local commit on `claude/post-de`, not pushed). Environment: macOS x86_64, Python 3.13.15, R 4.6.1, limma 3.68.5, glmnet 5.1, pROC 1.19.1, e1071 1.7-17, dynamicTreeCut 1.63.1, `.r-lib/`, `LANG=en_US.UTF-8`.

## Gates executed

| Command | Exit | Result | Log |
|---|---|---|---|
| `.venv/bin/python -m pytest tests/integration/test_post_de_generality.py -q -rs -p no:cacheprovider` | 0 | 12 passed (V166 + 11 V167 cells) | `evidence/python-tests-1.log` |
| `.venv/bin/python -m pytest tests -q -p no:cacheprovider -rs` | 0 | 351 passed, 0 skipped | `evidence/full-python-suite.log` |
| all R testthat files (`testthat::test_dir(..., load_package = "installed")`) | 0 | 677 expectations, 0 failed, 0 skipped | `evidence/r-all-testthat.log` |
| `.venv/bin/python scripts/check_spec_kit.py` | 0 | PASS | `evidence/spec-kit-check.log` |

## Acceptance cases

| Case | Status | Positive oracle | Negative case |
|---|---|---|---|
| V166 | PASS | One run where sets and sensitivity complete, association is refused for aliasing (`E_PHENOTYPE_ALIASED`), the random-forest biomarker is NOT_RUN because the package is shadowed as absent (`E_ENGINE_NOT_AVAILABLE`), and co-abundance is refused for small n (`E_COABUNDANCE_SMALL_N`). `post_de/eligibility/eligibility.json` and the report table list every requested module with plan eligibility, actual state, typed reason, claim label and verified inputs, equal to the stage states of the run. `dependency_summary.tsv` lists exactly the primary-family discoveries, and their robustness fractions equal `robustness_summary.tsv`. | Refused modules carry no claim label, no inputs and no values: the report shows "none (no values)", the dependency summary shows the module state (`INAPPLICABLE`, `NOT_RUN`), never zeros, and the report-data sections have `values: null`. NOT_RUN and INAPPLICABLE have different state meanings and CSS state classes. |
| V167 | PASS | 11 cells (two-group, three-group, paired, repeated, continuous exposure, 3 vs 3, 15 vs 80, 200 units, human/mouse/rat identifiers) × 5 modules, every module declared identically. Every module × cell reaches the expected state written in the test before the runs. Completed cells pass reduced module oracles: set membership equals the primary q < 0.05 set; robustness fractions are in [0, 1] for the primary members; association statistics equal a direct `limma` fit on group + centred phenotype (with the frozen consensus correlation in blocked designs, and the exposure term in the exposure cell); out-of-fold AUC and the permutation P recomputed from the predictions and the null; the planted co-abundance module recovered and tested; PPI P = (k+1)/(B+1) with the snapshot of each species convention. | Every stage reason is a typed `E_` code and none is `E_STAGE_EXCEPTION`; every claim label is from the frozen vocabulary; the declared classifier, selection, B, association method and module rule are the executed ones (no silent switch); refused modules write no stage result. Expected refusals: co-abundance in paired and repeated designs (`E_COABUNDANCE_DESIGN_UNSUPPORTED`) and 3 vs 3 (`E_COABUNDANCE_SMALL_N`); biomarker in 3 vs 3 (`E_BIOMARKER_SMALL_N`) and its paired single-feature AUC (`E_BIOMARKER_PAIRED_AUC`). |

## Findings during the packet

- The first matrix run crashed the networks stage in the two-group cell (one co-abundance module; limma drops the row names of a one-row fit). This was fixed in the R14e follow-up commit 5327df4 with a fail-before/pass-after test.
- The first draft of the matrix asserted association *power* in the repeated cell (12 subjects). Power is not an oracle, so it was replaced by a direct limma oracle on every cell. A power check is kept only for independent cells with at least 24 units.
- The `continuous_exposure` cell of the inherited draft was an ordinary two-group design. It now has the exposure `dose` as a centred design term and the primary contrast is the exposure slope.

## Calibration

The R11-style calibration extension for R14d (nested-CV AUC under the null, permutation-P uniformity, leaky-reference inflation) and for the R14e connectivity null is recorded in `calibration/` and its own section of this receipt after it runs (separate commit).

## Acceptance

All gates exited 0 and both acceptance cases, including their negative cases, passed. R14f is accepted by the operator-authorized Claude route (self-verified; independent audit NOT_RUN). Synthetic fixtures only.
