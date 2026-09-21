# Validation strategy — phased gates

Status: 1.2.0-frozen. All maintained implementation tests remain pending. Document/schema checks performed during hardening are not evidence that a statistical engine exists.

## Numeric tolerances

IDs, masks, counts, ranks, declared family membership, state and integer permutation counts require exact equality. Deterministic transforms, effects and BH on identical inputs require abs_error≤1e-10 plus a fieldwise relative tolerance 1e-8 where magnitude comparison is meaningful. Backend SE/statistic/P references use fieldwise abs≤1e-8 or rel≤1e-6, prespecified before candidate comparison; near-zero comparisons use the absolute condition rather than dividing by zero. Ill conditioning gets a diagnostic/rejection, not an unannounced larger tolerance. Tiny positive P values use documented log-scale comparison rather than silent underflow coercion.

Stochastic references fix seed, RNG kind, threads, package and input. Where exact repeated output is not guaranteed, prescribe a Monte Carlo interval/precision gate before candidate results. Figure source keys/counts are exact and numerical coordinates use their table's tolerance; image metadata timestamps are not scientific differences. Actual rendered figures still require visual inspection.

## Independent oracle layers

Use structural/identity/roundtrip tests, hand arithmetic (scale, contrast covariance, finite permutations, BH, hypergeometric, TOST), separate direct package calls, real Python→R→files→HTML integration, independent locked offline reproduction and explicitly scoped simulation. A reference that imports the production adapter under test is not independent. Private historical work remains Maintainer-only and is not required to author public synthetic tests.

## Phase 1 must-pass

All V001–V050 and V091–V094 are required for v0.1-limma-core. R10 remains partially complete, not completed as a whole.

| Scenario | Exact gate / oracle | Cases |
|---|---|---|
| No second log | log2 zero preserved; requested log2→log2 hard-fails before transformation. | V012, V021 |
| Anti-pseudoreplication | Four injections of one specimen yield n=1; whole subjects retained for paired/repeated inference. | V013–V014, V034, V049 |
| Confounded batch | Independent QR/SVD rejects batch=treatment with aliases; no pseudoinverse effect claim. | V031–V033, V040 |
| Exact SE | Sparse/nonorthogonal/weighted/blocked contrast matches direct coefficient refit and independent covariance; approximate shortcut is explicitly insufficient. | V037, V042 |
| Primary missingness/QC | No new primary imputation, original mask counts, display-only PCA fill and no silent sample removal. | V024–V030, V041 |
| Hypothesis/family separation | TREAT distinct from zero-null and display filtering; pooled primary BH, secondary auxiliary contrasts exported. | V044–V048 |
| Zero discoveries | Actual completed fit with zero rejections and thin report is COMPLETED/0; no “no biological effect” inference. | V050, V091–V094 |
| Shared-control covariance | For independent C/U/T covariance diag(.25,.25,.25), calculate Cov(U−C,T−U)=−.25 and Var(T−C)=.5, not 1. This is a Phase 1 contrast-algebra fixture, not an R09 score module. | V035, V050 |
| No score P values | Phase 1 config rejects score testing and produces no score-test/P-value artifact, including empty score-P tables. | V050, V092 |
| Environment/state failure | Actual missing R, child failure, collision and interrupted stage stay nonzero/unavailable/failed and render honestly. | V003–V010, V091–V094 |
| Requested capability requiredness | Every model declares required or optional independently of installed software; primary is required. Eligible missing optional engine remains in the plan and yields PARTIAL/nonzero, while the same required absence yields FAILED/nonzero. | V003, V006, V038, V059 |

R03/R04 enforce adapter eligibility boundaries and deferred capability states; they do not have to run nonexistent R06 engines to pass Phase 1. V027 is executable in R03 against R01 canonical hashing with a complete synthetic plan envelope, then is rerun against R04's actual AnalysisPlan under V039. V029 uses R03's narrow detection-only Fisher/BH implementation and is rerun after R05 to prove family separation. These integration reruns are mandatory Phase 1 barriers without creating backward implementation dependencies. All other science assigned to R02–R05 remains in Phase 1, including explicit vendor mappings, TMT bridge/no-bridge eligibility, detection sensitivity and qualified paired/repeated limma. No 1,000-dataset calibration requirement is imposed on Phase 1. Actual solved environment and direct numerical references are still necessary; it cannot claim R11 lock/calibration acceptance.

## Phase 2 must-pass

V051–V090 cover real DEqMS/proDA references and failures, hashed offline snapshots, species/orthology/multi-gene mapping, label-independent representatives, finite pathway refits, CAMERA/ROAST design dispatch, fgsea null labeling, ORA all-set multiplicity, covariance-aware response classes, TOST/conjunction eligibility and independent-score selection/randomization. Repeat Phase 1 gates for regression without changing their ownership.

The exact k=2 score fixture uses [1,2,3,4,5,6,7,8] with the top four versus bottom four: enumerate all 70 allocations. A distinct k=6 fixture uses [1,2,4,8,16,32,64,128] with observed group [4,32,64,128]; enumerate all 70 and compute the absolute mean-difference tail. Expected arithmetic is 2/70 and 6/70. These are independently authored synthetic oracle cases, not instructions to return those numbers on other data. Score selection-overlap negatives must physically lack P/q columns.

## Phase 3 must-pass

V095–V120, with compatibility reruns of V091–V094, complete full report/figures, actual environment lock restoration, semantic offline/cache behavior, independent reference matrix, calibration, cross-platform CI, performance and Maintainer-only legacy reconciliation. R12 public-release authorization/ownership/license and private validation remain explicit external gates; absence is never PASS.

## R11 release calibration protocol

Freeze simulator configuration, seed bank and expected methods before reading candidate outputs. Fast core all-null and mixture scenarios each use at least 1,000 independent datasets; a separately labeled CI smoke profile may be smaller. Future release generators may create the required larger datasets (for example 2,000 features and 4–10 biological units per group); no such matrix is included or fabricated in this hardening pass. Slow specialized profiles may use a prespecified smaller count only with their Monte Carlo precision/limitations made explicit; insufficient precision cannot support the same validation claim.

At alpha=.05 under the supported independent Gaussian all-null scenario, measure P(any rejection); require its one-sided 95% exact binomial upper confidence limit≤.075. For mixtures use mean dataset-level FDP, FDP=0 when no rejections; require the prespecified one-sided 95% Monte Carlo/bootstrap upper bound≤.075. Do not replace FDR by the fraction of null proteins falsely rejected. Report power, effect bias, nonestimability and failures, not only passing summaries.

For nominal 95% effect intervals, report stratified coverage and dataset-level uncertainty with a core aggregate target .93–.97 at adequate Monte Carlo precision. Proteins correlated within a simulated dataset are not independent replication units for the coverage interval. Include true effects below/on/above the TREAT margin and residual effects inside/on/outside both equivalence margins.

Heavy tails, unequal variances, MAR/MNAR dropout and other misspecifications are prespecified stress/applicability scenarios, not a universal nominal-error promise. Preserve their failures and restrict claims accordingly; do not move a failed core scenario into stress-only after results. Correlated pathway simulations define competitive and self-contained nulls separately. fgsea remains exploratory gene-set-null analysis, not a sample-label permutation test.

## R11 benchmark and report gates

Measure the 20,000-feature×100-observation/eight-contrast/2,000-set QC+limma+CAMERA+report workload on documented four-core hardware: targets≤8 GiB peak process-tree memory and≤30 minutes. Record actual hardware/time/memory and assess failures; do not extrapolate from a smaller matrix. Specialized engines and expensive release calibration are timed separately. Execute real Windows/Linux Python/R examples. Inspect actual full/null/failed reports and source-linked figures, not just file existence.

## R12 historical boundary

Maintainer alone verifies the original ZIP, 165-file preserved evidence, 3,714×20 matrix, 1,345 missing entries and Method A/B semantic crosswalk against actual private inputs. Compare six axis and nine joint contrasts where compatible runtime/resources permit. Negr1/A0A8I6AJV3's auxiliary treated-disease-versus-control result and old permutation fractions are regression/audit cases, never production constants or rescue proof. Missing historical R/resources yields NOT_RUN and, for numerical reproduction, NOT_REPRODUCED; the recovery did not itself rerun R. Changed families/filters/mappings/nulls legitimately require explicit reconciliation, not forced numeric matching.

## Evidence contract

For each acceptance ID record requirement/task ID, actual reviewed tree/commit, command, fixture hashes, independent oracle/expected values, observed differences, exit code, software/R session, artifact hashes, timestamp and reviewer conclusion. A passing process exit is not enough without inspected assertions. PASS, FAIL, NOT_RUN, INAPPLICABLE and SKIPPED remain distinct; none of the latter three is a passing software acceptance gate. Private logs stay local. Candidate traceability contains empty evidence arrays until execution.
