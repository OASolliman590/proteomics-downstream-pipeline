# Validation strategy and acceptance thresholds

Validation is designed before implementation. Correctness, calibration, reproducibility, portability and scientific provenance are separate outcomes. A passing test suite is not proof the old study's biological conclusions are true.

## Layers and independent oracles

1. **Contract tests:** schemas/IDs/masks/grain and roundtrip behavior. Reference values are hand-calculated fixtures or fixed files authored independently from implementation.
2. **Mathematical tests:** means, contrast algebra, exact finite permutations, hypergeometric tails, BH/BY and TOST on small examples. Compare against independent Python/SciPy/R references as appropriate; a wrapper invoking itself is not an oracle.
3. **Backend golden tests:** call pinned limma/TREAT, DEqMS, proDA, CAMERA, mroast and fgsea directly in minimal reference scripts separate from the production adapter. Save input, output, version and purpose. Do not require unrelated engines to agree statistically.
4. **Scientific simulations:** assess error, effect/interval behavior and applicability across declared data-generating assumptions. Preserve expected limitations for misspecified models instead of tuning implementation to conceal them.
5. **Integration:** run canonical examples through the real Python→R→files→HTML path, including no discoveries, no eligible pathways and failed stages.
6. **Historical regression:** Maintainer locally runs the archive conversion/comparison and maintained analysis, without giving private data to Packet Implementer, Independent Reviewer or any documentation worker. Preserve historical ambiguities and exact expected intentional changes.
7. **Clean-environment reproduction and final visual inspection:** independently restore locks, execute offline cached resources, inspect actual report/figures and compare semantic outputs.

## Numeric tolerances

- ID alignment, masks, group counts, rank/contrast structure, output statuses, selected family membership and exact permutation counts: exact equality.
- Deterministic scalar transforms, logFC and BH on the same finite input: absolute ≤1e-10 and relative ≤1e-8 where appropriate. Underflow P values compare on log scale with a documented finite floor; do not silently turn tiny positive P into zero.
- Backend reference SE/statistic/P under identical versions and input: absolute ≤1e-8 or relative ≤1e-6, stated fieldwise. Ill-conditioned cases require explicit conditioning diagnostics, not blanket widened tolerances.
- Stochastic enrichment/rotation: seed/RNG/backend/thread/package settings must be fixed. Compare native deterministic output where the engine guarantees it; otherwise prespecify repeated-run Monte Carlo tolerances/intervals and retain numerical diagnostics. Never widen tolerance after a failure without diagnosis and an ADR.
- Figure data: exact keys/counts and numeric table tolerance. PDF/PNG metadata/timestamps are not scientific differences; rendered layout must still be inspected.

## Mandatory scenario matrix

| Scenario | Expected behavior |
|---|---|
| Complete independent Gaussian log2 protein data | limma estimates and uncertainty match direct reference; null calibration assessed |
| Unequal n + continuous covariate + batch | Exact identifiable contrast matches direct fit; coefficient meanings preserved |
| Missing values + nonorthogonal design + weights | Exact refit/covariance path required; demonstrate old approximate shortcut differs on a constructed case |
| Fixed paired subject and duplicateCorrelation repeated model | Correct design/block handling; whole biological-unit influence; CAMERA incompatibility tested |
| Confounded batch/treatment | Plan rejects unidentifiable estimand, no correction claiming recovery |
| Technical replicates per specimen | Known aggregation and biological n, no pseudoreplication |
| Upstream imputed versus true observed input | Masks and coverage language differ; ineligible dropout fit rejected |
| Qualified LFQ dropout including a group with no observations for some features | proDA reference behavior, observed counts and prior-driven uncertainty preserved; limma nonestimability is distinct |
| Count-dependent variance with peptide/PSM evidence | DEqMS uses strictly positive aligned count covariate and sca statistics; missing/proxy counts rejected |
| TMT with bridge and balanced no-bridge design | Explicit intended normalization/plex model, confounding test and no fabricated bridge |
| Correlated genes within sets | CAMERA/ROAST evaluated under their actual nulls; fgsea remains appropriately exploratory |
| Duplicate/multi-gene mappings and orthology | Label-independent representative rule, mapping loss and finite universe correctly reported |
| Zero significant proteins/pathways | Valid complete analysis and honest report with zero discoveries |
| Shared untreated group in disease/treatment contrasts | Negative noise covariance demonstrable; no current-selected confirmatory score route |
| RI near 0/1 and >2; disease denominator near zero | Descriptive categories/eligibility and uncertainty handle boundaries without full-rescue claims |
| Equivalence versus nonsignificance | A broad interval overlapping zero fails equivalence; narrow interval inside margin passes correct TOST |
| Four vs four exact score test | 70 assignments; cases with k=2 and k=6 yield 2/70 and 6/70, not 3/71 or 7/71 |
| Paired/restricted score randomization | Exact admissible permutations/sign flips only; independence evidence required |
| Resource mismatch, dependency missing, R crash, killed stage | Nonzero/honest state, preserved logs, no stale cache accepted and no report PASS fiction |

## Calibration protocol

Freeze simulator configuration, seed bank and expected methods before reading candidate results. Use independent seeds from the unit/golden fixtures. Release profile uses at least 1,000 independently generated datasets for fast core all-null and mixture scenarios (e.g., 2,000 features, 4–10 biological units per group), with a separately recorded smaller smoke profile for CI. Specialized slow engine profiles may use a prespecified smaller number only when their Monte Carlo intervals are correspondingly reported; insufficient precision cannot support the same validation claim.

Under an eligible all-null independent Gaussian scenario at declared alpha=.05, measure P(any rejection), because all discoveries are false. Require its one-sided 95% binomial upper confidence limit ≤0.075 for the release core profile. For prespecified mixtures, measure the mean false-discovery proportion across simulations, setting FDP=0 when there are no rejections; require a prespecified bootstrap/Monte Carlo one-sided 95% upper bound ≤0.075 under the supported model assumptions. Report achieved power, bias, uncertainty coverage, nonestimability and failures as well. Do not use empirical false-positive rate among proteins as a substitute for FDR.

For nominal 95% effect CIs, report coverage by effect/abundance/missingness strata and require aggregate coverage within the prespecified Monte Carlo acceptance interval around 0.95 (release core target 0.93–0.97 with enough independent simulation units). Correlated features do not create independent biological replicates for CI coverage uncertainty; estimate uncertainty by simulated dataset/block. Effect-threshold tests include true effects at 0, below, at and above the margin. Equivalence simulations include values inside, on and outside both margins.

Additional stress scenarios include heavy-tailed outliers, unequal variances, MAR/MNAR dropout and correlation. These are applicability/robustness diagnostics; an engine need not be guaranteed under every misspecified scenario. Document where validation does not support nominal error and restrict the claimed capability or default accordingly. Do not silently relabel a failed core assumption test “stress-only” after observing it.

Pathway calibration uses constructed membership and correlation structures with stated competitive versus self-contained nulls. ROAST and CAMERA do not test interchangeable hypotheses. Gene-set permutation limitations are retained for fgsea rather than demanding or advertising sample-randomization guarantees it does not provide.

## Historical validation boundaries

Verify the original ZIP hash, all 165 preserved files, matrix 3,714×20, 1,345 missing entries, Method A/B input equivalence and source contrast identities. Legacy-mode expected core comparisons include six axis and nine joint contrast effects/BH, the Negr1 auxiliary joint contrast exception and exact fixed-score arithmetic. These are fixture expectations in regression, never production constants. New primary families/filtering/models/pathway nulls may change results; compare with a table explaining each expected change. Historical moderated P values and fgsea exact values cannot be declared reproduced without an executed compatible environment/resource version.

## Evidence and completion

Each V001…V120 artifact records acceptance ID, requirement/task ID, git commit/tree hash, exact command, fixture hashes, expected behavior/oracle source, actual result, exit code, versions, artifact paths, timestamp, status and reviewer conclusion. PASS requires inspection of assertions and actual output. NOT_RUN, SKIPPED, INAPPLICABLE and external study-provenance limitations are separate. Maintainer independently reruns gates after each Packet Implementer receipt and records its review; Independent Reviewer never receives test logs.

Final acceptance includes clean setup, supported examples, all eligible engine paths, report content/visual inspection, scientific release profile, source/data separation and SSD checksum delivery. A code-only handoff may truthfully list remaining external study limitations, but missing core implementation/testing cannot be hidden behind them.
