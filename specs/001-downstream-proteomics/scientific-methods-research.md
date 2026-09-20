# Scientific method contracts for the maintained downstream proteomics pipeline

Research date: 2026-09-12. Scope: protein/protein-group abundance matrices and associated experimental metadata. Raw-MS search, peptide identification, protein inference, and PTM-site analysis are outside this capability boundary. No private data were downloaded and no historical pipeline was executed. Recommendations labeled **Contract** are proposed engineering/scientific decisions, not claims that a cited author mandates that exact implementation.

Sources below were opened using the web tool. Official package documentation/source takes precedence for API details. Some package source links point to a moving development branch: the implementation must pin a tested release and its exact APIs, record versions and arguments, and avoid adopting newly added options merely because they appear on `master`. There is no universally best downstream method across assays and preprocessing histories.

The release manuals retrieved on this date report limma 3.68.5 and DEqMS 1.30.0 on Bioconductor 3.23. These are candidate release targets, not evidence of the versions used by the historical project, and not a substitute for testing the selected dependency lock.

## 1. Capability policy and mandatory core

**Contract:** implement all named capability adapters and their tests, while making execution conditional on documented eligibility. An ineligible explicitly requested method is an actionable validation error, never an automatic switch to another backend. Unrequested methods may be `not_requested`; an eligible zero-discovery result is `completed` with zero discoveries. Missing dependencies, unsupported designs, nonconvergence, and empty eligible universes are distinguishable from scientific null findings. Never fabricate a successful output or substitute historical output tables for execution.

| Capability | Mandatory implementation | Execution/default contract |
|---|---|---|
| Input/provenance, QC, contrast construction, estimability, multiplicity, output/report schemas | Core | Always |
| limma differential abundance, zero-null and TREAT | Core | Default model for declared processed log2 protein abundances; threshold inference is selected explicitly |
| Batch covariates; paired fixed-subject model; one repeated-measures block with limma consensus correlation | Core | Only when identifiable from supplied metadata |
| LFQ median/reference normalization and unchanged-input mode | Core | Prespecified from input history; no significance-driven normalization selection |
| TMT protein-level channel/run/reference handling | Core | Only for declared TMT/iTRAQ protein abundance data with compatible channel and plex provenance |
| CAMERA and ROAST/mroast | Core adapters | CAMERA for eligible fixed-effects matrix analyses; ROAST for separately specified self-contained tests, including eligible blocked designs |
| fgsea and exploratory ORA | Core adapters | Explicit exploratory competitive analyses with their own universes/families |
| DEqMS | Assay-qualified adapter | Requires genuine compatible quantified peptide/PSM counts |
| proDA | Assay-qualified adapter | LFQ data with preserved un-imputed observation/missingness information; fixed-effects design |
| Reversal descriptors; independent fixed signatures; restricted permutation; equivalence | Core adapters | Descriptors default. Inferential extensions require explicit eligibility and prespecification |
| Synthetic calibration and backend parity | Core | Required release evidence for every enabled inferential capability |

**Mandatory provenance fields:** assay and acquisition type; feature-unit definition; abundance scale and units; normalization status/method; imputation status/method and original missingness-mask availability; biological-unit and technical-replicate IDs; condition/time/covariates; batch/plex/channel roles where applicable; known upstream identification filtering; organism/taxonomy; mapping and gene-set provenance. Unknown values must be represented as unknown, not inferred from file names or intensity ranges.

**Contract:** processed input with incomplete upstream history can still undergo clearly labeled conditional secondary analysis. It cannot be silently represented as originally un-imputed LFQ, or as validated identification/quantification output. Uncertainty in tissue, sample origin, randomization, and upstream processing must survive into the final report.

## 2. Limma: estimand, variance model, and exact contrast calculation

The [official eBayes/TREAT documentation](https://raw.githubusercontent.com/bioc/limma/master/man/ebayes.Rd) distinguishes intensity-dependent prior-variance moderation (`trend`) from robust estimation of variance hyperparameters (`robust`). Robust empirical-Bayes moderation is not robust regression against individual outlying animals, batch correction, or a substitute for biological replication. TREAT tests against a nonzero fold-change interval, whereas eBayes tests a zero effect; a posthoc absolute-logFC cutoff does not turn a zero-null P value into threshold-null inference. The ordinary and threshold tests must remain separately labeled.

**Contract:**

- Prespecify a formula, contrast vector, model scope, sample IDs, feature universe, and test family. Default processed-input model is limma with `trend=TRUE, robust=TRUE`; expose explicit configuration and diagnostics for other scientifically justified choices. Record mean-variance trend and residual/variance summaries.
- Use log2 effect units only when the input scale supports that interpretation. Log2 ratios are already log2 quantities; do not transform them again. Reject incompatible raw-count/statistic matrices masquerading as protein intensities.
- When an effect threshold is scientifically specified, use standard TREAT and `topTreat`; record the null interval and threshold rationale. Do not use an observed cutoff to claim that the true effect exceeds that cutoff. Do not automatically enable a new `upshot` variant or a changed empirical-Bayes algorithm; such choices need separate versioned validation.
- Export native effect, unscaled SE, moderated SE, residual df, total df, raw P, adjusted P, test kind, and estimability/coverage status. Ordinary confidence intervals must be identified as such; they are not simultaneous post-selection confidence intervals.
- Prespecified omnibus tests and stage-by-treatment interactions are distinct estimands. “Significant in one stage, nonsignificant in another” is not a test of stage specificity; use an interaction contrast.

The [official contrasts.fit documentation](https://raw.githubusercontent.com/bioc/limma/master/man/contrasts.fit.Rd) states an easily missed limitation: nonorthogonal designs combined with missing values or precision weights can yield approximate unscaled standard errors. One-way group-means models without blocking are an exact exception. **Contract:** for general weighted/missing-data/batch/paired models, fit a reparameterized design in which each requested contrast is a coefficient, or compute the exact feature-specific covariance. Validate against direct least-squares/GLS fits. Do not silently accept the approximation while claiming exact general contrast support.

## 3. Experimental design, replication, and estimability

The [limma manual](https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf) describes `duplicateCorrelation` as a consensus within-block correlation obtained from gene-wise REML, requiring sufficient observations beyond the design rank. The [current duplicateCorrelation source](https://raw.githubusercontent.com/bioc/limma/master/R/dups.R) explicitly warns and returns zero correlation if every block is a singleton, or if the block factor is already encoded in the fixed design. That behavior must not become an unnoticed successful repeated-measures analysis.

**Contract:**

1. Build and validate the actual model matrix before fitting. Report rank, aliases, residual degrees of freedom, sample exclusions, empty factor levels, contrast vectors, and condition/batch/subject cross-tabs. Require a full-rank parametrization and estimable declared contrasts; never drop an aliased covariate automatically.
2. Repeat estimability/coverage checks on each feature's observed samples. A globally full-rank design can become locally unidentifiable under missingness. Record nonestimable feature-contrast rows as such, rather than P=1 or effect=0. A documented minimum-residual-df and biological-replicate policy is required; four-per-group is a property of the historical example, not a universal requirement.
3. For paired within-subject contrasts, support a fixed subject effect plus the within-subject factor. This cannot identify a between-subject condition effect that is fully absorbed by subject indicators.
4. Separately support one exchangeable repeated-measures block using `duplicateCorrelation`, then pass the estimated `block` and `correlation` into the actual `lmFit`. Do not include the same subject factor twice. Save the consensus value and available gene-wise diagnostics; fail explicit repeated-block requests on singleton/no-estimable-correlation data.
5. Multiple crossed/nested random effects and random slopes are outside this adapter's supported design class. If future scope requires them, add and validate a real suitable backend; do not flatten them into a single block or silently fit independent samples.
6. Biological units determine replication. Technical replicates must either be aggregated by a declared rule before modeling, or explicitly modeled. Pooled TMT references and empty channels are never biological replicate observations.
7. Batch is a model covariate when identifiable. Complete confounding is not solved by `removeBatchEffect`, centering, ComBat, or adding an aliased formula term. Any batch-adjusted PCA is a separate visualization derivative; inferential fitting should retain the measurement matrix and explicit design.

Acceptance fixtures must cover balanced/unbalanced independent groups, matched pairs, repeated measures, batch-confounded failure, missing feature-specific cells, zero residual df, malformed weights, and changed sample ordering. Positive sample-correlation estimates must not be confused with within-protein residual block correlation.

## 4. Normalization and missingness contracts

The [proDA introduction](https://const-ae.github.io/proDA/articles/Introduction.html) demonstrates log2 LFQ data with NA-coded missing values and median normalization. It also explains why censoring can distort observed intensity distributions, so equalizing entire distributions can be inappropriate. This supports assay-aware diagnostics, not a universal command to normalize every matrix. The [proDA median-normalization reference](https://const-ae.github.io/proDA/reference/median_normalization.html) defines its operation relative to a pseudo-reference; it is not interchangeable with subtracting each sample's raw median.

**Contract for generic/LFQ input:**

- Require explicit linear/log2/other scale. Diagnostics can flag a contradictory declaration, but must not silently choose a scale by numerical range. Convert zeros to missing only when the exporter/assay declaration says zero is a missing-value sentinel. A zero log2 ratio is a valid value.
- Preserve the original matrix and original missingness mask. Reject infinities and invalid negative linear intensities. Record all transforms and matrices by stable IDs/hashes.
- Implement `none` for already processed input and prespecified median/reference-based normalization for suitable unnormalized LFQ. Define the exact formula, reference-feature set, and missing-data behavior; do not label different formulas with one generic “median normalized” name.
- Normalize against prespecified stable/spike-in standards when that is the scientific design. Median-centering requires an approximately stable reference proteome or balanced changes; it can erase genuine global abundance shifts. Quantile/VSN options, if included, must be opt-in sensitivity methods with their own assumptions and units.
- Primary available-case limma does not impute; it also does not model nonrandom missingness. QC must report missingness by sample, condition, abundance, and observed-count patterns, including all-missing groups. Separate observed estimates from missingness-driven evidence.
- Visualization imputation must have its own derivative matrix and never leak into inference. Any explicitly requested inferential imputation sensitivity must retain method/seed/mask and avoid presenting filled values as independent measured information. Unknown or previously imputed input is ineligible for a preserved-dropout proDA claim unless the exact original missingness information and required measurements are restored.

**Contract for protein-level TMT/iTRAQ input:** require unique plex/run/channel identity, sample/biological-unit mapping, channel role, reference identity/composition, and the declared upstream normalization and isotopic-correction state. Implement within-plex protein-level normalization and cross-plex reference scaling only when requested and applicable, with exact formulas and QC. A bridge/reference rule must check connected comparability across plexes; disconnected treatment/plex groups cannot be made identifiable by normalization. Exclude reference/empty channels from biological contrasts after normalization. Preserve plex/batch effects in the model where identifiable. Repeated technical runs do not increase biological n.

The [official MSstatsTMT workflow](https://bioconductor.org/packages/release/bioc/vignettes/MSstatsTMT/inst/doc/MSstatsTMT.html) starts from feature-level records with Protein, PSM, TechRepMixture, Mixture, Run, Channel, Condition, BioReplicate, and Intensity, and distinguishes global normalization, protein summarization, and reference normalization. **Contract:** do not claim to reproduce that complete method by fabricating PSM rows from protein summaries. A protein-only adapter may perform an explicitly named protein-level reference transformation, but cannot recover peptide-level normalization, interference correction, or ratio compression. It must not treat TMT missingness as ordinary independent LFQ dropout without validation.

**Implementable protein-level reference option (proposed formula):** after any explicitly chosen within-plex channel offset correction, let `b[g,p]` be the median log2 abundance in designated equivalent bridge channels for protein g in plex p; let `a[g]` be the median available bridge value across the prespecified connected plex set. For biological sample s in plex p, set `z[g,s]=y[g,s]-b[g,p]+a[g]`. Record each b and a. A missing bridge value makes that protein/plex transformation ineligible; do not replace it silently. A within-plex median option can use `m[s]=median_g(y[g,s])` over a declared reference-feature set and subtract `m[s]-median_s_in_plex(m[s])`. These formulas are explicit protein-level engineering choices whose assumptions and calibration must be tested; they are not a claim to reproduce all MSstatsTMT steps. Absence of a bridge does not automatically invalidate a balanced, identifiable within-plex treatment comparison with plex effects: the no-bridge route must be explicitly modeled and cannot claim reference-based cross-plex normalization.

## 5. DEqMS adapter: genuine counts and correct statistics

The [DEqMS 1.30.0 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DEqMS/inst/doc/DEqMS-package-vignette.html) uses `lmFit -> contrasts.fit -> eBayes`, adds a protein-aligned `fit$count`, and then calls `spectraCounteBayes`. Its TMT example uses minimum quantified PSM count across experiments. Its LFQ example uses minimum quantified unique-plus-razor peptide counts, adding one to every count when zeros require a documented pseudocount. It requires at least two observed values per comparison group in that example. Output contains both limma and DEqMS statistics.

The [official DEqMS manual](https://bioconductor.org/packages/release/bioc/manuals/DEqMS/man/DEqMS.pdf) specifies an eBayes fit plus count covariate as input, and native `sca.t`, `sca.p`, `sca.dfprior`, `sca.priorvar`, and `sca.postvar`. `outputResult` exposes count-adjusted fields separately from the original limma columns. These facts support the following adapter contract; they do not authorize arbitrary count surrogates.

**Contract:**

- Accept quantified PSM or peptide counts with a declared source/definition, measured feature unit, sample/plex scope, aggregation rule, and optional explicit pseudocount. Missing count annotations are not zero. Do not use protein intensity, amino-acid length, annotation count, number of replicates, or arbitrary ones as substitutes.
- Join counts by unique feature ID after filtering and verify an exact one-to-one alignment. Model covariates must be finite and strictly positive after an explicit supported adjustment; preserve original counts. Counts with no useful variability or failed variance-count fits must produce a diagnostic failure, never silent limma fallback.
- The canonical adapter's statistic/P/q must come from `sca.t`, `sca.P.Value`/`sca.p`, and recomputed count-adjusted BH; original `P.Value` and `adj.P.Val` must have an unambiguous `limma_` prefix in native exports. Preserve variance-count fit diagnostics.
- Do not call the combined result TREAT or assume that DEqMS supplies TREAT's threshold-null test. Nonzero threshold inference is outside this backend unless separately derived and validated. Do not combine its count prior with other priors merely by chaining method names.
- Initial supported model class should be fixed-effect designs whose exact contrasts and coverage pass validation. Any blocked extension requires explicit parity/calibration evidence before being declared supported.

Acceptance tests: shuffled count-table rows preserve results; duplicate/missing/zero/negative count failures; a real synthetic count-variance relationship produces native DEqMS output; every exported P/q matches count-adjusted native fields, not the adjacent limma fields.

## 6. proDA adapter: dropout likelihood has an applicability boundary

The [official proDA model reference](https://const-ae.github.io/proDA/reference/proDA.html) accepts a protein-by-sample matrix with NA missingness and a fixed linear design. It estimates dropout curves and location/variance priors; location moderation handles proteins unobserved in one condition. The returned fit includes convergence information. This is not a generic mixed-model interface.

The [official test_diff reference](https://const-ae.github.io/proDA/reference/test_diff.html) accepts an expression or string linear contrast tested against zero, or a reduced model for an F-test. Wald output includes `diff`, `se`, `df`, `t_statistic`, `pval`, `adj_pval`, `n_obs`, and `n_approx`; the latter is inferred information, not a biological replicate count. There is no documented TREAT threshold argument in this interface.

**Contract:**

- Eligible input is LFQ abundance with authentic observation/missingness information and a declared normalization history compatible with the model. Assay qualification must be explicit for DDA or DIA rather than assuming acquisition types have identical censoring.
- Fit after the chosen normalization. Preserve NA entries; never impute first. Already imputed input without a recoverable original mask fails this adapter. Do not attempt to infer a mask from repeated low values.
- Do not reuse limma's requirement for observations in every group: all-missing-in-one-condition proteins are part of proDA's intended use. Apply an explicitly documented minimal total-observation/information policy; exclude globally unobserved features from substantive claims and flag heavily prior-driven results. Report per-condition observations alongside inferred uncertainty.
- Support declared fixed batch/subject covariates only when the design is identifiable. Reject requested random block/correlation, random slopes, and unsupported assay combinations. Do not discard those terms and proceed.
- Export convergence, iteration limits, dropout parameters, uncertainty fields, and warnings. Distinguish a no-missing-data case where dropout parameters may be NA from a nonconverged fit. Nonconverged parameters cannot yield a completed inferential run.
- Preserve this backend's native scale/statistic/null. A zero-null result plus effect filtering is not a TREAT result; keep backend-specific result fields in addition to the common schema.

Acceptance tests should use the package's synthetic generator or independently generated LFQ-like data, with complete data, intensity-dependent dropout, an all-missing condition, and broken input/design cases. Verify that exact missing-mask restoration differs from treating imputed values as measurements.

## 7. Gene mapping, finite universes, and pathway backends

**Contract for mapping/universes:** retain protein groups as measured units in protein analysis. Never split one ambiguous group into multiple independently tested genes. Build a separate pathway matrix/rank table with one row per unambiguous target gene. Default representative selection must be prespecified and independent of contrast significance, for example highest observed coverage, then median abundance, then stable feature ID; record excluded/merged members. Expression-level aggregation is a separately named, validated option. Never select the largest absolute t/P rank as an unnoticed default.

For each contrast/backend save the input feature list, mapping table/version, unambiguous gene list, eligibility exclusions, pathway membership after intersection, eligible pathway sizes, and hashes. A finite DE-statistic universe can differ from a complete-matrix universe; that difference must be explicit. Do not claim sensitivity analyses use the same universe when they do not. Require at least one finite outside-set gene for competitive tests, and nonempty within-set membership. Prespecify size cutoffs; 10-500 is a reasonable configurable engineering default, not a universal scientific constant.

The [msigdbr API reference](https://igordot.github.io/msigdbr/reference/msigdbr.html) separates source database species from output species. The [msigdbr introduction](https://igordot.github.io/msigdbr/articles/msigdbr-intro.html) exposes database version and orthology provenance and cautions that predicted gene orthology does not establish complete pathway conservation. **Contract:** explicitly set both species; snapshot full memberships and source/output IDs, orthology sources/support, release/version, terms of use, and checksums. Deduplicate target-gene membership within each set, report one-to-many/many-to-one mapping, and require an explicit mapping policy. Rat output from human MSigDB must be described as ortholog-mapped rather than rat-native curated pathways. Network/resource failure cannot silently remove a required collection.

### CAMERA

The [official CAMERA documentation](https://raw.githubusercontent.com/bioc/limma/master/man/camera.Rd) defines a competitive test. The [implementation](https://raw.githubusercontent.com/bioc/limma/master/R/geneset-camera.R) has default `inter.gene.cor=0.01`; this is a fixed assumption. `NA`/`NULL` requests estimation. Extra arguments are warned as disregarded, so passing `block`/`correlation` does not create repeated-measures support.

**Contract:** use an explicitly finite complete-row gene matrix and fixed-effect design. The correlation-aware primary configuration is `inter.gene.cor=NA`, with `allow.neg.cor=FALSE`, explicit contrast and `trend.var`; retain estimated correlations. A fixed-correlation sensitivity, e.g. 0.01, must say it is fixed. `robust=TRUE` is not a supported CAMERA argument; do not pretend its internal moderation is identical to robust limma. Fail blocked design requests. Export `NGenes`, `Direction`, `PValue`, estimated `Correlation` when available, native FDR, and centrally defined q values. Test gene-set indices, universe sizes, complete finite matrix, rank, and residual df before invocation.

### ROAST/mroast

The [official ROAST documentation](https://raw.githubusercontent.com/bioc/limma/master/man/roast.Rd) explicitly forbids NA/Inf input, identifies a self-contained null, and recommends at least 9,999 rotations for publication. It distinguishes directional and mixed alternatives and notes native mroast's mid-P behavior. The [mroast source](https://raw.githubusercontent.com/bioc/limma/master/R/geneset-roast.R) passes `block` and `correlation` into model-effect calculation.

**Contract:** accept the documented complete finite pathway matrix; preserve blocked structure and the correct correlation estimate. Use `mroast` with explicit `midp=FALSE` so primary q values can be reproduced from ordinary P values under the declared central adjustment. Directional and mixed tests are distinct registered hypotheses/families; do not choose the smaller after fitting. Preserve native `NGenes`, `PropDown`, `PropUp`, `Direction`, `PValue`, `FDR`, `PValue.Mixed`, and `FDR.Mixed`. Fix RNG seed/streams, rotations, statistic, and approximation setting. Do not call rotation simulation an exhaustive sample permutation.

When a matrix pathway backend requires complete rows, selecting that complete-row universe must be an explicit pipeline operation with attrition/coverage reporting. If the eligible universe is inadequate, fail the requested pathway analysis. Never fill missing values just to make the backend run or silently switch the null to fgsea.

### fgsea and ORA

The [official fgsea tutorial](https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html) documents preranked enrichment, Monte Carlo calculation, `eps=0`, `log2err`, and finite measured backgrounds for ORA. **Contract:** default rank to a finite signed native zero-null test statistic, not a threshold-filtered hit list; preserve backend identity. Check ties and make ordering deterministic without inventing random biological evidence. Retain all fgsea diagnostics/warnings and leading edges, and state that the competitive rank null does not model animal-label randomization. ORA uses a prespecified eligible measured universe and a declared candidate rule; it remains exploratory when foreground selection is nominal or posthoc. Never use the whole genome by default or convert pathway significance into protein-level significance.

## 8. Multiplicity and canonical result schema

**Contract:** a test registry must identify the estimand, backend, null kind, contrast, gene/protein universe, pathway collection(s), directionality, and family membership before fitting. Export raw P, native package q, and pipeline q with separate names. Primary family adjustment must include all eligible tests specified for that family, not only plotted/nonzero-overlap/significant tests. Technical failure remains NA with a failure reason; it is never a P=1 observation and never silently changes a declared family. A family with unexplained failed tests is incomplete.

BH is a declared conventional method, not proof of study-wide FDR over every exploratory comparison or correlated method choice. Support explicit within-contrast and combined prespecified families; provide Holm/BY when chosen for their respective control goals. Do not select the correction that produces the most hits. Native mroast mid-P and R's missing-P default semantics must not cause a false adapter-parity failure: compare the actual configured hypothesis set and algorithm.

Required common columns include feature/set ID, contrast ID, backend, test kind, estimate and units where defined, SE/df where defined, statistic name/value, raw P, family ID, q method/value, observation counts, and eligibility/status/reason. Not every backend has every field: undefined SE or effect for a pathway is explicit null, not a zero-valued fabricated estimate.

## 9. Reversal, signatures, exact randomization, and restoration

**Mandatory descriptive reversal:** define `d=untreated-control`, `t=treated-untreated`, and `r=treated-control=d+t`. Report all three and their uncertainty. `RI=-t/d` is undefined/unstable near a zero disease effect; use a prespecified denominator guard and do not hide the excluded count. Report directional change, residual magnitude, and overshoot separately. RI=1 means the point estimates meet control; RI>1 is overshoot; RI>2 moves farther from control than the untreated point estimate. These are algebraic descriptions, not statistical rescue claims. Disease and treatment effects share untreated observations and have negatively correlated sampling errors.

**Independent fixed signature:** protein IDs, directions, weights, centering/scaling reference, mapping policy, and missing-feature policy must be frozen from independent training data or prior biological knowledge. Test only on distinct biological units. Do not use disease-selected features/signs from the validation animals and then pretend that ordinary score tests validate a fixed signature. Plot-oriented signatures from current data remain descriptive.

The [Winkler et al. primary paper](https://orbi.uliege.be/bitstream/2268/210172/1/Winkler%20et%20al.%20-%202014%20-%20Permutation%20inference%20for%20the%20general%20linear%20model.pdf) distinguishes exchangeability, nuisance effects, and within/whole-block transformations. The [official SciPy permutation documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.permutation_test.html) distinguishes complete enumeration from randomized sampling and documents tail conventions and tie precision.

**Restricted permutation contract:**

- Require a declared randomization/exchangeability scheme, unit, block restrictions, statistic, tail, finite-data rule, and enumeration cap. Independent two-group assignments preserve group sizes; paired swaps preserve pairs; clustered assignment permutes whole valid clusters rather than individual observations. An arbitrary batch column is not sufficient to establish exchangeability.
- For a finite admissible set, include the observed assignment exactly once, enumerate unique assignments, and use `extreme/N`. For a separately specified Monte Carlo scheme, include the observed statistic via its appropriate `(extreme+1)/(B+1)` construction. Never mix denominators. Save N/B, extreme count, ties/tolerance, seed, and resolution.
- In a nested treatment-randomization test of a selected disease axis, keep control samples fixed and exchange treated/untreated labels only within admissible disease-stage/randomization units. Refit disease effects, repeat feature selection and direction/weight construction, and recompute the score statistic for every assignment. Cache only provably label-invariant operations. A prespecified zero-selected-feature branch must return a defined neutral statistic for that assignment, not silently skip the assignment and shrink the null space.
- If labels were not randomized and identical-distribution exchangeability is not justified, do not call the result an exact causal treatment test. Heteroscedastic equal means alone do not justify raw mean-difference permutation. Freedman-Lane residual permutation with nuisance covariates is not universally exact; it requires a separate explicitly approximate supported contract, not automatic fallback.

### Equivalence of treated and control

The [TOSTER authors' introduction](https://aaroncaldwell.us/TOSTERpkg/articles/IntroductionToTOSTER.html) distinguishes nonsignificance from equivalence and defines two one-sided tests against prespecified bounds. The [current tsum_TOST reference](https://aaroncaldwell.us/TOSTERpkg/reference/tsum_TOST.html) provides a reference implementation for simple independent/paired test fixtures.

**Contract:** equivalence bounds are scientific inputs in log2 units, supplied before inspection of results. For residual effect estimate `r`, SE `s`, df `v`, and bounds `L<U`, compute `p_lower=P(t_v >= (r-L)/s)`, `p_upper=P(t_v <= (r-U)/s)`, and `p_equivalence=max(p_lower,p_upper)`. At nominal alpha, the associated `100*(1-2*alpha)%` interval must lie wholly inside the bounds. Report ordinary 95% intervals separately. A BH-adjusted equivalence decision does not correspond to blindly applying an unadjusted 90% interval rule.

Apply equivalence inference across its prespecified eligible family, not only proteins that happened to pass disease screening. For limma-moderated SE/df this is a model-based extension whose calibration must be tested; do not imply it is the identical procedure to unmoderated TOSTER. Support simple validated limma-based residual contrasts first, and reject backend combinations without an established SE/df contract. Disease change, opposing treatment change, and equivalence remain separate evidence fields unless a separately specified conjunction test is implemented and calibrated. A small equivalence P value alone is not evidence that disease was originally present or that treatment caused restoration.

## 10. Calibration and release validation

[Morris, White and Crowther](https://discovery.ucl.ac.uk/10066118/1/2019%20-%20Morris%20-%20simulation%20studies%20tutorial%20-%20stat%20med.pdf) recommend explicit simulation aims, data-generating mechanisms, estimands, methods, and performance measures (ADEMP), with Monte Carlo uncertainty and reproducible random-number handling. The following concrete validation grid is a proposed pipeline contract, not a guarantee established by that paper.

### Deterministic tests and native-backend parity

1. Hand-computable independent, paired, batch, and interaction contrasts; matrix/data permutations; exact feature-specific SE; TREAT zero threshold compared with zero-null behavior; genuine positive thresholds with effects just below/above the boundary.
2. Counts shuffled/duplicated/missing for DEqMS; canonical sca rather than limma fields; proDA missing-mask/nonconvergence cases; all adapter eligibility errors.
3. Mapping ambiguity, gene duplicates, gene-set overlap and finite outside-set complement; no pathways/genes; full/partial collection failure; no hidden backend substitution.
4. Exhaustive 4+4 test with 70 allocations; exact pairs with `2^n` swaps; tied statistics; restricted blocks; impossible assignments; zero-selected-set nested branches. Exact results must agree with independently written enumeration using the same tail.
5. Equivalence at and beyond both bounds, equivalence without zero-null significance, zero-null significance without equivalence, overshoot including RI>2, and zero denominators.
6. Direct pinned-package parity for all scientifically meaningful exported native quantities. A successful mocked function call is not backend validation; every adapter needs a real executable synthetic example.

### Simulation grid

Generate truth independently of the estimator. Include a well-specified base scenario for each model and separate stress scenarios. Vary biological n, effect size, abundance-dependent variance, peptide/PSM count variance, sample-level heterogeneity, correlated pathways, batch balance/confounding, repeated-unit correlation, random missingness, LFQ-style dropout, all-missing groups, and TMT plex/reference connectivity. Stress cases must be labeled stress cases; a method is not expected to be universally calibrated under incompatible assumptions.

Measure empirical type-I error, per-dataset FDP and its mean (FDR), discovery count, power, effect bias, interval coverage, equivalence false-positive rate at both boundary nulls, nonconvergence/failure frequency, and sensitivity to preprocessing. Under the global null, FDP is 1 when any discovery occurs and 0 otherwise; pooling all protein calls across simulations is not the same as empirical FDR. For mixed null/non-null data compute FDP within each simulated dataset first. Correlated genes mean feature-level outcomes are not independent binomial replicates: assess Monte Carlo uncertainty across independently simulated datasets.

**Proposed acceptance design:** freeze scenarios and acceptance bounds before implementing the release. Use a quick deterministic/small-simulation CI suite for gross regressions and a separate reproducible calibration suite with enough independently simulated datasets for useful MC precision (e.g. 2,000 for inexpensive base models; expensive backends sized by prespecified MCSE). A nominal 0.05 test should be assessed with a declared tolerance and an uncertainty interval, not a brittle equality check. One practical release criterion is an upper 95% MC bound no greater than 0.075 for base-scenario type-I error/FDR, with a lower 95% bound of at least 0.925 for nominal 95% coverage; these tolerances are engineering choices to record explicitly in the validation plan and must not be widened after seeing failures. Power is descriptive and must never be optimized at the expense of error control.

Require the old invalid fixed-selected-score procedure as a labeled negative-control demonstration, while the independently fixed-score and correctly nested-randomization implementations pass their own assumptions. Preserve scenario seeds/streams, DGP parameters, truth, method configurations, failures, summary intervals, and plots. Report every prespecified scenario, including failures. Calibration passing establishes behavior on those scenarios, not universal validity or validation of the historical biological conclusion.

## 11. High-risk implementation traps to include in the one-shot specification

- No universal “best method” selector based on the largest discovery count.
- No inferred scale, unknown-as-false provenance, sample-label parsing masquerading as metadata, or silent missing-value conversion.
- No batch removal that conceals nonestimability; no technical channels counted as biological n.
- No TREAT claim from an absolute effect filter; no DEqMS result using adjacent ordinary-limma P columns; no proDA call on silently imputed observations.
- No unsupported CAMERA block arguments; no claim that its fixed 0.01 correlation was estimated; no assumption that CAMERA repeats robust limma identically.
- No ROAST/Camera imputation solely to satisfy a complete-data requirement; no unnoticed universe changes or native mroast mid-P adjustment.
- No protein-group duplication into gene-level pseudo-replicates; no maximum-significance representative default; no unversioned orthology resources.
- No selected-score P values advertised as independent validation; no exhaustive +1 correction; no skipping inadmissible/failed permutation statistics after enumeration; no nonsignificance-as-restoration.
- No adapter declared complete without native execution, negative applicability tests, numerical parity, and its prespecified calibration evidence.
