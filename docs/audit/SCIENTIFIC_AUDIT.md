# Scientific audit of the recovered study

The principal disease/direct-treatment contrasts show no protein at BH-FDR < 0.05. The archive's broader claim across all tested contrasts needs correction: **Negr1 / A0A8I6AJV3** has log2FC -1.07835 and q=0.0169815 in the joint PreDM+dapa-versus-CTRL contrast. This is an auxiliary sensitivity result, not evidence of treatment rescue. The corresponding axis-specific treated-versus-control contrast has no q<0.10 proteins.

| Primary axis contrast | Proteins tested | Raw P<0.05 | Raw-P/effect candidates | q<0.10 | Minimum q |
|---|---:|---:|---:|---:|---:|
| PreDM vs CTRL | 3,635 | 34 | 26 | 0 | 0.99923 |
| PreDM+dapa vs PreDM | 3,635 | 430 | 316 | 0 | 0.11490 |
| DM vs CTRL | 3,601 | 105 | 80 | 0 | 0.99954 |
| DM+dapa vs DM | 3,601 | 195 | 146 | 0 | 0.47857 |

The stage-by-treatment interaction in the joint fit has minimum q=0.98498. A stronger-looking PreDM result than DM result does not establish a difference in treatment effects between stages. Separate axis fits change both the tested universe and variance pooling; their choice as primary is described in the archive but no prespecified analysis plan was supplied.

Within-collection/contrast pathway FDR produces 18 significant axis rows across 4,302 tests and 10 joint rows across 4,308 tests. As an additional audit diagnostic, pooling all four contrasts and three collections within each model family gives 2 and 1 significant rows respectively. This is a sensitivity to the definition of the testing family, not a mandatory replacement analysis or evidence that the original per-family BH arithmetic is wrong. All 24 stored within-family BH adjustments agree after matching R's NA handling. One joint DM-treatment Reactome pathway has a nonfinite P value, which should be explicit in reporting.

The strongest conclusion is that this is a useful exploratory secondary analysis with unvalidated rescue inference and incomplete experimental/computational provenance. It cannot establish upstream peptide/protein identification quality or tissue-specific mechanisms. The detailed methodological assessment and primary-source citations follow.

## Detailed methods review

Audit date: 2026-09-12. Scope: static inspection of the archived abundance-refit R pipeline, supported by primary papers and official software documentation. Historical scripts were not executed or modified. Archive comments and generated notes were treated as evidence of intended interpretation, not instructions. Findings below do not establish that any particular biological result is true or false.

All script references use lines in the preserved directory:
`pipeline/` in the recovered repository (R source is byte-identical to the archive).

## Overall appraisal

The core limma analysis is a defensible starting point for secondary analysis of processed protein abundances. The strongest shortcomings concern the inferential interpretation of reversal scores, incomplete pathway provenance/diagnostics, and the inability of this downstream code to validate upstream mass-spectrometry processing. Exploratory labels already present in the scripts are valuable and should be retained prominently. Publication of the historical code is reasonable with the audit attached; presentation as a scientifically validated, reusable end-to-end proteomics pipeline would require further work.

## Findings requiring correction or explicit limitation

### S1. High: disease-axis score tests reuse selection data

**Observed:** `02_reversal.R:93-105` selects proteins and their signs from disease-versus-control results, constructs scores from those same animals, and tests treated versus untreated scores at lines 110-123. `05_axis_specific_limma.R:127-165` repeats this procedure. The permutation functions only shuffle already constructed scores; they do not refit disease effects or repeat selection and weighting. The fallback top-50 selection in `02_reversal.R:95-97` is also data-dependent.

**Consequence:** these score-test P values are not calibrated confirmatory evidence for a treatment effect after this selection. Selection and testing share the untreated group, so independence under the null is not established. The general requirement for independent selection and subsequent inference is demonstrated in [Kriegeskorte et al., 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2841687/). Application to these scripts is the audit's inference from the code, not a claim studied directly in that paper.

**Repair:** retain scores as descriptive, or learn the protein set/signs on independent data and validate the fixed score prospectively. A permutation-based alternative must respect the experimental randomization and repeat the entire selection and scoring procedure for every admissible treatment assignment; this is valid only under the corresponding exchangeability/randomization assumptions. Merely fixing the denominator in S3 does not repair this problem. With few biological replicates, splitting the current data may be too unstable to be useful.

The archived script already calls these scores exploratory and mentions regression to the mean (`02_reversal.R:148-149`). Those statements are scientifically appropriate, but exporting ordinary P-value columns can still invite overinterpretation.

### S2. High interpretation risk: shared untreated group induces apparent opposition

**Observed:** disease effect is `U-C`, and treatment effect is `T-U` (`02_reversal.R:30-32`; `05_axis_specific_limma.R:27-31`). Reversal is then derived from their ratio and signs.

**Audit derivation:** assuming independent group-mean sampling errors,

`Cov(U-C, T-U) = -Var(U)`.

Thus noisy positive disease estimates tend to accompany negative treatment estimates even without a true treatment effect. Selecting large disease estimates intensifies this regression-to-the-mean concern. This is an algebraic property of the actual contrasts, not evidence that a specific reported reversal is false. Disease-vs-control and direct treatment contrasts can still be estimated in the same model; the error is interpreting their opposition as independent corroboration. Direct treatment contrasts with their uncertainty remain the appropriate first assessment.

The two stage-specific models also share CTRL (`05_axis_specific_limma.R:28,30`), so calling them “independent” at line 213 means separately fitted, not statistically independent evidence.

### S3. Medium: exhaustive permutation calculation adds the observed case twice

**Observed:** `02_reversal.R:35-41` and `05_axis_specific_limma.R:58-63` enumerate every distinct allocation with `combn`, including the observed allocation, then return `(extreme + 1)/(all + 1)`.

**Consequence:** for this exhaustive absolute-tail test the exact probability is `extreme/all`. The extra correction makes values conservative but not the exact exhaustive probabilities advertised by the column name. With four animals per group, there are `choose(8,4)=70` allocations; if the observed split and its complement are the only extremes, the code gives `3/71` instead of `2/70`. This example follows directly from enumeration; it is not a claim about the archived data.

The +1 correction is appropriate to randomized sampling arrangements that need to include the observed statistic, as discussed in [Phipson and Smyth, 2010](https://arxiv.org/abs/1603.05766). Official [SciPy permutation-test documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.permutation_test.html) explicitly distinguishes exhaustive proportions from randomized +1 adjustments and notes that two-sided conventions differ. Preserve this pipeline's stated absolute-tail convention when repairing it, and use a numerical tolerance for floating-point ties.

The broader audit's independent `data_checks.json` alongside this audit now confirms this in the archived stage-specific scores: PreDM is stored as 0.0422535 (=3/71), whereas exhaustive fixed-score enumeration gives 0.0285714 (=2/70); DM is stored as 0.0985915 (=7/71), whereas fixed-score enumeration gives 0.0857143 (=6/70). These recalculations remain conditional on selected scores and must not be presented as repaired, selection-valid treatment P values.

### S4. Medium: “Full reversal” includes unlimited overshoot

**Observed:** `05_axis_specific_limma.R:113-120` labels every reversal index of at least 0.80 as “Full reversal”, with no upper bound. The formal-rescue rule at lines 122-123 requires opposing disease/treatment effects but does not test restoration to control.

**Audit derivation:** let `d=U-C`, `t=T-U`, and `RI=-t/d`. The treated residual is `T-C=d+t=d*(1-RI)`. Therefore RI=3 is called full reversal even though the treated mean is twice as far from control as the untreated mean. RI>2 exceeds the original disease deviation in the opposite direction.

**Repair:** label directional reversal and overshoot separately; report the treated-versus-control effect and interval, already fitted at lines 82-88. If “restoration” means practical equivalence to control, prespecify an equivalence margin and suitable inference. A nonsignificant treated-versus-control test alone does not demonstrate equivalence. This recommendation follows from the contrast algebra and the different null hypothesis needed for a restoration claim.

The broader audit's data checks find four of ten PreDM and four of fourteen DM “Full reversal” candidates with RI>1.2, but none with RI>2. Thus the unbounded classification is an observed code defect; a treated deviation exceeding the original deviation is a hypothetical failure mode, not an observed result in these candidates. Neither stage has a formal statistical rescue under the archived rule.

### S5. Medium: pathway findings require a precise null and correlation sensitivity

**Observed:** both enrichment scripts give `fgseaMultilevel` a fixed vector of signed `-log10(P)` ranks (`03_enrichment.R:45-68`; `06_axis_specific_enrichment.R:31-50`). No expression matrix or animal labels enter the full enrichment test.

**Interpretation:** this is preranked enrichment against random gene sets, not biological-sample randomization and not an independent replication of the protein analysis. [Official GSEA documentation](https://docs.gsea-msigdb.org/GSEA/GSEA_User_Guide/) explains the gene-set permutation null for preranked analysis and the importance of rank weighting and ties. [Wu and Smyth's CAMERA paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC3458527/) demonstrates that competitive tests based on gene permutation can be anti-conservative under inter-gene correlation and provides a correlation-aware alternative. This motivates a sensitivity check; it does not prove these archived pathways are false positives.

**Repair:** describe the competitive rank-based null, and assess important findings with a correlation-aware method such as CAMERA on a suitably prepared gene-level abundance matrix and the prespecified contrasts. Compare signed moderated-t ranks as a documented sensitivity analysis; signed P-value weighting is a choice, not automatically an implementation error. Standard weighted GSEA uses metric magnitudes, so these choices can change results even when rank order is similar.

### S6. Medium: gene-set provenance is discarded; rat sets are ortholog mapped

**Observed:** `03_enrichment.R:28-29` and `06_axis_specific_enrichment.R:16-18` request rat gene symbols from msigdbr and retain only gene symbol and set name. They do not explicitly set database species or export the database version, original identifiers, orthology evidence, full gene-set snapshot, or a snapshot hash.

Official [msigdbr reference documentation](https://igordot.github.io/msigdbr/reference/msigdbr.html) shows a human database default and rat ortholog mapping; the [msigdbr introduction](https://igordot.github.io/msigdbr/articles/msigdbr-intro.html) documents `db_version`, source identifiers, `ortholog_sources`, and `num_ortholog_sources`, and cautions that predicted gene orthology does not establish conservation of complete pathways. Current documentation is not evidence of the exact database release used in the old run.

**Repair:** describe results as human MSigDB collections mapped to rat if confirmed against the historical installed version; record the exact package/database version and export a checksum-verified snapshot plus mapping coverage. Report unmapped measured proteins, duplicate gene mappings, and pathways retained after overlap filtering. Preserve external-resource attribution and redistribution terms.

### S7. Medium: warning suppression and omitted Monte Carlo diagnostics weaken auditability

**Observed:** `03_enrichment.R:68-71,132-136` and `06_axis_specific_enrichment.R:49-54` suppress fgsea warnings and omit `log2err` from exported results. Stage 03 also catches a missing collection and proceeds with other collections (`25-43`). Its summary does not distinguish a deliberately excluded collection from a retrieval failure.

**Consequence:** downstream readers cannot assess the discarded numerical diagnostics from these tables or reliably distinguish complete from partial enrichment runs. `eps=0` itself is a valid fgsea option: the [official fgsea tutorial](https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html) documents it and shows `log2err` in returned output. It permits smaller estimated P values, not biological certainty.

**Repair:** preserve warnings, errors, `log2err`, finite-result counts, and expected/loaded collection manifests; fail a required collection or explicitly mark the run incomplete. Record the RNG seed, package versions, and parallel backend.

### S8. Moderate caveats: duplicate collapse, testing families, and selective stability

**Observed:** `03_enrichment.R:49-50` and `06_axis_specific_enrichment.R:35-36` retain the most extreme protein rank for each gene. This is deterministic for a fixed order but gives genes with multiple entries more opportunities to supply an extreme score. The broader audit's data checks report 877 repeated nonblank gene-symbol rows in the input, so the collapse affects a substantial part of the mapping. Use a prespecified aggregation/representative rule and report sensitivity of important pathways; the exact effect on pathways has not been recalculated here.

Each fgsea call adjusts within one collection and contrast (`03:63-71`; `06:45-55`). ORA does BH before removal of zero-overlap rows (`03:75-84`), which appropriately includes those tested sets in that adjustment. The exported adjusted P values do not establish a single FDR guarantee over all collections, stages, model variants, or analysis choices. Report the actual testing family, or calculate a broader adjustment if making a study-wide family claim.

**BH verification detail:** independent reimplementations must match R's NA handling. The [R source for p.adjust](https://svn.r-project.org/R/trunk/src/library/stats/R/p.adjust.R) removes NA entries from local `p` before first evaluating default `n=length(p)`. Because default arguments are evaluated lazily in the function's local frame ([R Language Definition, argument evaluation](https://stat.ethz.ch/CRAN/doc/manuals/R-lang.html#Argument-evaluation)), the default count is the nonmissing length; NA positions are restored in the result. Supplying an explicit original-vector length would be different. Current [fgseaMultilevel source](https://raw.githubusercontent.com/ctlab/fgsea/master/R/fgseaMultilevel.R) calls default `p.adjust(pval, method="BH")` after combining result rows (lines 216-219). Thus a check that counts failed/NA pathways in its multiplier can spuriously flag an adjustment discrepancy. Compare using nonmissing P values before alleging a code defect; exact historical package provenance remains desirable.

The treatment-pathway leave-one-out analysis starts with full-run significant pathways (`03:96-100,114-119`) and reruns selected sets, reporting nominal P-value stability (`143-146`). The code explicitly identifies this as sensitivity analysis. That is a valid descriptive use, not held-out validation or independent confirmation. It does not repair the pathway correlation assumptions or selection issue.

## Defensible choices and remaining scientific boundary

**Limma model.** `01_qc_dea.R:120-140` uses a group-means model and explicit contrasts; `05_axis_specific_limma.R:80-85` separately fits each stage. `trend=TRUE` allows an intensity-dependent prior variance; `robust=TRUE` robustifies empirical-Bayes variance estimation, not individual-sample regression or confounder adjustment. These are valid limma options, subject to design/residual diagnostics and the supplied experimental structure. See the [official limma manual, eBayes section](https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf). Choosing a joint or stage-specific fit changes variance pooling and the tested universe; neither is inherently wrong, but a primary strategy should be fixed before choosing the most favorable result.

**Effect-size cutoff.** The code combines zero-null BH-adjusted P values with observed absolute log2FC cutoffs (`01:157-159`; `05:41-43`). That is a transparent descriptive selection rule, but it does not test that the true effect exceeds 0.25. For that inferential claim use a threshold-null method such as TREAT, with a biologically justified threshold and a new BH adjustment. [McCarthy and Smyth, 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2654802/) establishes this distinction. Do not silently change the historical method, and do not describe raw-P candidates as FDR-confirmed proteins.

**Missingness and input scale.** `01:34-47` reads an already processed log2 protein matrix, retains proteins observed at least three times in any one group, and keeps missing values in primary limma fitting (`136-140`). PCA-only median replacement is separate (`59-76`). Stage-specific fitting instead requires at least three observations in every relevant group (`05:71-74`), improving coverage while changing scope. No blanket imputation in the primary analysis is defensible, but available-case estimation does not resolve abundance-dependent censoring. [Lazar et al., 2016](https://pubmed.ncbi.nlm.nih.gov/26906401/) shows why missing-value methods depend on mechanism and processing stage. [Karpievitch et al., 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2723007/) explicitly models both random and abundance-dependent missingness. Assess missingness by group and abundance, contrast estimability, and sensitivity of important proteins; do not choose imputation solely to increase significant hits.

The broader audit quantifies 3,714 proteins by 20 samples (four per group), 1.81% missing overall, and one sample with 11.25% missing. Retention is 3,714 in the joint model, 3,635 for PreDM, and 3,601 for DM. These rates bound the observed scope of missingness; they neither establish its mechanism nor justify automatic sample exclusion.

The sensitivity routine called “MinDet” draws random values from a down-shifted normal distribution (`01:209-221`), so the name does not accurately describe a deterministic minimum replacement. Document its actual formula (sample median minus 1.8 sample SD, width 0.3 SD); assess its assumptions rather than relying on the label.

**Upstream validation boundary.** These R stages do not establish sample provenance, tissue, randomization, biological versus technical replication, batch effects, search database and version, peptide/protein identification FDR, protein grouping, contaminant/decoy handling, intensity normalization, or any prior imputation. Whether other archive files supply these is for the broader project audit. The code explicitly says tissue is unavailable (`03:181`) and no batch was supplied (`01:120,296`); this is an uncertainty to document, not permission to invent covariates. [HUPO-PSI MIAPE](https://www.psidev.info/miape) specifies the need for sample and analytical metadata sufficient to interpret and reassess proteomics experiments. A processed abundance refit must not be called a validated raw-MS-to-result pipeline without the missing upstream evidence.

## Priorities for a maintained successor

1. Keep the historical snapshot immutable; make scientific changes in a separately versioned maintained implementation.
2. Disable confirmatory interpretation of selected-axis score P values; correct exhaustive arithmetic and separate overshoot labels in the successor.
3. Declare the primary estimand, contrasts, universe/filter, testing families, and exploratory status before new analysis.
4. Restore provenance for processed input and animal metadata; pin dependencies and gene-set snapshots.
5. Retain all enrichment diagnostics and assess consequential pathway conclusions with correlation-aware sensitivity analysis.
6. Validate the successor on synthetic edge cases (selection under a null, exact permutation enumeration, overshoot, missing groups, duplicate mappings) and then reproduce the archived data outputs with documented expected changes.
