# Methods and capabilities of the maintained pipeline

This page lists what the maintained successor implements today, which declarations it accepts, and where each method is specified. The normative rules are the scientific methods contract ([SM01–SM41](../../specs/001-downstream-proteomics/contracts/scientific-methods.md), v1.3.0). Every typed refusal and warning the code can emit is listed in [reason-codes.md](reason-codes.md). Exact, tested commands are in the [usage guide](../user-guide/usage.md). Release and versioning rules are in the [release guide](../user-guide/release.md).

**Status.** Phases 1–3 (R01–R11), PERMANOVA (R13) and post-differential analysis (R14a–R14f) are implemented. Their acceptance cases pass on **synthetic fixtures**, self-verified by the operator-authorized Claude route; the independent review of 2026-10-05 covered R14a–R14f, the calibration and R12 (findings fixed); the later adaptive-policy redesign (ADR 0010) is self-verified. R12 private regression and release gates are Maintainer-only and **NOT_RUN**. This is **not a v1.0 release**: v1.0 also needs the Maintainer-only private regression (V111, V114) and the release decisions below. Publication and disclosure review are operator decisions (AGENTS.md); the license is stated below.

**License.** The software is available under the [PolyForm Noncommercial License 1.0.0](../../LICENSE.md) (decided by the Maintainer, Omar A. Solliman, on 2026-10-08; D-63). Academic, research, personal and non-profit use is permitted; any commercial use requires a separate written license from Omar A. Solliman.

## Capabilities

| Capability | Packet | Method (contract rules) |
|---|---|---|
| `foundation.io_roundtrip` | R01 | Python ↔ R typed request/result exchange, hashed artifacts, offline execution (SM25). |
| `intake` | R02 | Protein-level wide/long TSV, mapped protein tables and the recovered legacy workbook through a declared crosswalk; scale states and masks retained (SM01–SM03). |
| `preprocessing` | R03 | Declared normalization (`preserve`, `median`, `reference`, `tmt_loading`), technical-replicate aggregation, coverage masks and QC and declared exclusions; no silent second normalization (SM02, SM04–SM06). |
| `design` | R04 | Frozen design grammar: groups, centred continuous covariates, categorical covariates, interactions, subject blocking; exact rank, aliasing and contrast estimability, planned before any fit (SM07–SM08). |
| `limma` | R05 | Observed-data limma with exact sparse/weighted contrasts, `trend`/`robust` empirical Bayes, TREAT, prespecified BH families (SM09, SM11–SM12). |
| `assay_engines` | R06 | DEqMS (declared peptide/PSM counts) and proDA (native dropout likelihood) as declared adapters compared with the primary model (SM10, SM13). |
| `resources` | R07 | Hashed, versioned offline snapshots (mapping, gene sets, orthology, PPI) prepared by `proteomics resources prepare`; protein-to-gene mapping; no runtime downloads (SM14–SM15). |
| `pathways` | R08 | Design-valid enrichment: CAMERA, ROAST, fgsea on moderated statistics, ORA with the tested universe, enrichment summaries (SM16–SM18). |
| `response` | R09 | Descriptive treatment response, equivalence, and restricted independent score inference (SM19–SM24). |
| `report_stub`, `report_full` | R10a, R10b | Thin and full offline HTML reports built only from hash-verified stage outputs; every figure has a source table (SM25). |
| `permanova` | R13 | PERMANOVA/PERMDISP with restricted permutation schemes, exact enumeration when the admissible set is small (SM26–SM30). |
| `post_de_sets` | R14a | Set grammar over completed families, UpSet/Venn, concordance, overlap P only for unit-disjoint contrasts (SM31). |
| `post_de_sensitivity` | R14b | Re-planned covariate/subgroup/interaction models, matched-n draws, leave-one-unit-out influence (SM32). |
| `post_de_association` | R14c | Protein–phenotype association with the primary engine or permutation correlation; aliased pooled analyses refused (SM33). |
| `post_de_biomarker` | R14d | Leakage-safe nested cross-validation, whole-procedure permutation test, fixed panels, locked-model external validation (SM34–SM38). |
| `post_de_networks` | R14e | Opt-in co-abundance modules (signed weighted correlation or correlation-distance clustering), module–trait tests, PPI connectivity against a measured-universe null (SM39–SM40). |
| `post_de_eligibility` | R14f | Integrated eligibility/dependency report: state, typed reason, claim label and inputs of every module (SM41). |

R11 reproduction tooling (environment locks, offline reproduction, resume, calibration, benchmark and the validation ledger) is a library and command set, not a pipeline stage; see the [usage guide](../user-guide/usage.md). The `legacy` capability slot is not a pipeline stage either. Maintainer-local regression against the original archive is the separate script `scripts/maintained/legacy_regression.py` (R12, see the [Maintainer runbook](../validation/013-release/MAINTAINER_RUNBOOK.md)).

## Accepted declarations

- **Assays:** `lfq_dda`, `lfq_dia`, `tmt`, `processed_protein`. **Input formats:** `wide_tsv`, `long_tsv`, `mapped_protein`, `legacy_workbook`.
- **Source scales:** `linear_positive`, `log2`, `ratio_log2`, `standardized_unknown`. A second logarithm is refused (`E_SCALE_SECOND_LOG`).
- **Designs:** independent groups, paired and repeated observations with subject blocking (`none`, `fixed_subject`, `duplicate_correlation`), centred continuous exposures, categorical covariates and interactions.
- **Engines:** `limma` (primary), `deqms`, `proda`. Missing R packages give NOT_RUN, never INAPPLICABLE.
- **Coverage policies:** `available_case`, `native_dropout`. **Sensitivity imputations** (never for the primary model): `min_deterministic`, `left_shifted_gaussian`, `knn`, `complete_case`.
- **Post-DE (`post_de` block, Phase 4 only):** set rules (`union`, `intersect`, `difference`, `complement_within_tested`); association `model` or `correlation`; biomarker selection `top_k_auc`, `top_k_t`, `elastic_net`, `lasso` and classifiers `penalized_logistic` (glmnet), `svm_linear`, `svm_polynomial`, `svm_radial` (e1071), `random_forest` (optional); co-abundance rules `wgcna_signed`, `hclust_correlation` (base R plus dynamicTreeCut; the WGCNA package is not used, D-44). The block and its defaults are specified in [post-de.md](../../specs/015-post-de-analysis/contracts/post-de.md).

## Notes on specific post-DE options

- **Adaptive policy (D-59, ADR 0010).** Post-DE modules and PERMANOVA adapt to the data and never hold up the DE results. Only an invalid or non-computable piece is refused, and the valid alternative runs automatically. Every adaptation (requested, used, reason) is listed in the module's `eligibility.json`, a `W_POST_DE_ADAPTED`/`W_PERMANOVA_ADAPTED` warning, and the report's adaptation table. The refusal audit is D-60.
- **Association correlations.** `pearson` and `spearman` are plain correlations. With an adjustment (pooled analyses adjust for group by default, plus any `adjust_for`), their partial form runs: Pearson on residuals, or Spearman as Pearson on ranks residualised on the adjustment columns. The rows record `correlation_requested` and `correlation_used`. `adjusted_for` lists only the adjustment actually used.
- **Aliasing.** The pooled estimate is refused (`E_PHENOTYPE_ALIASED`) when the phenotype is determined by group, even if the group term is not fitted (D-56), and the within-group analyses run instead (D-59). The Simpson flag compares the within-group slopes with the unadjusted pooled slope (D-57).
- **Co-abundance.** It runs normally at `min_units`, as exploratory from 4 units, and is refused below 4. Repeated designs use subject means when group is constant within subject; when group varies within subject, subject effects are removed with the group term kept. The `row_transform` and `exploratory` columns record which applied.
- **Biomarker compute budget.** `compute_policy: adapt` (default) scales repeats, then B, down to `compute_guard_hours`, with minimums of 10 repeats and B = 99. The attainable P floor is recorded. `refuse` is the strict mode. Below 4 units per class, nested CV is not computable and is refused.
- **Single-feature direction `train_only`.** Despite its name (kept because it is part of the frozen contract), this rule picks each feature's direction on the same data that the AUC is computed on. Such AUCs are labelled `in_sample` and flagged direction-optimistic. Use `prespecified` for a direction fixed in advance (D-57).
- **PPI connectivity null.** Degree-preserving draws from the measured, mapped universe, one derived L'Ecuyer stream per declared set (`stream` column). The degree-binned null is approximate; in calibration it was mildly conservative on random sets.

## Claim labels for post-DE results

`descriptive`, `exploratory_raw_p`, `in_sample`, `cross_validated_nested`, `fixed_panel_cv`, `independently_validated`, `module_level` (frozen vocabulary). A module that did not complete has no label and shows no values. NOT_RUN (software or implementation unavailable) and INAPPLICABLE (scientifically ineligible, typed reason) are shown as different states.

## What is out of scope

Raw-MS search and quantification, PTM localization, single-cell data, web UI or databases, network deployment, clinical or diagnostic claims, and causal drug-mechanism claims. Network and module results make no mechanism claim; biomarker results are never called diagnostic or validated unless a locked model was evaluated once on a disjoint declared cohort.

## Calibration evidence

- Core limma calibration (R11, V105/V106): [012-validation/calibration](../validation/012-validation/calibration/).
- Post-DE calibration (nested-CV AUC and permutation P under the null, leaky reference, connectivity null, and a subject-blocked scenario with whole-subject permutation and an ungrouped reference): [015-post-de-analysis/R14f/calibration](../validation/015-post-de-analysis/R14f/calibration/). Only fixed ridge λ with top-k AUC selection and repeated k-fold CV are calibrated; LOOCV, lasso/elastic-net, SVM, random forest and tuned λ grids are not.
