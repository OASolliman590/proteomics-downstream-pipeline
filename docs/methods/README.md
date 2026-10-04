# Methods and capabilities of the maintained pipeline

This page lists what the maintained successor implements today, which declarations it accepts, and where each method is specified. The normative rules are the scientific methods contract ([SM01–SM41](../../specs/001-downstream-proteomics/contracts/scientific-methods.md), v1.3.0). Every typed refusal and warning the code can emit is listed in [reason-codes.md](reason-codes.md). Exact, tested commands are in the [usage guide](../user-guide/usage.md). Release and versioning rules are in the [release guide](../user-guide/release.md).

**Status.** Phases 1–3 (R01–R11), PERMANOVA (R13) and post-differential analysis (R14a–R14f) are implemented. Their acceptance cases pass on **synthetic fixtures**, self-verified by the operator-authorized Claude route; the independent audit of the post-DE work has not been run. R12 private regression and release gates are Maintainer-only and **NOT_RUN**. This is **not a v1.0 release**: v1.0 also needs the Maintainer-only private regression (V111, V114) and the release decisions below. No license has been selected; licensing and publication are operator decisions (AGENTS.md).

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

## Claim labels for post-DE results

`descriptive`, `exploratory_raw_p`, `in_sample`, `cross_validated_nested`, `fixed_panel_cv`, `independently_validated`, `module_level` (frozen vocabulary). A module that did not complete has no label and shows no values. NOT_RUN (software or implementation unavailable) and INAPPLICABLE (scientifically ineligible, typed reason) are shown as different states.

## What is out of scope

Raw-MS search and quantification, PTM localization, single-cell data, web UI or databases, network deployment, clinical or diagnostic claims, and causal drug-mechanism claims. Network and module results make no mechanism claim; biomarker results are never called diagnostic or validated unless a locked model was evaluated once on a disjoint declared cohort.

## Calibration evidence

- Core limma calibration (R11, V105/V106): [012-validation/calibration](../validation/012-validation/calibration/).
- Post-DE calibration (nested-CV AUC and permutation P under the null, leaky reference, connectivity null): [015-post-de-analysis/R14f/calibration](../validation/015-post-de-analysis/R14f/calibration/).
