# Feature Specification: Post-differential-abundance analysis (result structure, robustness, phenotype association, biomarker evaluation, co-abundance and networks)

**Phase:** 4 (v1.1-post-de). **Packets:** R14a–R14f (one slice, six dispatch units, following the R10a/R10b precedent). **Status:** operator-authorized scope amendment of 2026-10-02 ([ADR 0009](../../docs/adr/0009-post-de-scope-amendment.md)); not yet implemented. It adds FR-131–FR-167, T131–T167 and V131–V167 and does not renumber any existing identity.

## Scope

Deliver FR-131–FR-167 for two new journeys:
- **US8:** the analyst structures and stress-tests differential results. This covers declared set logic, concordance, covariate/subgroup/resampling/influence robustness, and protein–phenotype association.
- **US9:** the analyst evaluates candidate biomarkers and co-abundance structure with leakage-safe procedures whose claims match their evidence.

Normative statistics are SM31–SM41 in the [scientific contract](../001-downstream-proteomics/contracts/scientific-methods.md#post-differential-analysis). Configuration, outputs and vocabulary are in [contracts/post-de.md](contracts/post-de.md).

Every module consumes frozen plan artifacts and completed R05/R08/R13 outputs read-only and verifies their hashes. No module changes the primary matrix, the primary families or any earlier result.

### Generality principles (normative for every packet in this slice)

1. **Declarative, study-agnostic.** No study, species, group, covariate or panel name appears in code. Groups, contrasts, covariates, subgroups, set rules, phenotypes and panels are declared in configuration and validated against the frozen plan. There is no free-form code evaluation.
2. **Declared eligibility.** Every module publishes an eligibility rule over the frozen design: design type, biological-unit n per class, balance, covariate aliasing and the identifier/species resources required. The planner records ELIGIBLE or INAPPLICABLE with a typed reason before any computation. Missing software is NOT_RUN, never INAPPLICABLE. An ineligible request never silently degrades to a weaker method.
3. **Design coverage.** Each module states its behaviour for:
   - independent two-group and k-group designs;
   - multi-factor designs;
   - paired and repeated designs with subject blocks;
   - continuous exposures;
   - severely unbalanced designs.

   Where a design is unsupported, the module refuses with a typed error.
4. **Biological units.** All resampling, cross-validation, permutation and influence operations act on biological units (SM03), never on technical injections or repeated observations as if they were independent.
5. **Small-sample-safe defaults.** Exact or permutation procedures where asymptotics are unreliable. Small-n warnings are part of the outputs, not just the logs.
6. **Claims vocabulary.** Each output carries a frozen claim label, such as `descriptive`, `in_sample`, `cross_validated_nested` or `independently_validated`. The report renders claims only at the level the evidence supports.

## Requirements

### R14a — Result structure and concordance

- **FR-131 — Declarative set algebra:** The system MUST build feature sets from completed differential families using a finite grammar over named contrasts and models (`union`, `intersect`, `difference`, `complement_within_tested`). It MUST use declared membership criteria (family q threshold, or TREAT; raw-P criteria allowed only with the label `exploratory_raw_p`) and optional direction, and export the full membership matrix.
- **FR-132 — Region provenance and refusals:** Every exported set or region MUST carry its rule, criterion, threshold, family IDs and null type. The system MUST refuse rules that reference undefined contrasts, mix zero-null with TREAT memberships, or combine families from different primary matrices.
- **FR-133 — Exact region partition:** The system MUST compute the exclusive regions of up to the declared number of sets (UpSet for any k; circle-based Venn for k ≤ 3). Region counts MUST sum to the size of the union, and protein-group versus unique-gene counts MUST be reported separately.
- **FR-134 — Direction-aware overlap:** For any two sets the system MUST report concordant-up, concordant-down and discordant counts, plus the features significant in only one.
- **FR-135 — Effect concordance across contrasts and models:** The system MUST report descriptive concordance of log2 effects on common tested features: Pearson and Spearman correlation, sign agreement among features significant in either, and an origin-constrained slope. The label is descriptive.
- **FR-136 — Overlap inference eligibility:** A hypergeometric overlap P is reported only when the two contrasts share no biological units. Otherwise overlap is descriptive and the output states why (shared control/arm dependence).

### R14b — Robustness and sensitivity

- **FR-137 — Covariate-imbalance diagnostics:** The system MUST cross-tabulate declared covariates against groups. It MUST detect aliasing (e.g. a covariate level present in only one group) and classify each requested adjustment as estimable-additive, estimable-interaction or non-estimable, before any refit.
- **FR-138 — Covariate-adjusted sensitivity model:** The system MUST re-plan the frozen design with declared additive covariates as a named secondary sensitivity model through R04 rank/estimability checks, fit it with the primary engine settings, and keep its own family.
- **FR-139 — Subgroup re-analysis:** The system MUST re-plan the design on a declared observation subset (e.g. one covariate level) as a named sensitivity model. It MUST refuse with a typed reason when any required group falls below the coverage rules, and never redefine the primary families.
- **FR-140 — Primary-versus-sensitivity comparison:** For every sensitivity model the system MUST report retained, lost and gained members per declared criterion, effect correlation, origin-constrained attenuation slope and Jaccard index. It MUST warn that differing significance is not evidence of an interaction. Where estimable, it offers the formal group × covariate interaction contrast as its own family.
- **FR-141 — Matched-n resampling sensitivity:** The system MUST repeatedly subsample biological units of larger groups to declared target sizes. With a recorded seed and draw count, it refits and reports the distribution of discovery counts and per-feature selection frequency (descriptive).
- **FR-142 — Unit-level influence on discoveries:** Building on SM13, the system MUST omit each biological unit (with all its observations) in turn, refit the fixed model, and report per-unit changes in discovery membership and the maximum absolute effect change.
- **FR-143 — Robustness summary:** The system MUST produce a per-feature robustness table across all requested sensitivities: fraction of analyses in which the feature meets the criterion, sign stability and influence flags. It is labelled descriptive and never re-adjusts primary P values.

### R14c — Protein–phenotype association

- **FR-144 — Continuous and ordinal phenotype models:** The system MUST test association of each feature with a declared numeric or ordinal phenotype using the primary engine. Proteins are the response; the phenotype is a design term, adjusted for declared covariates and, by default, for group. Results go in a separate family per phenotype.
- **FR-145 — Correlation alternative:** On request, the system MUST report Pearson or Spearman (and partial, given covariates) correlations with permutation P values on biological units, and BH within the phenotype family.
- **FR-146 — Missing phenotype values:** Phenotype-missing units are excluded per phenotype (complete case) with n reported. Phenotypes are never imputed.
- **FR-147 — Confounding guard:** The system MUST refuse a pooled association when the phenotype is aliased with group (e.g. measured in one group only, or perfectly separated), unless within-group analysis is declared. It MUST warn when pooled and within-group directions disagree (Simpson pattern).
- **FR-148 — Association outputs:** The system MUST export feature × phenotype effect, CI, P and q tables and a clustered heatmap with source tables, and report them offline.

### R14d — Biomarker discrimination evaluation

- **FR-149 — Scope and claims vocabulary:** Biomarker evaluation covers binary declared contrasts only. Every result carries one claim label: `in_sample`, `cross_validated_nested`, `fixed_panel_cv` or `independently_validated`. The system never emits diagnostic, clinical-utility or deployment claims.
- **FR-150 — Single-feature discrimination:** The system MUST report per-feature AUC with DeLong and stratified-bootstrap CIs on biological units. Direction is either prespecified or chosen on training data only. Auto-direction in-sample AUC is labelled `in_sample` and flagged as direction-optimistic.
- **FR-151 — Minimum-size refusal:** The system MUST refuse classifier evaluation when any class has fewer than the declared minimum of biological units (default 5) or when cross-validation folds cannot contain both classes. The error is typed.
- **FR-152 — Leakage-safe preprocessing:** Every learned transform MUST be fitted on the training portion of each fold only and applied to the held-out portion. This includes scaling, centring, feature filtering, feature selection, imputation where declared, and hyperparameters.
- **FR-153 — Nested cross-validation:** The system MUST estimate performance by repeated stratified outer CV (LOOCV permitted for tiny n), with an inner CV for selection and tuning. Folds are grouped by subject for repeated designs and seeds are recorded. It reports pooled out-of-fold AUC, per-repeat AUC mean ± SD and the out-of-fold predictions.
- **FR-154 — In-fold feature selection and stability:** Supported selectors are top-k by training-fold AUC or moderated t, and penalised regression (elastic net/lasso). They run inside each training fold, with k or penalty tuned in the inner loop. The system MUST report per-feature selection frequency across all outer folds.
- **FR-155 — Classifier families and score orientation:** Supported classifiers are penalised logistic regression (default), SVM (linear, polynomial, radial; decision values) and random forest (optional dependency). Score orientation MUST be derived from the fitted model's class coding, never from outcome labels at evaluation time.
- **FR-156 — Permutation test of the whole procedure:** The system MUST permute class labels across biological units, respecting declared blocks, and rerun the entire nested procedure, including selection. It reports P = (k+1)/(B+1) with B recorded and a declared compute guard.
- **FR-157 — Fixed-panel evaluation and circularity:** A declared panel MUST carry provenance (`independent`, `same_data` or `unknown`).
  - It is evaluated by CV with only the classifier trained in-fold.
  - `same_data` and `unknown` panels are labelled `fixed_panel_cv` (selection not cross-validated) and reported beside a nested-reselection estimate of the same size.
  - Eligibility follows SM23.
- **FR-158 — Calibration and threshold metrics:** For probabilistic classifiers the system MUST report out-of-fold calibration (curve, Brier score) and sensitivity/specificity at a threshold chosen inside training folds by a declared rule, with bootstrap CIs. Thresholds tuned on test folds are forbidden.
- **FR-159 — External validation:** When an independent validation cohort is declared and passes the SM23 disjointness audit, the system MUST lock the discovery-trained model (features, transforms, weights, threshold) and evaluate it once on the validation cohort. Without one, the report states "not externally validated".
- **FR-160 — Biomarker outputs:** The system MUST export ROC curves with CI bands, CV performance distributions, selection-stability tables, permutation nulls, calibration plots and out-of-fold prediction tables with source data, and report them with the claim label.

### R14e — Co-abundance modules and protein interaction networks

- **FR-161 — Co-abundance module eligibility:** Co-abundance analysis runs only when declared, with at least the declared minimum biological units in the analysed set (default 20) and complete genuinely observed features. Otherwise it is refused with a typed reason. It is never on by default.
- **FR-162 — Module construction and stability:** The system MUST build modules from the declared correlation and network rule (e.g. signed weighted correlation with recorded soft-threshold selection, or hierarchical clustering of correlation distance). It reports module membership, eigengenes and bootstrap (unit-resampled) module stability.
- **FR-163 — Module–trait association:** The system MUST test module eigengenes against groups and phenotypes with the primary engine on eigengenes, in a separate family, labelled as module-level inference.
- **FR-164 — Offline interaction resources:** Protein–protein interaction analysis MUST read only a hashed, versioned local snapshot under SM14 rules (source, release, species, score type), prepared by a separate explicit action. There are no runtime downloads.
- **FR-165 — Network connectivity with a measured-universe null:** The system MUST report the induced interaction subnetwork of a declared set and test its connectivity against degree-preserving random sets drawn from the measured, mapped universe, not the whole genome. Hub and degree results are descriptive.

### R14f — Integration, eligibility reporting and generality qualification

- **FR-166 — Post-DE eligibility and dependency report:** The report MUST show, for every requested post-DE module, its eligibility decision and reason, claim label and inputs. It also shows a "what each conclusion depends on" summary linking discoveries to their robustness and sensitivity outcomes.
- **FR-167 — Generality qualification matrix:** Each R14 module MUST either complete or refuse with its declared typed reason across a synthetic design matrix. The matrix covers: two-group independent; three-group; paired; repeated with subject blocks; continuous exposure; 3-vs-3 small n; 15-vs-80 unbalanced; 200-unit large n; and human, mouse and rat identifier sets. No module may crash, silently degrade or emit an ineligible claim.

## Acceptance Scenarios

<a id="V131"></a>

### V131: Declarative set algebra

**Fixture:** Three completed synthetic contrasts with known q values and directions.

**Oracle:** Python set operations written in the test over features with q < 0.05.

**Exact assertion:** `union`, `intersect`, `difference` and `complement_within_tested` memberships equal the oracle. A raw-P criterion is labelled `exploratory_raw_p`. The membership matrix has one row per tested feature.

**Negative case:** An unknown operator or arbitrary expression string fails E_SETRULE_GRAMMAR.

**Contract:** SM31; **owner:** R14a.

<a id="V132"></a>

### V132: Region provenance and refusals

**Fixture:** The V131 contrasts plus a TREAT family and a family from a sensitivity matrix.

**Oracle:** Declared rule metadata.

**Exact assertion:** Every set carries its rule, criterion, threshold, family IDs and null type.

**Negative case:**
- Mixing zero-null with TREAT fails E_SETRULE_MIXED_NULL.
- An undefined contrast fails E_SETRULE_UNKNOWN_CONTRAST.
- Mixing primary-matrix and sensitivity-matrix families fails E_SETRULE_MIXED_MATRIX.

**Contract:** SM31; **owner:** R14a.

<a id="V133"></a>

### V133: Exact region partition

**Fixture:** Three and five overlapping synthetic sets; protein groups that share genes.

**Oracle:** Brute-force enumeration of the 2^k − 1 exclusive regions.

**Exact assertion:**
- Region counts equal the oracle and sum to the size of the union.
- Venn output is produced only for k ≤ 3.
- Protein-group and unique-gene counts are reported separately and differ where genes are shared.

**Negative case:** Requesting a Venn for k = 4 fails E_VENN_K; UpSet is still produced.

**Contract:** SM31; **owner:** R14a.

<a id="V134"></a>

### V134: Direction-aware overlap

**Fixture:** Two contrasts with planted concordant and discordant features.

**Oracle:** Hand-labelled feature list.

**Exact assertion:** Concordant-up, concordant-down, discordant and one-only counts and memberships equal the oracle.

**Negative case:** Discordant features counted as overlap fail.

**Contract:** SM31; **owner:** R14a.

<a id="V135"></a>

### V135: Effect concordance

**Fixture:** Two contrasts whose effects are a known linear transform plus noise.

**Oracle:** `cor()` (Pearson, Spearman) and `lm(y ~ 0 + x)` written in the test.

**Exact assertion:** The statistics match within 1e-10 on common tested features, sign agreement matches, and the output is labelled `descriptive`.

**Negative case:** No P value column is emitted for concordance.

**Contract:** SM31; **owner:** R14a.

<a id="V136"></a>

### V136: Overlap inference eligibility

**Fixture:** Contrasts A-vs-C and B-vs-C (shared C) and contrasts on disjoint cohorts.

**Oracle:** `phyper` upper tail for the disjoint case.

**Exact assertion:** The disjoint case reports the hypergeometric P (the universe is the features tested in both). The shared-C case has no P column and the reason `shared_biological_units`.

**Negative case:** A hypergeometric P for the shared-control case fails.

**Contract:** SM31; **owner:** R14a.

<a id="V137"></a>

### V137: Covariate-imbalance diagnostics

**Fixture:** Groups × sex with one group entirely male (as in confounded clinical designs) and a balanced variant.

**Oracle:** Cross-tabulation and design rank written in the test.

**Exact assertion:**
- Confounded design: the additive adjustment is classified estimable when rank allows, and the interaction non-estimable with the empty cell named.
- Balanced design: both estimable.

**Negative case:** Fitting the interaction in the confounded design fails E_SENSITIVITY_NONESTIMABLE before any fit.

**Contract:** SM32; **owner:** R14b.

<a id="V138"></a>

### V138: Covariate-adjusted sensitivity model

**Fixture:** A synthetic effect that is partly explained by a covariate.

**Oracle:** A direct limma fit with `~ group + covariate` and the frozen trend/robust settings, written in the test.

**Exact assertion:** Estimates, moderated statistics and family q equal the oracle. The model is named `sensitivity` with its own family. The primary family is byte-identical before and after.

**Negative case:** Writing sensitivity rows into the primary family fails.

**Contract:** SM32; **owner:** R14b.

<a id="V139"></a>

### V139: Subgroup re-analysis

**Fixture:** Three groups × sex, with a subgroup filter `sex == M`; a variant where one group has one male.

**Oracle:** A direct limma fit on the subset.

**Exact assertion:** The subset model matches the oracle and coverage rules are re-applied on the subset.

**Negative case:** The one-male variant is refused with E_SUBGROUP_COVERAGE and the group named; no partial fit is reported.

**Contract:** SM32; **owner:** R14b.

<a id="V140"></a>

### V140: Primary-versus-sensitivity comparison

**Fixture:** The V138 and V139 outputs.

**Oracle:** Set operations, `cor`, `lm(y ~ 0 + x)` and Jaccard written in the test.

**Exact assertion:**
- Retained, lost and gained sets, correlation, attenuation slope and Jaccard equal the oracle.
- The non-interaction warning is present.
- Where estimable, the interaction contrast appears in its own family.

**Negative case:** Labelling "significant in subgroup only" as an interaction effect fails.

**Contract:** SM32; **owner:** R14b.

<a id="V141"></a>

### V141: Matched-n resampling

**Fixture:** A 31-vs-11 synthetic design with known effects.

**Oracle:** The same seed and draw count reproduced by a direct loop in the test.

**Exact assertion:**
- Draws subsample biological units only.
- The discovery-count distribution and per-feature selection frequencies equal the oracle.
- The seed and draws are recorded.
- The output is labelled `descriptive`.

**Negative case:** Subsampling technical injections independently fails.

**Contract:** SM32; **owner:** R14b.

<a id="V142"></a>

### V142: Unit-level influence

**Fixture:** A design with one planted high-leverage biological unit carrying two technical injections.

**Oracle:** Refits omitting each unit, with all its injections, written in the test.

**Exact assertion:** Per-unit membership changes and the maximum |Δ effect| equal the oracle, and the planted unit ranks first.

**Negative case:** Omitting a single injection as if it were a unit fails.

**Contract:** SM32 (builds on SM13); **owner:** R14b.

<a id="V143"></a>

### V143: Robustness summary

**Fixture:** The outputs of V138–V142.

**Oracle:** Per-feature fractions computed in the test.

**Exact assertion:** Robustness fractions, sign stability and influence flags equal the oracle. The output is labelled `descriptive` and primary q values are unchanged.

**Negative case:** Any adjusted-P column derived from robustness fractions fails.

**Contract:** SM32; **owner:** R14b.

<a id="V144"></a>

### V144: Continuous phenotype association

**Fixture:** A phenotype linearly related to 10 planted features, with a group effect and a covariate.

**Oracle:** A direct limma fit with `~ group + covariate + phenotype` written in the test.

**Exact assertion:** Phenotype coefficients, moderated t, P and family q equal the oracle, and there is one family per phenotype.

**Negative case:** Pooling phenotype rows into the group-contrast family fails.

**Contract:** SM33; **owner:** R14c.

<a id="V145"></a>

### V145: Correlation alternative with permutation

**Fixture:** n = 8 units, with Spearman and partial-correlation requests.

**Oracle:** `cor.test` and an exact/Monte-Carlo permutation P over units with a recorded seed, written in the test.

**Exact assertion:** Correlations equal the oracle within 1e-12, permutation P equals (k+1)/(B+1), and BH is applied within the phenotype family.

**Negative case:** Permuting technical observations independently fails.

**Contract:** SM33; **owner:** R14c.

<a id="V146"></a>

### V146: Missing phenotype values

**Fixture:** A phenotype missing for 3 of 20 units.

**Oracle:** A complete-case fit on 17 units.

**Exact assertion:** Results match the oracle, n = 17 is reported, and no imputed phenotype values exist.

**Negative case:** Any imputation of the phenotype fails E_PHENOTYPE_IMPUTATION.

**Contract:** SM33; **owner:** R14c.

<a id="V147"></a>

### V147: Confounding guard

**Fixture:** (a) A phenotype measured only in one group. (b) A Simpson fixture where the within-group slopes are negative and the pooled slope is positive.

**Oracle:** Design rank and the signs of within-group and pooled slopes.

**Exact assertion:**
- (a) The pooled request is refused; the within-group request runs.
- (b) Both are reported, with the Simpson warning.

**Negative case:** The pooled association in (a) fails E_PHENOTYPE_ALIASED.

**Contract:** SM33; **owner:** R14c.

<a id="V148"></a>

### V148: Association outputs

**Fixture:** The V144 output.

**Oracle:** Source-table re-derivation of every heatmap cell.

**Exact assertion:** Tables, heatmap and source data agree cell-for-cell, and the report shows family, n and adjustment.

**Negative case:** A heatmap without a source table fails.

**Contract:** SM33, SM25; **owner:** R14c.

<a id="V149"></a>

### V149: Scope and claims vocabulary

**Fixture:** Every R14d output on synthetic data.

**Oracle:** The frozen claim vocabulary.

**Exact assertion:** Every result row and report section carries exactly one allowed claim label, and no output contains the words "diagnostic", "clinical utility" or "validated" unless the label is `independently_validated`.

**Negative case:** A forbidden claim string or a missing label fails E_BIOMARKER_CLAIM.

**Contract:** SM34; **owner:** R14d.

<a id="V150"></a>

### V150: Single-feature discrimination

**Fixture:** Two classes of 15 units with planted up- and down-shifted features.

**Oracle:** `pROC::roc` / `ci.auc` (DeLong and stratified bootstrap with a recorded seed), written in the test.

**Exact assertion:**
- AUC and CIs equal the oracle.
- Prespecified-direction AUC is below 0.5 for an opposite-direction feature.
- Auto-direction AUC is labelled `in_sample` with the direction-optimism flag.

**Negative case:** Unlabelled auto-direction AUC fails.

**Contract:** SM34; **owner:** R14d.

<a id="V151"></a>

### V151: Minimum-size refusal

**Fixture:** Classes of 4 vs 30 units; 5-fold CV with 4 minority units.

**Oracle:** Declared minimum = 5; fold-composition check.

**Exact assertion:** Evaluation is refused, and both reasons are reported.

**Negative case:** Any AUC emitted for the refused fixture fails E_BIOMARKER_SMALL_N / E_BIOMARKER_FOLDS.

**Contract:** SM35; **owner:** R14d.

<a id="V152"></a>

### V152: Leakage-safe preprocessing

**Fixture:** Pure-noise features (no class signal) with n = 20 per class and p = 2,000.

**Oracle:** Two runs of the same procedure:
1. A deliberately leaky reference implemented in the test (selection and scaling on all data, then CV).
2. The pipeline.

**Exact assertion:**
- The leaky reference shows inflated mean AUC (> 0.75 on this fixture and seed).
- The pipeline's nested estimate is centred near 0.5 (95% interval across repeats includes 0.5).
- Audit logs show every transform fitted per training fold.

**Negative case:** Any transform fitted on held-out samples fails E_BIOMARKER_LEAKAGE.

**Contract:** SM35; **owner:** R14d.

<a id="V153"></a>

### V153: Nested cross-validation

**Fixture:** A planted 5-feature signal; a repeated-measures variant with 2 observations per subject.

**Oracle:**
- An independent nested loop written in the test with the same seeds, folds and selector.
- For the repeated variant, fold-membership checks.

**Exact assertion:**
- Out-of-fold predictions, pooled AUC and per-repeat mean ± SD equal the oracle.
- In the repeated variant, no subject has observations in both training and test of any fold.

**Negative case:** Subject leakage across folds fails E_BIOMARKER_GROUP_LEAKAGE.

**Contract:** SM35; **owner:** R14d.

<a id="V154"></a>

### V154: In-fold selection and stability

**Fixture:** The V153 fixture.

**Oracle:** Per-fold selections recomputed in the test.

**Exact assertion:** Selection frequencies equal the oracle, planted features have frequency ≥ 0.9, and noise features are mostly < 0.2 on this seed.

**Negative case:** A selection computed on the full data fails.

**Contract:** SM35; **owner:** R14d.

<a id="V155"></a>

### V155: Classifier families and score orientation

**Fixture:** A linearly separable fixture, with class factor levels supplied in reversed order.

**Oracle:**
- Direct glmnet and e1071 fits in the test.
- An AUC computed from decision values with orientation derived from the model coding.

**Exact assertion:** AUC is identical regardless of factor-level order. The SVM uses decision values, not predicted classes. Random forest is NOT_RUN when not installed.

**Negative case:** An orientation flip that yields AUC = 1 − true fails.

**Contract:** SM35; **owner:** R14d.

<a id="V156"></a>

### V156: Permutation test of the whole procedure

**Fixture:** The V152 noise fixture and the V153 signal fixture; small B for tests.

**Oracle:** Label permutations over units that rerun the full nested procedure in the test.

**Exact assertion:**
- Noise gives a non-small P.
- Signal gives P at the floor 1/(B+1).
- P = (k+1)/(B+1).
- Block-restricted permutations are used for blocked designs.

**Negative case:** Permuting only the final classifier (with selection kept fixed) fails E_BIOMARKER_PERMUTATION_SCOPE.

**Contract:** SM36; **owner:** R14d.

<a id="V157"></a>

### V157: Fixed-panel evaluation and circularity

**Fixture:** A panel selected on the full noise fixture (`same_data`) and a panel declared `independent` on a separate synthetic cohort.

**Oracle:** In-fold-classifier-only CV and nested reselection, both in the test.

**Exact assertion:**
- The same-data panel is labelled `fixed_panel_cv`, shows optimistic AUC, and is reported beside the nested estimate (≈ 0.5).
- The independent panel is labelled `fixed_panel_cv` without the optimism flag.

**Negative case:** A same-data panel reported without the nested comparison or the flag fails.

**Contract:** SM37 (with SM23); **owner:** R14d.

<a id="V158"></a>

### V158: Calibration and threshold metrics

**Fixture:** The V153 fixture with logistic classifier probabilities.

**Oracle:** Brier score and binned calibration computed in the test from the out-of-fold probabilities; threshold rule applied on training folds only.

**Exact assertion:** Brier score, calibration bins and sensitivity/specificity (with bootstrap CI) equal the oracle.

**Negative case:** A threshold chosen on test folds fails E_BIOMARKER_THRESHOLD_LEAKAGE.

**Contract:** SM35; **owner:** R14d.

<a id="V159"></a>

### V159: External validation

**Fixture:**
- A discovery cohort plus a disjoint validation cohort with the same planted signal.
- A variant sharing 3 subjects between the two.

**Oracle:** A locked model applied once in the test; the SM23 disjointness audit.

**Exact assertion:** The disjoint case reports a validation AUC labelled `independently_validated`, with the locked model hash.

**Negative case:** The overlapping variant fails E_SCORE_SELECTION_OVERLAP. Re-tuning on the validation cohort fails E_VALIDATION_RETUNED.

**Contract:** SM38 (with SM23); **owner:** R14d.

<a id="V160"></a>

### V160: Biomarker outputs

**Fixture:** The V150–V158 outputs.

**Oracle:** Source-table re-derivation of every plotted ROC/CV/calibration point.

**Exact assertion:** Every figure has a source table, and the report shows the claim labels, n per class, CV scheme, seeds and permutation B.

**Negative case:** A figure without source data, or a missing claim label, fails.

**Contract:** SM34, SM25; **owner:** R14d.

<a id="V161"></a>

### V161: Co-abundance eligibility

**Fixture:** n = 12 and n = 60 synthetic cohorts.

**Oracle:** The declared minimum of 20.

**Exact assertion:** n = 12 is refused with E_COABUNDANCE_SMALL_N; n = 60 runs.

**Negative case:** Running by default without a declaration fails.

**Contract:** SM39; **owner:** R14e.

<a id="V162"></a>

### V162: Module construction and stability

**Fixture:** n = 60, with three planted correlated modules plus noise.

**Oracle:** Planted membership; a direct correlation/clustering computation written in the test.

**Exact assertion:**
- Recovered modules match the planted ones (adjusted Rand index ≥ 0.9 on this seed).
- Eigengenes equal the first principal component of each module's scaled data, up to sign.
- Bootstrap stability is reported with its seed.

**Negative case:** Feature-level resampling instead of unit resampling fails.

**Contract:** SM39; **owner:** R14e.

<a id="V163"></a>

### V163: Module–trait association

**Fixture:** One planted module associated with group.

**Oracle:** A direct limma fit on eigengenes.

**Exact assertion:** Statistics equal the oracle, the result is a separate family, and the label is `module_level`.

**Negative case:** Pooling module rows into protein families fails.

**Contract:** SM39; **owner:** R14e.

<a id="V164"></a>

### V164: Offline interaction resources

**Fixture:** A synthetic PPI snapshot with a manifest, and a variant with a hash mismatch.

**Oracle:** The manifest hash.

**Exact assertion:** A valid snapshot loads with its source, release, species and score type recorded.

**Negative case:** The hash mismatch fails E_RESOURCE_HASH. Any network call during analysis fails E_RESOURCE_DOWNLOAD.

**Contract:** SM40 (with SM14); **owner:** R14e.

<a id="V165"></a>

### V165: Network connectivity with a measured-universe null

**Fixture:** A synthetic network where the declared set is a planted dense cluster, and a random set.

**Oracle:** Degree-preserving sampling from the measured universe, with a recorded seed, written in the test.

**Exact assertion:**
- The planted set has P at the floor and the random set a non-small P.
- The null draws come only from measured mapped features.
- Hub and degree results are labelled descriptive.

**Negative case:** A null drawn from the whole snapshot (genome) fails E_NETWORK_UNIVERSE.

**Contract:** SM40; **owner:** R14e.

<a id="V166"></a>

### V166: Post-DE eligibility and dependency report

**Fixture:** A run where some modules complete and some are refused (small n, aliasing, missing optional package).

**Oracle:** The stage states and reasons from the run.

**Exact assertion:**
- The report lists every requested module with its state, typed reason and claim label.
- The dependency summary links each primary discovery to its robustness and sensitivity outcomes.
- NOT_RUN and INAPPLICABLE are rendered distinctly.

**Negative case:** A refused module shown as zero results, or NOT_RUN shown as INAPPLICABLE, fails.

**Contract:** SM41, SM25; **owner:** R14f.

<a id="V167"></a>

### V167: Generality qualification matrix

**Fixture:** The synthetic design matrix listed in FR-167, each cell with known planted structure.

**Oracle:** Per-cell expected state (COMPLETED or the specific typed refusal) and, for completed cells, the module's own V-case oracle on that cell.

**Exact assertion:** Every module × design cell reaches its expected state. Completed cells pass their oracle. No cell crashes, emits an ineligible claim or silently switches method.

**Negative case:** Any untyped exception, unexpected completion or unexpected refusal fails.

**Contract:** SM41; **owner:** R14f.
