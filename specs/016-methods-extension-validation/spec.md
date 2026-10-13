# Feature Specification: Methods extension and multi-study validation (Welch engine, offline pathway libraries, calibration extension, multi-study framework)

**Phase:** 5 (v1.2-methods-validation). **Packets:** R15a–R15d (four serial dispatch units). **Status:** Maintainer scope decision of 2026-10-08 ([ADR 0011](../../docs/adr/0011-methods-extension-and-multi-study-validation.md)); drafted, not yet implemented. It adds FR-168–FR-194, T168–T194 and V168–V194 and does not renumber any existing identity.

## Scope

Deliver FR-168–FR-194 for two journeys:
- **US10:** the analyst runs a declared Welch t-test engine with its paper-compatible rule labelled, enriches results against user-supplied offline pathway libraries, and reads calibration evidence for each classifier and resampling scheme the pipeline offers.
- **US11:** the Maintainer validates the pipeline across five design slots and one reference re-run slot, using a study manifest, with private studies kept local and public datasets fetched only after approval.

Normative statistics are SM42–SM46 in the [scientific contract](../001-downstream-proteomics/contracts/scientific-methods.md#methods-extension). Configuration, file formats, error codes and the study manifest are in [contracts/methods-extension.md](contracts/methods-extension.md).

Every module consumes frozen plan artifacts and completed outputs read-only and verifies their hashes. No module changes the primary matrix, the primary families, the default engine (limma) or any earlier result. This slice is additive: a configuration that does not declare `methods_extension` produces the same outputs as before.

### Generality and governance principles (normative for every packet in this slice)

1. **Declarative, study-agnostic.** No study, group, organism, cohort or file name appears in code. Engines, contrasts, libraries, classifiers, calibration settings and study designs are declared in configuration and validated against the frozen plan. There is no free-form code evaluation.
2. **Declared eligibility.** Every module publishes an eligibility rule over the frozen design and records ELIGIBLE or INAPPLICABLE with a typed reason before computation. Missing optional software is NOT_RUN, never INAPPLICABLE.
3. **Adapt, never hold up ([ADR 0010](../../docs/adr/0010-adaptive-post-de-policy.md)).** When a declared analysis cannot be honoured as stated, the invalid piece is refused as narrowly as possible, the valid alternative runs, and the adaptation is recorded with the analysis, item, requested value, used value and reason. A `required` declaration remains strict.
4. **Biological units.** Every test, resampling, permutation, cross-validation fold and calibration replicate operates on biological units (SM03), never on technical injections as if they were independent.
5. **Claims vocabulary.** Each output carries a frozen claim or method label. A label states what was computed. For example, a paper-compatible DEP set is labelled `not_multiplicity_controlled`, and a calibration result states that it applies only to the families and designs it tested.
6. **Confidentiality.** No private data, private identifiers, private file names, sample or subject names, or study numbers appear in the repository. Private results are committed only as non-identifying PASS/FAIL statuses ([SM46](../001-downstream-proteomics/contracts/scientific-methods.md#methods-extension)). All fixtures are synthetic.

## Requirements

### R15a — Welch t-test engine

- **FR-168 — Welch per-contrast test:** The system MUST compute, for each declared pairwise group contrast on log2 abundance, a two-sample Welch t statistic with unequal variances and Welch–Satterthwaite degrees of freedom. The effect MUST equal the difference of the two group means, and its confidence interval MUST be computed at the declared level.
- **FR-169 — Engine declaration and placement:** The system MUST run `welch` only when it is declared under `methods_extension.welch`, either as the single primary engine for a hypothesis or as a named secondary model with its own family. Without a declaration the default primary engine (limma, [ADR 0002](../../docs/adr/0002-primary-engine-limma.md)) and its outputs MUST be unchanged. A missing or failed Welch result MUST NOT be replaced by a limma result under the Welch name ([ADR 0005](../../docs/adr/0005-no-silent-method-fallback.md)).
- **FR-170 — Multiplicity and paper-compatible rule labelling:** The system MUST apply BH within the declared family and MAY additionally export a paper-compatible DEP rule that selects features with raw P below a declared cutoff. Every paper-compatible selection MUST carry the label `not_multiplicity_controlled`, and BH MUST be reported alongside it in the same table.
- **FR-171 — Paired variant:** When the frozen design declares paired or subject-blocked observations for a contrast, the system MUST run the paired t-test on within-subject differences with n − 1 degrees of freedom, in place of the two-sample Welch test for that contrast.
- **FR-172 — Adaptation for unhonourable declarations:** When a Welch request declares covariates or blocking that a Welch test cannot honour, the system MUST run only the valid unadjusted or within-subject form and record the adaptation (analysis, item, requested, used, reason) in the plan decision, the module's `eligibility.json` and a `W_WELCH_ADAPTED` warning. It MUST NOT drop the declaration silently.
- **FR-173 — Welch and limma on the same fixture:** For the same contrast and features, the system MUST report identical effects from Welch and limma and MUST keep their P values distinct, because limma moderates variances. Each engine's P value MUST match its own independent oracle.

### R15b — Offline pathway libraries

- **FR-174 — Hashed offline GMT libraries:** The system MUST read user-supplied local GMT library files only. Each library MUST be recorded in a manifest with its file SHA-256, source, version, species, identifier type, licence identifier and a redistributable flag. Analysis MUST NOT open a network connection or download a library ([SM14](../001-downstream-proteomics/contracts/scientific-methods.md), [SM44](../001-downstream-proteomics/contracts/scientific-methods.md#methods-extension)).
- **FR-175 — Licensing and redistribution rule:** KEGG and MSigDB library files MUST NOT be redistributed. They MUST NOT be committed to the repository, copied into an artifact or published by the release tooling. Reactome, GO and WikiPathways libraries MAY be used under their recorded open licences. The public repository MUST contain only synthetic, marked GMT fixtures.
- **FR-176 — Over-representation with the measured background:** The system MUST run hypergeometric over-representation for each declared library using as its universe the measured, mapped and eligible features of the analysed set (SM17). A whole-genome or library-wide background MUST NOT be used.
- **FR-177 — Pre-ranked enrichment on Welch statistics:** The system MUST support pre-ranked enrichment on a declared statistic, including the Welch t of FR-168. It MUST reuse the existing fgsea adapter with an exploratory label (FR-075) or document a named alternative, and its ranks MUST be tie-free under the SM15 tie-break rule.
- **FR-178 — Family discipline across libraries:** Each combination of library, method and contrast MUST be a separate family. BH MUST NOT be applied across libraries, and no pooled library family is created.

### R15c — Calibration extension

- **FR-179 — Penalised families under the null:** The system MUST calibrate lasso and elastic net classifiers, with the elastic-net penalty tuned by inner cross-validation, within the nested procedure of FR-153. Under the null, each family MUST meet the pre-registered null gates of FR-183.
- **FR-180 — SVM and random forest under the null and planted signal:** The system MUST calibrate linear and polynomial SVM classifiers, and the optional random forest classifier, under the null and under a planted signal. A missing random forest package MUST yield NOT_RUN for that family only.
- **FR-181 — LOOCV outer resampling:** The system MUST calibrate the nested procedure with leave-one-out outer resampling and inner 5-fold tuning, and report its pooled out-of-fold AUC under the same null gates as the repeated-stratified scheme.
- **FR-182 — Subject-blocked units:** The system MUST calibrate repeated-measures designs with subject-blocked folds and block-restricted permutations. The ungrouped reference MUST be shown to inflate the AUC, which demonstrates that the blocked form is required.
- **FR-183 — Pre-registered gates:** Before the first calibration run, the system MUST record its gate values, the seeds and their SHA-256 in the R15c evidence. Any later change to the gates MUST be detected. The pre-registered gates are: (a) null nested AUC, |mean − 0.5| ≤ 0.03 and mean + 1.96 SE ≤ 0.53; (b) null permutation P, binomial P(K ≥ k | rate 0.05) ≥ 0.01 and chi-square goodness-of-fit P over the B + 1 values ≥ 0.01; (c) leaky reference, mean AUC > 0.75; (d) subject-blocked nulls and the ungrouped reference, as in R14f; (e) planted signal (proposed, pending Maintainer confirmation before the first run), mean nested AUC ≥ 0.70 and permutation P ≤ 0.05 in at least 90 % of replicates. Gates (a)–(d) are carried over unchanged from R14f. The full table and the seeds are in the [contract](contracts/methods-extension.md#calibration-gates).
- **FR-184 — Calibration evidence and reproducibility:** Every calibration run MUST record its seeds, replicate counts, package versions, classifier configurations and gate evaluation. Replicate counts MUST NOT be reduced silently. Calibration results apply only to the families and designs they tested ([SM45](../001-downstream-proteomics/contracts/scientific-methods.md#methods-extension)).

### R15d — Multi-study validation framework

- **FR-185 — Study manifest:** The system MUST validate a versioned study manifest before any run. It records study identifier, origin (`private` or `public`), design type from the five slots or the reference slot, input format, organism, subject, batch and phenotype columns, the SHA-256 of every input file, the licence, and, for a public study, its approval reference.
- **FR-186 — Validate, run and verify per study:** For every study the system MUST perform three steps: `validate` (manifest and inputs), `run` (the declared pipeline) and `verify` (artifact hashes and stage states), followed by a design-specific acceptance check. Each step MUST report its own state.
- **FR-187 — Public datasets only after approval:** The system MUST download a public dataset only when the Maintainer has approved an accession shortlist, and the approval record matches the accession, URL and expected SHA-256. The approval, the URL and the observed SHA-256 MUST be recorded, and the data MUST be stored outside the repository. Analysis MUST NOT trigger a download.
- **FR-188 — Private studies stay local:** Private study inputs and outputs MUST remain on the local machine. The only private information that may be committed is a non-identifying status record containing the slot, design type, PASS or FAIL state and check counts. A committed status record that contains any sample, subject, file or group name MUST be refused.
- **FR-189 — Slot 1, label-free DDA:** A two-group label-free data-dependent-acquisition study in MaxQuant-style protein-group input MUST validate its required columns, apply the declared contamination and reverse-hit filters, and recover planted contrasts at the declared threshold.
- **FR-190 — Slot 2, paired or repeated measures:** A paired or repeated-measures study MUST preserve the subject structure throughout, run the paired variant (FR-171) or subject-blocked analysis, and refuse to treat repeated observations as independent.
- **FR-191 — Slot 3, TMT multiplex with batches:** A TMT study with plex or batch structure MUST check group-by-batch confounding before any batch adjustment. A confounded batch MUST be handled as an adaptation under ADR 0010 and recorded, and an aliased batch term MUST NOT be fitted.
- **FR-192 — Slot 4, DIA:** A DIA study in DIA-NN-style or Spectronaut-style report input MUST map every required column through a recorded mapping. An unmapped column MUST NOT be dropped silently, and the mapping hash MUST be recorded.
- **FR-193 — Slot 5, multi-group with continuous phenotype in a non-human organism:** A study with three or more groups and a continuous phenotype in a declared non-human organism MUST run the R14c association with complete-case handling and MUST map identifiers only through the declared organism resource ([SM14](../001-downstream-proteomics/contracts/scientific-methods.md)).
- **FR-194 — Reference re-run of a three-group serum study with Welch t:** The system MUST support an explicit reference slot that re-runs a three-group serum study with the Welch engine, the paper-compatible rule and BH reporting, and records only a non-identifying PASS/FAIL status when the Maintainer runs it locally. This slot is a sixth explicit slot, outside the five design slots, and its private run is Maintainer-only.

## Acceptance Scenarios

<a id="V168"></a>

### V168: Welch per-contrast test

**Fixture:** Two groups of 6 and 9 synthetic log2 features (200 features) with known means and unequal variances, including one feature with no difference.

**Oracle:** `stats::t.test(x, y, var.equal = FALSE)` called per feature in the test.

**Exact assertion:** t, Welch–Satterthwaite df, P, effect (mean(x) − mean(y)) and CI endpoints equal the oracle within 1e-10. The effect sign follows the declared numerator group.

**Negative case:** Pooled-variance output (`var.equal = TRUE`) or a reversed effect sign fails.

**Contract:** SM42; **owner:** R15a.

<a id="V169"></a>

### V169: Engine declaration and placement

**Fixture:** (a) A configuration with no `methods_extension` block. (b) Welch declared as a secondary model. (c) Welch declared as the single primary engine. (d) Two primary engines for one hypothesis. (e) Welch requested but not declared, with a failed Welch fit.

**Oracle:** The unchanged R05 limma outputs for (a); the single-primary rule (ADR 0002) for (d).

**Exact assertion:** (a) outputs are identical to the R05 baseline. (b) and (c) produce a Welch family labelled by its declared name. (c) makes Welch the only primary family for that hypothesis.

**Negative case:** (d) fails with the primary-count error. (e) emits no Welch result and the limma result is not relabelled as Welch (E_WELCH_NOT_DECLARED).

**Contract:** SM42; **owner:** R15a.

<a id="V170"></a>

### V170: Multiplicity and paper-compatible rule labelling

**Fixture:** 200 features in one family with planted up and down effects and a known set of raw P values below 0.05.

**Oracle:** `p.adjust(p, method = "BH")` and `p < 0.05` written in the test.

**Exact assertion:** BH q equals the oracle within 1e-12. The paper-compatible selection equals the raw-P set exactly. Every paper-compatible row carries the label `not_multiplicity_controlled`, and BH q is present in the same table.

**Negative case:** A paper-compatible selection without the label, or a q column computed from raw P and labelled as multiplicity-controlled, fails.

**Contract:** SM43; **owner:** R15a.

<a id="V171"></a>

### V171: Paired variant

**Fixture:** Eight subjects measured under two conditions, paired within subject, with planted within-subject differences.

**Oracle:** `stats::t.test(x, y, paired = TRUE)` per feature in the test.

**Exact assertion:** t, df (= 7), P, mean within-subject difference and CI equal the oracle within 1e-10.

**Negative case:** A two-sample Welch result emitted for the declared paired design fails (E_WELCH_PAIRED_UNIT).

**Contract:** SM42; **owner:** R15a.

<a id="V172"></a>

### V172: Adaptation for unhonourable declarations

**Fixture:** A balanced two-group design with a declared additive covariate (sex), requested for the Welch engine.

**Oracle:** The expected adaptation record (analysis `welch`, item `covariate:sex`, requested `adjusted`, used `unadjusted`, reason `welch_not_covariate_adjusting`) written in the test.

**Exact assertion:** The Welch result equals the unadjusted two-sample Welch oracle. The record appears in the plan decision, `eligibility.json` and the `W_WELCH_ADAPTED` warning.

**Negative case:** A Welch result with no adaptation record, or a run that drops the covariate without any record, fails.

**Contract:** SM42, ADR 0010; **owner:** R15a.

<a id="V173"></a>

### V173: Welch and limma on the same fixture

**Fixture:** Two groups of 12 with heterogeneous variances on 200 features, one planted effect.

**Oracle:** Effects from `limma::lmFit` with design `~ group` and no moderation, and the moderated P from the frozen limma settings, both written in the test. The Welch P values come from the V168 oracle.

**Exact assertion:** Effects from Welch and limma are identical within 1e-10. Each engine's P values equal its own oracle. At least one feature has Welch and limma P values that differ by more than 1e-6.

**Negative case:** Identical P values for both engines on this heterogeneous fixture fail, because that would indicate that Welch was implemented by the pooled limma model.

**Contract:** SM42; **owner:** R15a.

<a id="V174"></a>

### V174: Hashed offline GMT libraries

**Fixture:** A synthetic GMT library with three marked sets and its manifest record.

**Oracle:** The SHA-256 of the file bytes computed in the test with `hashlib`.

**Exact assertion:** The manifest records the hash, source, version, species, identifier type, licence identifier and redistributable flag. The run reads the file and records the same hash.

**Negative case:** A one-byte edit fails E_RESOURCE_HASH. A socket call during the run fails E_RESOURCE_DOWNLOAD (the test patches the socket layer).

**Contract:** SM44 (with SM14); **owner:** R15b.

<a id="V175"></a>

### V175: Licensing and redistribution rule

**Fixture:** The repository tree with its GMT fixtures, plus variants of a release artifact containing (a) a KEGG-licensed library, (b) an MSigDB-licensed library and (c) a non-synthetic GMT with redistributable set to true.

**Oracle:** The licence table written in the test (KEGG and MSigDB non-redistributable; open-licence libraries redistributable only with a recorded identifier).

**Exact assertion:** Every GMT file in the repository carries the synthetic marker. The release artifact in the valid case contains no KEGG or MSigDB file.

**Negative case:** Each variant (a), (b) and (c) fails E_LIBRARY_REDISTRIBUTION and no artifact is written.

**Contract:** SM44; **owner:** R15b.

<a id="V176"></a>

### V176: Over-representation with the measured background

**Fixture:** A measured universe of 500 mapped features, a library set of 30 features with 8 foreground features, and a second set with k = 0 overlap.

**Oracle:** `stats::phyper(k − 1, K, N − K, n, lower.tail = FALSE)` written in the test, with N the measured universe, K the set size within the universe, and n the foreground size.

**Exact assertion:** P values equal the oracle within 1e-12 relative. BH is applied within the family. The k = 0 set is retained through BH.

**Negative case:** A whole-genome background (N = 20,000) fails E_ORA_UNIVERSE.

**Contract:** SM17, SM44; **owner:** R15b.

<a id="V177"></a>

### V177: Pre-ranked enrichment on Welch statistics

**Fixture:** Ten tie-free ranks taken from the V168 Welch statistics and two sets of size 3 and 4.

**Oracle:** A brute-force running-sum enrichment score (weights |t|, gene-set weight 1) and an exact permutation P, from enumeration over all subsets of the set size, written in the test.

**Exact assertion:** ES and exact P equal the oracle within 1e-12. The output carries the exploratory label.

**Negative case:** A duplicated feature ID fails E_RANK_DUPLICATE. Ties not broken by the SM15 rule fail E_RANK_TIE.

**Contract:** SM44, SM16; **owner:** R15b.

<a id="V178"></a>

### V178: Family discipline across libraries

**Fixture:** Two libraries with overlapping sets and two contrasts.

**Oracle:** The family membership computed in the test (one family per library, method and contrast).

**Exact assertion:** Family IDs are distinct and each q is computed within its family. No pooled family exists.

**Negative case:** A BH adjustment spanning libraries fails E_ENRICHMENT_FAMILY_POOLED.

**Contract:** SM44, SM16; **owner:** R15b.

<a id="V179"></a>

### V179: Penalised families under the null

**Fixture:** Pure-noise design with 15 versus 15 biological units and 200 features, nested 2 × 5-fold outer CV with 3-fold inner tuning, the R14f seeds.

**Oracle:** A direct nested loop written in the test with `glmnet::cv.glmnet` (lasso, alpha = 1; elastic net, alpha = 0.5) on the same folds and seeds.

**Exact assertion:** Per-replicate AUC equals the oracle within 1e-10 for the first 10 replicates. The mean AUC meets the null gates of FR-183.

**Negative case:** Lambda chosen on the outer test folds fails E_BIOMARKER_LEAKAGE. A missing glmnet yields NOT_RUN.

**Contract:** SM45, SM35; **owner:** R15c.

<a id="V180"></a>

### V180: SVM and random forest under the null and planted signal

**Fixture:** The V179 null design, and a planted-signal design of the same size with five features shifted by 1.5 standard deviations.

**Oracle:** Direct `e1071::svm` and `randomForest::randomForest` fits in the test, with AUC from decision values oriented by the model class coding.

**Exact assertion:** Null gates hold for both SVM kernels. The planted-signal gates of FR-183 hold. SVM AUC equals the oracle within 1e-10. The random forest is NOT_RUN with E_DEPENDENCY_UNAVAILABLE when randomForest is absent.

**Negative case:** AUC computed from predicted classes instead of decision values fails.

**Contract:** SM45, SM35; **owner:** R15c.

<a id="V181"></a>

### V181: LOOCV outer resampling

**Fixture:** Twelve versus twelve null units with leave-one-out outer resampling and 5-fold inner tuning.

**Oracle:** A direct leave-one-out loop written in the test with the same selector and classifier.

**Exact assertion:** The pooled out-of-fold AUC equals the oracle. The mean across replicates meets the null gates.

**Negative case:** A selector fitted on all units before the leave-one-out loop fails E_BIOMARKER_LEAKAGE and is reported as the leaky reference.

**Contract:** SM45, SM35; **owner:** R15c.

<a id="V182"></a>

### V182: Subject-blocked calibration

**Fixture:** Fifteen subjects with two observations each, pure noise and a planted signal, with subject-blocked folds and block-restricted permutations.

**Oracle:** Grouped fold membership checks in the test (no subject appears in both training and test of any fold), and the ungrouped reference computed from the same data.

**Exact assertion:** The blocked null gates of R14f hold. The blocked permutation P is uniform. The ungrouped reference exceeds the grouped mean AUC by more than 0.05.

**Negative case:** A declaration with `group_by_subject: false` in a blocked design runs the grouped form with a recorded adaptation (ADR 0010), and the ungrouped reference is labelled leaky. A run that uses ungrouped folds without that record fails E_BIOMARKER_GROUP_LEAKAGE.

**Contract:** SM45, SM35; **owner:** R15c.

<a id="V183"></a>

### V183: Pre-registered gates

**Fixture:** The gate file and its SHA-256 record written before the first run, and a run summary.

**Oracle:** The SHA-256 of the gate file recomputed in the test, and the gate evaluation recomputed from the summary.

**Exact assertion:** The recomputed hash equals the recorded hash, and the summary's pass or fail decisions equal the recomputed decisions.

**Negative case:** Any change to a gate value after the run changes the hash and fails E_CALIBRATION_GATE_CHANGED.

**Contract:** SM45; **owner:** R15c.

<a id="V184"></a>

### V184: Calibration evidence and reproducibility

**Fixture:** A calibration run with recorded seeds, replicate counts and package versions for each family.

**Oracle:** A re-run with the same seeds, reproduced in the test for one family.

**Exact assertion:** The re-run reproduces the raw AUC values bit-for-bit. Every family has its seeds, declared replicate count, versions and classifier configuration recorded.

**Negative case:** A run with fewer replicates than declared and no recorded adaptation fails E_CALIBRATION_REPLICATES.

**Contract:** SM45; **owner:** R15c.

<a id="V185"></a>

### V185: Study manifest

**Fixture:** A valid synthetic manifest with two input files and matching SHA-256 values, plus variants with a missing hash, an unknown design type, and a public study without an approval reference.

**Oracle:** The manifest schema and the file hashes computed in the test.

**Exact assertion:** The valid manifest loads with every required field recorded.

**Negative case:** A missing hash fails E_MANIFEST_HASH. An unknown design type fails E_MANIFEST_DESIGN. A public study without an approval reference fails E_PUBLIC_DATA_UNAPPROVED.

**Contract:** SM46; **owner:** R15d.

<a id="V186"></a>

### V186: Validate, run and verify per study

**Fixture:** One small synthetic study for each of the five design slots and the reference slot.

**Oracle:** The stage states and the design-specific acceptance checks computed in the test from the synthetic planted structure.

**Exact assertion:** Every study reaches PASS for validate, run, verify and its design acceptance check.

**Negative case:** Editing an artifact after the run changes verify to FAIL with E_ARTIFACT_HASH.

**Contract:** SM46; **owner:** R15d.

<a id="V187"></a>

### V187: Public datasets only after approval

**Fixture:** A synthetic accession shortlist, a local file location standing in for a public source, and an approval record. The test patches the network layer.

**Oracle:** The SHA-256 of the local file computed in the test.

**Exact assertion:** The download occurs only when the approval record names the same accession, URL and expected SHA-256. The observed SHA-256 equals the oracle and the data are written outside the repository root.

**Negative case:** A download without an approval record fails E_PUBLIC_DATA_UNAPPROVED and writes no file. An approval for a different accession fails the same way.

**Contract:** SM46; **owner:** R15d.

<a id="V188"></a>

### V188: Private studies stay local

**Fixture:** A private-style output tree with synthetic sample and subject names, and the committed status record.

**Oracle:** The allowed field list of the status record (slot, design type, state, check counts) written in the test.

**Exact assertion:** The committed status record contains only allowed fields, and the repository scan finds no sample, subject, file or group name from the private tree.

**Negative case:** A status row containing a sample name fails E_PRIVATE_IDENTIFIER.

**Contract:** SM46; **owner:** R15d.

<a id="V189"></a>

### V189: Slot 1, label-free DDA

**Fixture:** A synthetic MaxQuant-style protein-group table with two groups of 10 units, reverse-hit and contaminant rows, and planted contrasts.

**Oracle:** The declared filter applied in the test, and the planted feature list.

**Exact assertion:** The parsed matrix has the expected rows after filtering and n = 10 per group. Planted features are recovered at the declared threshold.

**Negative case:** A missing required intensity column fails intake with its typed error and no run is produced.

**Contract:** SM46 (with SM03); **owner:** R15d.

<a id="V190"></a>

### V190: Slot 2, paired or repeated measures

**Fixture:** Eight subjects at two time points (paired) and a repeated variant with three time points per subject.

**Oracle:** The paired t-test of V171 and the blocked model fit written in the test.

**Exact assertion:** The subject identifier is preserved on every row. Paired results equal the oracle.

**Negative case:** Treating the repeated observations as independent (unpaired Welch) fails E_WELCH_PAIRED_UNIT.

**Contract:** SM46 (with SM03); **owner:** R15d.

<a id="V191"></a>

### V191: Slot 3, TMT multiplex with batches

**Fixture:** Synthetic TMT channels across three plex batches, a balanced design and a variant in which one group lies entirely in one batch.

**Oracle:** The group-by-batch cross-tabulation computed in the test.

**Exact assertion:** The balanced design runs with the batch covariate. The confounded variant runs unadjusted with a recorded ADR 0010 adaptation.

**Negative case:** A batch-adjusted result for the confounded variant fails E_CONTRAST_NONESTIMABLE.

**Contract:** SM46 (with SM07, ADR 0010); **owner:** R15d.

<a id="V192"></a>

### V192: Slot 4, DIA

**Fixture:** Synthetic DIA-NN-style and Spectronaut-style report tables with the same underlying values.

**Oracle:** The column-to-canonical mapping table written in the test.

**Exact assertion:** Both mapped matrices equal the oracle and the mapping hash is recorded. The mapping is qualified only on these synthetic exports.

**Negative case:** An unmapped required column that is silently dropped fails E_MAPPING_UNMAPPED.

**Contract:** SM46; **owner:** R15d.

<a id="V193"></a>

### V193: Slot 5, multi-group with continuous phenotype in a non-human organism

**Fixture:** Four groups of 10 units, a continuous phenotype with 3 missing values, and a declared non-human organism resource.

**Oracle:** A direct limma fit with `~ group + phenotype` on complete cases, written in the test.

**Exact assertion:** Phenotype coefficients, moderated t and P equal the oracle, n = 37 complete cases is reported, and identifiers are mapped only through the declared organism resource.

**Negative case:** An imputed phenotype fails E_PHENOTYPE_IMPUTATION. A human-only mapping used for the non-human organism fails E_RESOURCE_TAXONOMY.

**Contract:** SM46 (with SM33, SM14); **owner:** R15d.

<a id="V194"></a>

### V194: Reference re-run of a three-group serum study with Welch t

**Fixture:** A synthetic three-group study generated by the synthetic design factory, with the same design shape as the reference analysis and no study data.

**Oracle:** Per-contrast Welch statistics, the paper-compatible raw-P selection and BH computed in the test.

**Exact assertion:** Welch effects, CIs and P values equal the oracle. The paper-compatible set equals the oracle and is labelled `not_multiplicity_controlled`. BH q is reported alongside it.

**Negative case:** A paper-compatible set labelled as multiplicity-controlled fails. A private identifier in any output fails E_PRIVATE_IDENTIFIER. The Maintainer's local private re-run is NOT_RUN in automated tests and is recorded only as a non-identifying PASS/FAIL status.

**Contract:** SM42, SM43, SM46; **owner:** R15d.
