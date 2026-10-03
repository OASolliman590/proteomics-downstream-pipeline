# Feature Specification: Multivariate PERMANOVA and dispersion testing

**Phase:** 2 (v0.2-qualified-methods). **Packet:** R13. **Status:** operator-authorized scope amendment of 2026-10-01 ([ADR 0006](../../docs/adr/0006-permanova-scope-amendment.md)); adds FR-121–FR-130, T121–T130 and V121–V130 without renumbering any existing identity.

## Scope

Deliver only FR-121–FR-130 for US7 (multivariate group separation as a descriptive-inferential complement to protein-wise limma). Normative rules are SM26–SM30 in the [scientific contract](../001-downstream-proteomics/contracts/scientific-methods.md#multivariate-permanova). Configuration and outputs are in [contracts/permanova.md](contracts/permanova.md). Ownership and dispatch are in the [packet index](../001-downstream-proteomics/packet-index.md) and [ownership record](../001-downstream-proteomics/packet-ownership.json). PERMANOVA consumes the frozen plan's primary matrix, observed mask and design grouping (after R05) and, for DEP-derived feature sets, the completed R05 differential table; it never changes them.

## Requirements

- **FR-121 — Distance and feature-set construction:** The system MUST compute the declared distance (default Euclidean on per-protein z-scored log2 abundance, equivalent to PCA space) on genuinely observed, complete-case features of the analysed observations, for all-protein, declared-panel and DEP-derived feature sets, recording every excluded feature and reason.
- **FR-122 — PERMANOVA term reporting:** The system MUST run vegan::adonis2 with a declared permutation count (default 9,999) and a recorded seed and report Df, SS, R², pseudo-F and P, with P floored at 1/(nperm+1) and displayed as "< floor" at the floor.
- **FR-123 — Global and pairwise tests:** The system MUST test the global group term across all declared groups and every group pair post hoc, adjusting pairwise P values within each feature set by the declared method (Holm default).
- **FR-124 — Covariate-adjusted margins with restricted permutations:** The system MUST use marginal sums of squares, test covariate terms with permutations restricted within the primary grouping factor, test the group term only from an unrestricted model, and never report a group term from a blocked model.
- **FR-125 — Interaction estimability:** The system MUST fit a group×covariate interaction only when every cell is non-empty and refuse otherwise with a typed error.
- **FR-126 — PERMDISP alongside every PERMANOVA:** The system MUST report vegan::betadisper + permutest with the same permutation scheme next to each group/covariate test and state the interpretation rule.
- **FR-127 — Per-feature contribution and mean-R² identity:** The system MUST report each feature's univariate group R² and, under z-scored Euclidean distance, assert that the multivariate R² equals their mean.
- **FR-128 — Selection-circularity context:** The system MUST compare panels selected on the same data, and DEP-derived sets, with a null of random equal-size feature sets and with the best possible same-size set, and flag in-sample R² as optimistic.
- **FR-129 — Outputs, figures and report integration:** The system MUST export machine-readable tables, figure source tables and figures (ordination with centroids and spider lines, PERMDISP distances, permutation nulls, R² bars, per-feature R², random-set null) and show them in the offline report.
- **FR-130 — Synthetic qualification fixtures and reproducibility:** The system MUST behave as specified on null, strong location-shift and dispersion-only synthetic fixtures and reproduce identical results for the same seed.

## Acceptance Scenarios

<a id="V121"></a>

### V121: Distance and feature-set construction

**Fixture:** Six-feature, three-group synthetic matrix with one feature that has a missing cell, one constant feature, a declared two-feature panel and a DEP-derived set from a completed R05 table.

**Oracle:** Base-R `dist(scale(t(X)))` on the complete, non-constant features, and Euclidean distance between PCA scores of the same z-scored matrix.

**Exact assertion:** Distances equal the oracle within 1e-10 and equal the PCA-space distances; the incomplete and constant features are excluded with reasons; panel and DEP-derived sets contain exactly their declared/selected members; a configurable metric (Manhattan) matches `vegdist`.

**Negative case:** A panel naming an unknown feature fails E_PERMANOVA_FEATURE_SET; a set with fewer than two usable features is refused rather than reported.

**Contract:** SM26; **owner:** R13.

<a id="V122"></a>

### V122: PERMANOVA term reporting

**Fixture:** Two-group 6+6 synthetic matrix.

**Oracle:** Hand-computed between/within sums of squares from the Gower-centred distance matrix, pseudo-F = (SS_between/df_between)/(SS_within/df_within), and a direct `adonis2` call with the same seed.

**Exact assertion:** Df, SS, R², pseudo-F match the oracle within 1e-10; P equals the direct call; nperm and seed are recorded; a P at 1/(nperm+1) is displayed "< 1/(nperm+1)".

**Negative case:** A production profile with fewer than 999 permutations fails E_PERMANOVA_RESOLUTION; P is never displayed as 0.

**Contract:** SM27; **owner:** R13.

<a id="V123"></a>

### V123: Global and pairwise tests

**Fixture:** Three-group synthetic matrix with one strongly shifted group.

**Oracle:** Direct `adonis2` on all groups and on each pair subset (re-scaled within the subset), then `p.adjust(..., "holm")` over the pairwise P values of a feature set.

**Exact assertion:** One global row and choose(k,2) pairwise rows per feature set; P and adjusted P equal the oracle; the adjustment method is recorded and configurable (BH also tested).

**Negative case:** Adjusting pairwise P values across feature sets, or reporting the global row as adjusted, fails.

**Contract:** SM27; **owner:** R13.

<a id="V124"></a>

### V124: Covariate-adjusted margins and restricted permutations

**Fixture:** Two groups crossed with a two-level covariate.

**Oracle:** `adonis2(D ~ group + cov, by = "margin")` unrestricted for the group term and `permutations = how(blocks = group)` for the covariate term, with the same seed; `shuffleSet` under the blocked design.

**Exact assertion:** The group row comes from the unrestricted model and the covariate row from the blocked model; every blocked permutation keeps each observation inside its group; no group row from the blocked model appears in any output.

**Negative case:** Reporting the group term from the blocked model, or permuting the covariate freely, fails.

**Contract:** SM28; **owner:** R13.

<a id="V125"></a>

### V125: Interaction estimability

**Fixture:** Groups × covariate with all cells filled, and a variant where one group has no level "F".

**Oracle:** Cell counts by group × level.

**Exact assertion:** The filled design reports the interaction term from the margin model; the empty-cell variant is refused.

**Negative case:** The empty-cell variant fails E_PERMANOVA_INTERACTION_NONESTIMABLE before any permutation; no interaction row is written.

**Contract:** SM28; **owner:** R13.

<a id="V126"></a>

### V126: PERMDISP and interpretation

**Fixture:** Dispersion-only synthetic matrix (same centroid, group B four-fold spread) and a location-shift matrix.

**Oracle:** Direct `permutest(betadisper(D, g))` with the same seed and permutation scheme.

**Exact assertion:** PERMDISP F/P match the oracle; the dispersion-only fixture has significant PERMDISP and is labelled `dispersion_difference_location_not_established`; the location fixture is labelled `location_shift`; the report states the interpretation rule.

**Negative case:** A PERMANOVA row without its PERMDISP result, or a `location_shift` label when PERMDISP is significant, fails.

**Contract:** SM29; **owner:** R13.

<a id="V127"></a>

### V127: Per-feature contribution and mean-R² identity

**Fixture:** Strong location-shift fixture.

**Oracle:** For each feature, `anova(lm(z ~ group))` SS_group/SS_total on the z-scored feature.

**Exact assertion:** Per-feature R² match the oracle within 1e-10 and their mean equals the global PERMANOVA R² within 1e-8 under z-scored Euclidean distance; the identity check is recorded as not_applicable for other metrics.

**Negative case:** A violated identity fails the stage with E_PERMANOVA_IDENTITY rather than being reported.

**Contract:** SM30; **owner:** R13.

<a id="V128"></a>

### V128: Selection-circularity context

**Fixture:** Location-shift matrix with a same-data panel of the most shifted features and an independent panel.

**Oracle:** Mean of univariate R² over random equal-size subsets drawn with the recorded seed, and the mean of the top-k univariate R².

**Exact assertion:** The null has the declared number of draws (default 2,000) and reproduces with the same seed; the best-possible R² equals the top-k mean; same-data and DEP-derived sets are flagged "in-sample R² is optimistic"; an independent panel is not flagged.

**Negative case:** A same-data panel reported without the null context, or labelled as independent validation, fails.

**Contract:** SM30; **owner:** R13.

<a id="V129"></a>

### V129: Outputs, figures and report integration

**Fixture:** Full synthetic run with PERMANOVA enabled.

**Oracle:** Figure source tables compared with the test tables; PNG signatures; report HTML references.

**Exact assertion:** Tables and figures exist with hashes in the stage manifest; each figure has a source table; the offline report shows the PERMANOVA section with tests, PERMDISP, interpretation rule and circularity flags and links only local files.

**Negative case:** A figure without a source table, a remote resource, or a PERMANOVA section shown as completed when the stage failed or did not run fails.

**Contract:** SM30; **owner:** R13.

<a id="V130"></a>

### V130: Synthetic qualification and reproducibility

**Fixture:** Null (no group structure), strong location-shift and dispersion-only matrices with known construction.

**Oracle:** Known generating structure; repeated runs with the same seed.

**Exact assertion:** Null: global P > 0.05 and R² below the location fixture's; location: global P at the permutation floor and large R²; dispersion-only: PERMDISP P < 0.05; identical seeds give byte-identical test tables.

**Negative case:** A different seed changing the documented P only within Monte Carlo resolution is allowed; identical seeds with different outputs fail.

**Contract:** SM27, SM29; **owner:** R13.

## Boundary

PERMANOVA describes multivariate separation of the analysed observations; it does not replace protein-wise inference, is not a classifier or validation of a panel, and makes no causal claim. In-sample selected sets are optimistic by construction. Private study data and private-derived numbers are never fixtures. Missing vegan/permute is NOT_RUN, not scientific inapplicability.
