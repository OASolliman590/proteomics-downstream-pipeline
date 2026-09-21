# Feature Specification: Preprocessing, missingness and quality control

**Phase:** 1. **Packet:** R03. **Status:** 1.2.0-frozen; implementation pending prerequisite acceptance and explicit authorization.

## Scope

Deliver only FR-021–FR-030 for US1. [Packet index](../001-downstream-proteomics/packet-index.md), [ownership](../001-downstream-proteomics/packet-ownership.json), [data model](../001-downstream-proteomics/data-model.md) and [scientific contract](../001-downstream-proteomics/contracts/scientific-methods.md) define the shared interfaces. No implementation dispatch is authorized yet.

## Requirements

- **FR-021 — Preserve processed abundance by default:** The system MUST leave the primary preserve input unchanged and issue a distinct output artifact for any requested transformation.
- **FR-022 — LFQ normalization policies:** The system MUST match analytic median/reference factors, preserve NA, and keep quantile normalization explicitly sensitivity-only.
- **FR-023 — TMT plex and bridge policies:** The system MUST recover the declared normalization/estimand in both qualified strategies, preserving reference-channel provenance.
- **FR-024 — Contrast-aware coverage masks:** The system MUST with minimum two and fraction .5, mark 4/4 and 2/4 eligible, 1/4 ineligible and 0/4 nonestimable for ordinary available-case inference; retain all reasons.
- **FR-025 — Missingness mechanism diagnostics:** The system MUST distinguish observed coverage, prior fill and numeric availability; show group/observation/abundance summaries without asserting MCAR/MAR/MNAR from a cutoff.
- **FR-026 — PCA and correlation diagnostics:** The system MUST match eigenvalues/explained variance and aligned scores up to sign; retain pairwise n; constant/empty inputs yield explanatory QC states without fake variance.
- **FR-027 — Prespecified exclusions and watchlists:** The system MUST flag the extreme sample without dropping it; honor only declared exclusions and change the plan hash when policy/reason changes.
- **FR-028 — Explicit imputation sensitivities:** The system MUST use exact algorithm names and secondary models; deterministic replacement equals the observed feature minimum, seeded Gaussian repeats, and primary values/masks remain unchanged.
- **FR-029 — Detection-only exploratory endpoint:** The system MUST match the exact detection P values and separate detection family; never substitute that endpoint for continuous abundance inference.
- **FR-030 — QC source data and stage report:** The system MUST complete the declared QC-only scope without any model fit or discovery requirement; every displayed numerical QC value has an exact source table.

## Acceptance Scenarios

<a id="V021"></a>

### V021: Preserve processed abundance by default

**Fixture:** The normalized log2 example with preserve/none and one separately requested median normalization.

**Oracle:** Exact input array/mask equality, independently calculated transformation hash.

**Exact assertion:** Leave the primary preserve input unchanged and issue a distinct output artifact for any requested transformation.

**Negative case:** A second log or hidden normalization fails its guard; primary data never become the PCA display matrix.

**Contract:** SM02; **owner:** R03; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V022"></a>

### V022: LFQ normalization policies

**Fixture:** Two columns offset by +2 log2 with complete references; then one column lacking required reference coverage.

**Oracle:** Observed median subtraction relative to median of column medians; declared reference medians.

**Exact assertion:** Match analytic median/reference factors, preserve NA, and keep quantile normalization explicitly sensitivity-only.

**Negative case:** Missing required reference coverage raises E_REFERENCE_COVERAGE; do not silently use all proteins instead.

**Contract:** SM04; **owner:** R03; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V023"></a>

### V023: TMT plex and bridge policies

**Fixture:** Two tiny plexes with known loading offsets and bridge labels; a separate balanced no-bridge design.

**Oracle:** Hand-calculated within-plex offsets and full-rank within-plex effect design.

**Exact assertion:** Recover the declared normalization/estimand in both qualified strategies, preserving reference-channel provenance.

**Negative case:** Missing bridge in bridge mode fails E_TMT_BRIDGE_REQUIRED; plex=treatment confounding fails rather than being corrected away.

**Contract:** SM04; **owner:** R03; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V024"></a>

### V024: Contrast-aware coverage masks

**Fixture:** Eight features with observed counts 4/4, 2/4, 1/4 and 0/4 in a required group, plus prior-imputed cells.

**Oracle:** Count genuine observed biological units from the original mask, not numeric availability.

**Exact assertion:** With minimum two and fraction .5, mark 4/4 and 2/4 eligible, 1/4 ineligible and 0/4 nonestimable for ordinary available-case inference; retain all reasons.

**Negative case:** An all-missing target group is not filled for limma. A native-dropout request is deferred to the separately qualified R06 route, not prematurely credited in Phase 1.

**Contract:** SM04; **owner:** R03; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V025"></a>

### V025: Missingness mechanism diagnostics

**Fixture:** Same numeric matrix with none_documented versus masked prior imputation, then unknown original mask.

**Oracle:** Explicit original-observed/imputed/numeric masks and missingness tallies.

**Exact assertion:** Distinguish observed coverage, prior fill and numeric availability; show group/observation/abundance summaries without asserting MCAR/MAR/MNAR from a cutoff.

**Negative case:** Unknown original mask blocks primary observed-coverage inference with E_ORIGINAL_MASK_REQUIRED; it cannot become all true.

**Contract:** SM04; **owner:** R03; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V026"></a>

### V026: PCA and correlation diagnostics

**Fixture:** A nonconstant 4×4 matrix, its one-NA variant, an all-constant matrix and zero retained features.

**Oracle:** Independent SVD after exactly declared display centering/fill; pairwise correlation on shared observed cells.

**Exact assertion:** Match eigenvalues/explained variance and aligned scores up to sign; retain pairwise n; constant/empty inputs yield explanatory QC states without fake variance.

**Negative case:** Hash primary input before/after PCA: unchanged; feeding display-fill artifact to the model is rejected.

**Contract:** SM06; **owner:** R03; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V027"></a>

### V027: Prespecified exclusions and watchlists

**Fixture:** One isolated extreme sample plus a config with an explicit exclusion/reason and a changed-reason replan.

**Oracle:** Retained observation IDs plus R01 `canonical_json_sha256` applied to two otherwise identical complete synthetic plan envelopes containing the R03-resolved exclusion fragment. R04/V039 later repeats the check on the actual AnalysisPlan before Phase 1 advancement.

**Exact assertion:** Flag the extreme sample without dropping it; honor only declared exclusions. R03 emits a canonical resolved exclusion/watchlist fragment and proves with the R01 canonical hash utility that changing policy/reason changes a complete plan-envelope hash. R03 does not claim to create the actual AnalysisPlan; after R04, V027 is rerun with V039 and the real plan hash must change.

**Negative case:** An exclusion without a nonblank reason fails E_EXCLUSION_REASON_REQUIRED; a post-fit exclusion cannot mutate the existing plan. Treating the R03 fragment hash as the final plan hash, or advancing Phase 1 without the R04/V039 integration rerun, fails.

**Contract:** SM06; **owner:** R03; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V028"></a>

### V028: Explicit imputation sensitivities

**Fixture:** Tiny matrix with internal NA, a known feature minimum and enough features for configured KNN; two identical seeded Gaussian runs.

**Oracle:** SM05 algorithms; direct pinned impute::impute.knn invocation for the named KNN configuration.

**Exact assertion:** Use exact algorithm names and secondary models; deterministic replacement equals the observed feature minimum, seeded Gaussian repeats, and primary values/masks remain unchanged.

**Negative case:** Name MinDet with stochastic parameters is rejected; zero SD Gaussian fill fails E_IMPUTATION_DISTRIBUTION instead of arbitrary replacement.

**Contract:** SM05; **owner:** R03; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V029"></a>

### V029: Detection-only exploratory endpoint

**Fixture:** Independent binary counts 3 observed/1 missing versus 1 observed/3 missing; two features in a declared detection family.

**Oracle:** Two-sided Fisher exact enumeration of fixed-margin 2×2 tables and independently authored BH. The production endpoint is the R03-owned detection-only Fisher/BH implementation in `detection.R`, not the later R05 general family-adjustment API.

**Exact assertion:** Match the exact detection P values and detection-only BH values in a separate declared family; never substitute that endpoint for continuous abundance inference. R05 later reruns V029 as an integration regression and must not merge detection P values into a protein abundance family.

**Negative case:** Paired or repeated-design detection request yields E_DETECTION_DESIGN_UNSUPPORTED with no detection P column. Calling a not-yet-existing R05 adjustment service, or merging detection and abundance P values during later integration, fails.

**Contract:** SM04; **owner:** R03; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V030"></a>

### V030: QC source data and stage report

**Fixture:** QC-only config, zero retained variable features and a second ordinary eight-feature example.

**Oracle:** Published QC dimensions/tables compared with the plotted input data.

**Exact assertion:** Complete the declared QC-only scope without any model fit or discovery requirement; every displayed numerical QC value has an exact source table.

**Negative case:** An unavailable QC computation is a typed NOT_RUN/FAILED diagnostic, not a fabricated zero-valued PCA panel.

**Contract:** SM06; **owner:** R03; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

## Boundary

No recovered source/audit changes, private-data transfer, invented license, fake backend or unsupported completion claim. A missing prerequisite is NOT_RUN, not a successful acceptance case. Only the Maintainer updates traceability after independent gate execution.
