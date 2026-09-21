# Feature Specification: Assay-qualified DEqMS and proDA backends

**Phase:** 2. **Packet:** R06. **Status:** 1.2.0-frozen; implementation pending Phase 1 acceptance and separate authorization.

## Scope

Deliver only FR-051–FR-060 for US2. [Packet index](../001-downstream-proteomics/packet-index.md), [ownership](../001-downstream-proteomics/packet-ownership.json), [data model](../001-downstream-proteomics/data-model.md) and [scientific contract](../001-downstream-proteomics/contracts/scientific-methods.md) define the shared interfaces. No implementation dispatch is authorized yet.

## Requirements

- **FR-051 — Count evidence schema:** The system MUST accept genuine declared count evidence only; preserve original zeros/missing values and require positive actual fit counts after any explicit justified policy.
- **FR-052 — DEqMS fit pipeline:** The system MUST match count-adjusted statistics and P values at frozen backend tolerance, with actual call/version evidence.
- **FR-053 — DEqMS statistic alignment:** The system MUST align counts by identity and retain sca.P.Value/statistic/variance distinctly from original limma fields; compute central q from sca raw P.
- **FR-054 — DEqMS hypotheses and design guard:** The system MUST qualify the sparse exact route numerically and reject unqualified block/precision-weight/TREAT combinations before backend execution.
- **FR-055 — LFQ dropout provenance gate:** The system MUST permit only the qualified LFQ independent fixed-model path; preserve actual original observations and native missingness.
- **FR-056 — proDA fitting and test_diff:** The system MUST match native effects, P and native uncertainty/diagnostics where reference results exist; retain n_obs and prior/dropout dependence for all-missing-group estimates.
- **FR-057 — Alternative-engine uncertainty:** The system MUST export only justified intervals/statistics with correct native names/df; unsupported quantities stay unavailable with reasons.
- **FR-058 — Method sensitivity comparison:** The system MUST compare matched effects and disclose unmatched universes; keep model/hypothesis families distinct and never promote the larger hit count to primary.
- **FR-059 — Dependencies and failure propagation:** The system MUST distinguish scientific inapplicability from missing installation and numerical failure; required failure blocks the run, optional operational failure is PARTIAL/nonzero.
- **FR-060 — Assay-engine golden/calibration fixtures:** The system MUST prespecify reference calls/expected statuses and save actual diagnostics, precision and eligibility outcomes; no tiny-fixture claim of calibrated FDR.

## Acceptance Scenarios

<a id="V051"></a>

### V051: Count evidence schema

**Fixture:** Eight protein IDs with peptide counts [1,2,3,4,5,6,7,8], a separate per-plex PSM table and a provenance manifest.

**Oracle:** Explicit source grain/type/aggregation and keyed feature alignment.

**Exact assertion:** Accept genuine declared count evidence only; preserve original zeros/missing values and require positive actual fit counts after any explicit justified policy.

**Negative case:** Observed-sample counts or abundance-derived proxies fail E_DEQMS_COUNT_EVIDENCE; hidden +1 is not an acceptable correction.

**Contract:** SM10; **owner:** R06; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V052"></a>

### V052: DEqMS fit pipeline

**Fixture:** ≤20 proteins with aligned strictly positive support counts and a complete independent log2 design.

**Oracle:** Separate direct official lmFit/eBayes/count/spectraCounteBayes/outputResult invocation in the pinned environment.

**Exact assertion:** Match count-adjusted statistics and P values at frozen backend tolerance, with actual call/version evidence.

**Negative case:** An ordinary eBayes table copied as DEqMS or a spectraCounteBayes call without aligned count evidence fails.

**Contract:** SM10; **owner:** R06; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V053"></a>

### V053: DEqMS statistic alignment

**Fixture:** Same count fixture with shuffled count-table rows and intentionally different native limma versus sca P fields.

**Oracle:** Join by feature_id and direct native field names/values.

**Exact assertion:** Align counts by identity and retain sca.P.Value/statistic/variance distinctly from original limma fields; compute central q from sca raw P.

**Negative case:** Duplicate/missing count IDs fail; silently using original P.Value as DEqMS P fails exact field assertions.

**Contract:** SM10; **owner:** R06; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V054"></a>

### V054: DEqMS hypotheses and design guard

**Fixture:** One eligible sparse independent fixture plus blocked, weighted and treat requests.

**Oracle:** SM10 allowed-design table plus direct exact-refit DEqMS reference for the eligible sparse case.

**Exact assertion:** Qualify the sparse exact route numerically and reject unqualified block/precision-weight/TREAT combinations before backend execution.

**Negative case:** A backend accepting and ignoring a block/weight argument cannot count as supported behavior.

**Contract:** SM10; **owner:** R06; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V055"></a>

### V055: LFQ dropout provenance gate

**Fixture:** Documented unimputed LFQ matrix with real NA mask; copies relabeled TMT, unknown-imputed and duplicate-correlation.

**Oracle:** Declared assay/mask/imputation/design eligibility, not result count.

**Exact assertion:** Permit only the qualified LFQ independent fixed-model path; preserve actual original observations and native missingness.

**Negative case:** TMT, unknown incompatible prior imputation or repeated block models produce the named proDA eligibility errors, never a limma fallback.

**Contract:** SM10; **owner:** R06; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V056"></a>

### V056: proDA fitting and test_diff

**Fixture:** ≤20 LFQ features with dropout, including a feature absent in one target group; full and nested reduced fixed designs.

**Oracle:** Separate pinned proDA/test_diff calls for the zero-null and reduced-model tests.

**Exact assertion:** Match native effects, P and native uncertainty/diagnostics where reference results exist; retain n_obs and prior/dropout dependence for all-missing-group estimates.

**Negative case:** Invented zero abundance, substituted residual df or a fake TREAT endpoint fails even if a result row is finite.

**Contract:** SM10; **owner:** R06; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V057"></a>

### V057: Alternative-engine uncertainty

**Fixture:** One DEqMS and one proDA direct-reference fit with their native SE/variance/df fields.

**Oracle:** Model-specific direct calculations permitted by the pinned official output semantics.

**Exact assertion:** Export only justified intervals/statistics with correct native names/df; unsupported quantities stay unavailable with reasons.

**Negative case:** Borrowing limma df or SE to fill missing proDA/DEqMS uncertainty is rejected.

**Contract:** SM10; **owner:** R06; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V058"></a>

### V058: Method sensitivity comparison

**Fixture:** Two predeclared engines with overlapping but nonidentical eligible feature universes and different rejection counts.

**Oracle:** Stable-ID intersection/union and the unchanged frozen primary declaration.

**Exact assertion:** Compare matched effects and disclose unmatched universes; keep model/hypothesis families distinct and never promote the larger hit count to primary.

**Negative case:** Shuffle engine result order or increase sensitivity hits: primary engine/model/plan hash must not change.

**Contract:** SM10; **owner:** R06; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V059"></a>

### V059: Dependencies and failure propagation

**Fixture:** A schema-valid Phase 2 config with installed eligible primary limma model declaring execution_requirement=required and a DEqMS sensitivity model declaring execution_requirement=optional whose package is absent; then the same eligible DEqMS model declaring required; then an installed child process that raises an error.

**Oracle:** Frozen model declarations plus actual process/dependency state and the CLI exit/state table, evaluated without deleting unavailable requests from the plan.

**Exact assertion:** Structural/scientific validation accepts the contractually eligible optional request independent of installation and freezes its requiredness. Its missing package yields NOT_RUN plus E_ENGINE_NOT_AVAILABLE and makes the otherwise successful run PARTIAL/exit 3. The same missing eligible engine declared required makes the run FAILED/exit 3 with dependants NOT_RUN. An actual optional child failure yields FAILED at that stage and PARTIAL/exit 4; an actual required child failure makes the run FAILED/exit 4.

**Negative case:** A primary model declared optional fails structural validation; omission of execution_requirement fails schema validation. A validator that rejects an eligible optional request solely because its package is absent, drops it from the frozen plan, emits a fallback table under its name or labels the operational absence INAPPLICABLE fails.

**Contract:** SM10; **owner:** R06; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V060"></a>

### V060: Assay-engine golden/calibration fixtures

**Fixture:** Small count-dependent-variance and LFQ-dropout fixtures with a fixed seed, plus unsuitable-input negatives.

**Oracle:** Independent direct package references; release-scale calibration is explicitly V105–V106, not this tiny fixture.

**Exact assertion:** Prespecify reference calls/expected statuses and save actual diagnostics, precision and eligibility outcomes; no tiny-fixture claim of calibrated FDR.

**Negative case:** Nonconvergence or missing R is retained as failure/NOT_RUN, not replaced with fixture-derived expected statistics.

**Contract:** SM10; **owner:** R06; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

## Boundary

No recovered source/audit changes, private-data transfer, invented license, fake backend or unsupported completion claim. A missing prerequisite is NOT_RUN, not a successful acceptance case. Only the Maintainer updates traceability after independent gate execution.
