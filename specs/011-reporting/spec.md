# Feature Specification: Unified offline report and streamlined workflow

**Phase:** 1 (R10a: 091–094); 3 (R10b: 095–100). **Packet:** R10. **Status:** 1.2.0-frozen; implementation pending exact prerequisites and authorization.

## Scope

Deliver only FR-091–FR-100 for US5. [Packet index](../001-downstream-proteomics/packet-index.md), [ownership](../001-downstream-proteomics/packet-ownership.json), [data model](../001-downstream-proteomics/data-model.md) and [scientific contract](../001-downstream-proteomics/contracts/scientific-methods.md) define the shared interfaces. No implementation dispatch is authorized yet.

## Requirements

- **FR-091 — Typed report data assembler:** The system MUST assemble only verified executed values; zero real discoveries remain zero, missing values remain unknown and failures remain failed; core methods/status summary derives actual metadata.
- **FR-092 — One-run command integration:** The system MUST execute the valid Phase 1 DAG and R10a report from one command; stop required dependants on failure, retain a useful partial report and correct nonzero exit.
- **FR-093 — Offline HTML report:** The system MUST render readable offline HTML with clear scope/status, navigation and accessible tables, and redact private absolute input paths in presentation.
- **FR-094 — QC and inclusion/exclusion sections:** The system MUST show biological versus technical n, scale, missingness, normalization and exclusion rationale with exact source links; no sample silently disappears.
- **FR-095 — Complete differential results views:** The system MUST display every planned contrast/hypothesis, effects/CI/q-family and estimability, including auxiliary comparisons without highlighting only significant subsets.
- **FR-096 — Pathway and response views:** The system MUST show null types, mapped universe, leading edges, exploratory labels and response restrictions with no cross-null confirmation language.
- **FR-097 — Publication figure and source exports:** The system MUST export requested PDF/SVG/PNG with readable units, independent n, thresholds and complete source data; no unexplained significance stars.
- **FR-098 — Methods and limitations generation:** The system MUST generate methods naming actual inputs/scale/models/hypotheses/families/resources and uncertainty limits, keeping unknown provenance explicit.
- **FR-099 — Null/empty/partial report robustness:** The system MUST render truthful full reports for all cases without upgrading partial/failed/cancelled analysis to COMPLETED; preserve valid zero-row tables where a stage actually completed.
- **FR-100 — Quickstart and user workflow:** The system MUST exercise validate, plan, run, inspect/report, verify and compare without manual R stage orchestration and document only implemented capabilities.

## Acceptance Scenarios

<a id="V091"></a>

### V091: Typed report data assembler

**Fixture:** Real stage manifests for completed QC/DEA with zero rejections, a missing DEA stage and a failed required stage.

**Oracle:** Manifest states, actual table row/count values and strict ReportData schema.

**Exact assertion:** Assemble only verified executed values; zero real discoveries remain zero, missing values remain unknown and failures remain failed; core methods/status summary derives actual metadata.

**Negative case:** Hardcoded PASS, absent DEA interpreted as healthy zero, or an unverified artifact included as completed fails.

**Contract:** SM25; **owner:** R10a; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V092"></a>

### V092: One-run command integration

**Fixture:** Phase 1 synthetic configuration through the actual Python→R stages, repeated with a genuine R failure.

**Oracle:** Observed stage order, plan timestamp/hash before first fit, actual child exits and resulting manifest.

**Exact assertion:** Execute the valid Phase 1 DAG and R10a report from one command; stop required dependants on failure, retain a useful partial report and correct nonzero exit.

**Negative case:** A Phase 1 score request fails eligibility rather than running a hidden score backend or emitting score P values.

**Contract:** SM25; **owner:** R10a; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V093"></a>

### V093: Offline HTML report

**Fixture:** Thin R10a report with embedded CSS, plain tables, local source links and no JavaScript/network dependency.

**Oracle:** Open with networking disabled and inspect DOM/resource references.

**Exact assertion:** Render readable offline HTML with clear scope/status, navigation and accessible tables, and redact private absolute input paths in presentation.

**Negative case:** Any external CSS/font/script request or a broken mandatory local table link fails.

**Contract:** SM25; **owner:** R10a; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V094"></a>

### V094: QC and inclusion/exclusion sections

**Fixture:** QC with four injections aggregated to one specimen, a justified exclusion and missing upstream metadata.

**Oracle:** Actual ObservationHierarchy, ExclusionDecision and QC source TSV values.

**Exact assertion:** Show biological versus technical n, scale, missingness, normalization and exclusion rationale with exact source links; no sample silently disappears.

**Negative case:** A display count based on injections or a missing exclusion reason rendered as confirmed removal fails.

**Contract:** SM06; **owner:** R10a; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V095"></a>

### V095: Complete differential results views

**Fixture:** Complete all-contrast tables including secondary treated-control, nonsignificant and nonestimable rows.

**Oracle:** Keyed comparison of rendered ReportData and full DEA tables.

**Exact assertion:** Display every planned contrast/hypothesis, effects/CI/q-family and estimability, including auxiliary comparisons without highlighting only significant subsets.

**Negative case:** Dropping an auxiliary contrast or substituting display-filtered counts for total tested family counts fails.

**Contract:** SM12; **owner:** R10b; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V096"></a>

### V096: Pathway and response views

**Fixture:** Mixed CAMERA/ROAST/fgsea/ORA artifacts plus descriptive response, equivalence and unavailable score inference.

**Oracle:** Each typed result/null/resource/universe and eligibility artifact.

**Exact assertion:** Show null types, mapped universe, leading edges, exploratory labels and response restrictions with no cross-null confirmation language.

**Negative case:** fgsea presented as sample-level replication or descriptive-score P columns in the report fails.

**Contract:** SM16; **owner:** R10b; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V097"></a>

### V097: Publication figure and source exports

**Fixture:** A small contrast plot and pathway/response plot with finite/NA and no discoveries.

**Oracle:** Exact plotted coordinates/labels versus saved figure-source tables plus actual file format headers.

**Exact assertion:** Export requested PDF/SVG/PNG with readable units, independent n, thresholds and complete source data; no unexplained significance stars.

**Negative case:** A figure with altered rounded data used for inference, fabricated points, unreadable labels or no source table fails visual/numeric review.

**Contract:** SM25; **owner:** R10b; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V098"></a>

### V098: Methods and limitations generation

**Fixture:** Two actual runs differing in hypothesis/engine and one with unknown tissue/batch provenance.

**Oracle:** Resolved configs, executed package inventory and method eligibility/status, not study-specific boilerplate.

**Exact assertion:** Generate methods naming actual inputs/scale/models/hypotheses/families/resources and uncertainty limits, keeping unknown provenance explicit.

**Negative case:** Mentioning an unexecuted engine, unverified tissue or canned biological rescue claim fails.

**Contract:** SM25; **owner:** R10b; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V099"></a>

### V099: Null/empty/partial report robustness

**Fixture:** No discoveries, empty pathway foreground, tiny QC-only cohort, optional method missing, required backend crash and interrupted stage.

**Oracle:** Stage/RunStatus and expected numerical section availability for each scenario.

**Exact assertion:** Render truthful full reports for all cases without upgrading partial/failed/cancelled analysis to COMPLETED; preserve valid zero-row tables where a stage actually completed.

**Negative case:** An exception-handling path that silently substitutes zeros or hardcoded healthy prose fails.

**Contract:** SM25; **owner:** R10b; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V100"></a>

### V100: Quickstart and user workflow

**Fixture:** Fresh documented setup and all three public synthetic examples; one run comparison through real later commands.

**Oracle:** Copy documented commands verbatim in a clean environment and inspect actual outputs.

**Exact assertion:** Exercise validate, plan, run, inspect/report, verify and compare without manual R stage orchestration and document only implemented capabilities.

**Negative case:** A command present only in prose, an undisclosed missing prerequisite or a README claiming unshipped phases fails.

**Contract:** SM25; **owner:** R10b; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

## Boundary

No recovered source/audit changes, private-data transfer, invented license, fake backend or unsupported completion claim. A missing prerequisite is NOT_RUN, not a successful acceptance case. Only the Maintainer updates traceability after independent gate execution.
