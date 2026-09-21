# Feature Specification: Canonical protein input and biological identity

**Phase:** 1. **Packet:** R02. **Status:** 1.2.0-frozen; implementation pending prerequisite acceptance and explicit authorization.

## Scope

Deliver only FR-011–FR-020 for US1. [Packet index](../001-downstream-proteomics/packet-index.md), [ownership](../001-downstream-proteomics/packet-ownership.json), [data model](../001-downstream-proteomics/data-model.md) and [scientific contract](../001-downstream-proteomics/contracts/scientific-methods.md) define the shared interfaces. No implementation dispatch is authorized yet.

## Requirements

- **FR-011 — Canonical wide and long intake:** The system MUST produce exactly the same ordered values and masks from both formats, with no duplicate or missing key silently accepted.
- **FR-012 — Explicit scale and missing encoding:** The system MUST transform positive linear values once; preserve log2 zero/negative values; distinguish explicitly encoded linear missing zero from invalid nonpositive observations.
- **FR-013 — Biological observation hierarchy:** The system MUST report injections=4 and biological n=1 for the first specimen; paired visits retain four subjects and do not become eight independent subjects.
- **FR-014 — Explicit technical aggregation:** The system MUST match the two different means within transform tolerance and retain one biological unit with source coverage/lineage.
- **FR-015 — Feature annotations and protein groups:** The system MUST keep three unique feature rows, complete constituent members and unknown flags; do not use gene symbols as unique protein IDs.
- **FR-016 — Versioned vendor mappings:** The system MUST validate the named profile/version and map only its declared protein-level columns without inferring a vendor version from similar names.
- **FR-017 — Historical workbook importer:** The system MUST wrap the existing parser without source edits; retain sheet identities and require explicit noninterchangeable Method A/B naming.
- **FR-018 — Immutable canonical bundle:** The system MUST publish a self-consistent canonical manifest; roundtrip values/masks exactly; refuse collisions and changed source bytes.
- **FR-019 — Intake validation report:** The system MUST export grain, scale, biological/technical n, missingness, errors and provenance gaps with no asserted confirmation of unknown biological facts.
- **FR-020 — Input edge-case fixtures:** The system MUST reject each invalid file with a specific field/row reason; correctly roundtrip a valid quoted delimiter and Unicode identifier.

## Acceptance Scenarios

<a id="V011"></a>

### V011: Canonical wide and long intake

**Fixture:** 8×12 independent fixture and equivalent shuffled long form, including its two NA cells.

**Oracle:** Pivot/alignment from explicit feature-observation keys.

**Exact assertion:** Produce exactly the same ordered values and masks from both formats, with no duplicate or missing key silently accepted.

**Negative case:** Add a duplicate feature-observation pair or duplicate header: E_ID_DUPLICATE and no canonical publication.

**Contract:** SM01; **owner:** R02; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V012"></a>

### V012: Explicit scale and missing encoding

**Fixture:** One-row values [1,2,4,0,-1] under separate linear/log2 configurations and explicit zero-encoding variants.

**Oracle:** log2([1,2,4])=[0,1,2] and the SM02 transition table.

**Exact assertion:** Transform positive linear values once; preserve log2 zero/negative values; distinguish explicitly encoded linear missing zero from invalid nonpositive observations.

**Negative case:** source_scale=log2 with transform=log2 fails E_SCALE_SECOND_LOG before a transform is executed.

**Contract:** SM02; **owner:** R02; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V013"></a>

### V013: Biological observation hierarchy

**Fixture:** Four injections of one specimen and a four-subject paired fixture with distinct visit specimens.

**Oracle:** Unique biological-unit counts within group and unique subject IDs.

**Exact assertion:** Report injections=4 and biological n=1 for the first specimen; paired visits retain four subjects and do not become eight independent subjects.

**Negative case:** Unaggregated injections with technical mode none fail E_TECHNICAL_REPLICATION_UNMODELED.

**Contract:** SM03; **owner:** R02; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V014"></a>

### V014: Explicit technical aggregation

**Fixture:** Two injections with linear abundances [4,16], plus a distinct visit that must remain separate.

**Oracle:** mean_linear→log2(10); mean_log2→3.

**Exact assertion:** Match the two different means within transform tolerance and retain one biological unit with source coverage/lineage.

**Negative case:** An aggregation group spanning before and after or different specimens fails instead of averaging the biological effect away.

**Contract:** SM03; **owner:** R02; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V015"></a>

### V015: Feature annotations and protein groups

**Fixture:** Three opaque features sharing an accession/symbol, one multi-gene group and unknown decoy flags.

**Oracle:** Explicit input members and tri-state flag values.

**Exact assertion:** Keep three unique feature rows, complete constituent members and unknown flags; do not use gene symbols as unique protein IDs.

**Negative case:** A duplicate feature_id fails even when its annotations differ; unknown decoy is not coerced to false.

**Contract:** SM01; **owner:** R02; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V016"></a>

### V016: Versioned vendor mappings

**Fixture:** Separate ≤8-feature synthetic protein-level exports with documented header sets for DIA-NN, MaxQuant, FragPipe and Spectronaut.

**Oracle:** The frozen per-profile column maps and canonical expected rows, authored before adapter output.

**Exact assertion:** Validate the named profile/version and map only its declared protein-level columns without inferring a vendor version from similar names.

**Negative case:** A peptide/precursor/phosphosite-grain table fails E_UNSUPPORTED_SCOPE; unknown vendor headers need an explicit mapping, not a guessed one.

**Contract:** SM01; **owner:** R02; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V017"></a>

### V017: Historical workbook importer

**Fixture:** A tiny synthetic four-sheet workbook and explicit Method A/Method B contrast/group crosswalk.

**Oracle:** Read-only recovered parser behavior and crosswalk by biological numerator/denominator, not D1/T1 ordering.

**Exact assertion:** Wrap the existing parser without source edits; retain sheet identities and require explicit noninterchangeable Method A/B naming.

**Negative case:** Swap D1/D2 or T1/T2 without updating the semantic crosswalk: fail contrast identity validation rather than relabeling results.

**Contract:** SM01; **owner:** R02; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V018"></a>

### V018: Immutable canonical bundle

**Fixture:** Canonical matrix/masks/metadata from V011 and an attempted same-output rerun.

**Oracle:** Roundtrip IDs/values plus independent SHA-256.

**Exact assertion:** Publish a self-consistent canonical manifest; roundtrip values/masks exactly; refuse collisions and changed source bytes.

**Negative case:** Delete a mask or change one observation ID: verify fails integrity and no downstream fit starts.

**Contract:** SM01; **owner:** R02; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V019"></a>

### V019: Intake validation report

**Fixture:** Small input with unknown tissue/upstream normalization and declared technical units.

**Oracle:** Input declarations, independent dimensions and missing-cell counts.

**Exact assertion:** Export grain, scale, biological/technical n, missingness, errors and provenance gaps with no asserted confirmation of unknown biological facts.

**Negative case:** Missing tissue/assay evidence remains unknown or blocks its dependent analysis; it does not become a hardcoded healthy/PASS statement.

**Contract:** SM01; **owner:** R02; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V020"></a>

### V020: Input edge-case fixtures

**Fixture:** Individually malformed tiny files: ragged row, duplicated ID, mismatched mask, Inf, quoted TAB, unknown encoding, unsupported assay.

**Oracle:** Parser error location and correct CSV/TSV escaping roundtrip.

**Exact assertion:** Reject each invalid file with a specific field/row reason; correctly roundtrip a valid quoted delimiter and Unicode identifier.

**Negative case:** Reject Inf/nonfinite observed numeric data, rather than silently treating it as missing or clipping it.

**Contract:** SM01; **owner:** R02; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

## Boundary

No recovered source/audit changes, private-data transfer, invented license, fake backend or unsupported completion claim. A missing prerequisite is NOT_RUN, not a successful acceptance case. Only the Maintainer updates traceability after independent gate execution.
