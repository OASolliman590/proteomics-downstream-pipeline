# Feature Specification: Versioned annotation, protein groups and gene sets

**Phase:** 2. **Packet:** R07. **Status:** 1.2.0-frozen; implementation pending Phase 1 acceptance and separate authorization.

## Scope

Deliver only FR-061–FR-070 for US3. [Packet index](../001-downstream-proteomics/packet-index.md), [ownership](../001-downstream-proteomics/packet-ownership.json), [data model](../001-downstream-proteomics/data-model.md) and [scientific contract](../001-downstream-proteomics/contracts/scientific-methods.md) define the shared interfaces. No implementation dispatch is authorized yet.

## Requirements

- **FR-061 — Local resource snapshot registry:** The system MUST verify every required snapshot before analysis; preserve terms and original source/target organism identity.
- **FR-062 — Explicit resource preparation command:** The system MUST acquire/build only during the requested preparation action; record source version/hash/terms, then consume a local verified snapshot offline.
- **FR-063 — Species-aware identifier mapping:** The system MUST map only the declared taxonomy/ID namespace and retain retired/unmapped records with reasons.
- **FR-064 — Orthology evidence preservation:** The system MUST retain projection evidence/version and label the collection ortholog-projected, not rat-native.
- **FR-065 — Label-independent gene representatives:** The system MUST select the exact winner at each tie stage and keep the winner unchanged under every group-label permutation.
- **FR-066 — Ambiguous protein-group policy:** The system MUST exclude the multi-gene group from the default gene matrix with a reason; retain an unambiguous same-gene group and separate any median-aggregation sensitivity universe.
- **FR-067 — Finite pathway matrix and model:** The system MUST use the finite measured subset without imputation, quantify lost genes/set coverage and fit the corresponding gene model consistently.
- **FR-068 — Gene-set overlap eligibility:** The system MUST retain full/mapped/eligible membership counts and select only the size-3 set before testing; test-only bounds remain labeled.
- **FR-069 — Measured gene background for ORA:** The system MUST use those six genes as ORA universe and document predictable universe changes under declared filtering.
- **FR-070 — Resource/mapping regression fixtures:** The system MUST reproduce mapping/representative/universe identities offline, retaining each loss and corrupted-resource failure.

## Acceptance Scenarios

<a id="V061"></a>

### V061: Local resource snapshot registry

**Fixture:** A tiny local mapping/gene-set snapshot with actual SHA-256, release/taxonomy/source/terms metadata.

**Oracle:** Independent hash of actual local bytes and explicit manifest fields.

**Exact assertion:** Verify every required snapshot before analysis; preserve terms and original source/target organism identity.

**Negative case:** One-byte mutation or missing file fails E_RESOURCE_HASH/E_RESOURCE_MISSING; an all-zero placeholder digest cannot pass.

**Contract:** SM14; **owner:** R07; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V062"></a>

### V062: Explicit resource preparation command

**Fixture:** Small permissibly redistributable source and an explicit prepare invocation, followed by a network-disabled production run.

**Oracle:** Preparation manifest and observable absence of analysis-time network calls.

**Exact assertion:** Acquire/build only during the requested preparation action; record source version/hash/terms, then consume a local verified snapshot offline.

**Negative case:** Running analysis with no cached snapshot fails rather than invoking msigdbr or downloading latest content.

**Contract:** SM14; **owner:** R07; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V063"></a>

### V063: Species-aware identifier mapping

**Fixture:** Identical gene symbols in two taxonomy namespaces, one retired ID and one unmapped feature.

**Oracle:** Hand-authored species-keyed mapping relations.

**Exact assertion:** Map only the declared taxonomy/ID namespace and retain retired/unmapped records with reasons.

**Negative case:** A cross-species symbol-only join fails E_RESOURCE_TAXONOMY instead of merging unrelated genes.

**Contract:** SM15; **owner:** R07; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V064"></a>

### V064: Orthology evidence preservation

**Fixture:** Tiny human→rat projection with one ambiguous ortholog and one unambiguous counterpart.

**Oracle:** Explicit source/target IDs and evidence from the fixture manifest.

**Exact assertion:** Retain projection evidence/version and label the collection ortholog-projected, not rat-native.

**Negative case:** Missing source taxonomy or an unsupported rat-native label fails provenance/report validation.

**Contract:** SM14; **owner:** R07; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V065"></a>

### V065: Label-independent gene representatives

**Fixture:** Three proteins per one gene: differing coverage; then equal coverage/differing median; then both tied/differing stable IDs.

**Oracle:** Lexicographic key (descending overall coverage, descending overall median, ascending feature_id), calculated without labels.

**Exact assertion:** Select the exact winner at each tie stage and keep the winner unchanged under every group-label permutation.

**Negative case:** A max-abs(t), minimum-P or disease-group-specific coverage winner fails even if it yields more significant pathways.

**Contract:** SM15; **owner:** R07; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V066"></a>

### V066: Ambiguous protein-group policy

**Fixture:** One unique single-gene feature, one multi-accession same-gene group and one group mapping to two genes.

**Oracle:** Complete mapping cardinality and explicitly declared aggregation sensitivity.

**Exact assertion:** Exclude the multi-gene group from the default gene matrix with a reason; retain an unambiguous same-gene group and separate any median-aggregation sensitivity universe.

**Negative case:** Duplicating a multi-gene protein into both genes in the default matrix fails.

**Contract:** SM15; **owner:** R07; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V067"></a>

### V067: Finite pathway matrix and model

**Fixture:** Gene matrix with one required-observation NA and a distinct complete-case subset.

**Oracle:** Exact complete-row mask and a direct refit on that same finite gene matrix.

**Exact assertion:** Use the finite measured subset without imputation, quantify lost genes/set coverage and fit the corresponding gene model consistently.

**Negative case:** Using maximum protein ranks from discarded representatives or a silently imputed CAMERA/ROAST matrix fails.

**Contract:** SM15; **owner:** R07; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V068"></a>

### V068: Gene-set overlap eligibility

**Fixture:** A 12-gene universe and sets with intersections of sizes 0,1,3 and 6, using declared test bounds 2–5.

**Oracle:** Hand intersections of each set with the eligible universe.

**Exact assertion:** Retain full/mapped/eligible membership counts and select only the size-3 set before testing; test-only bounds remain labeled.

**Negative case:** Selecting eligible sets after inspecting their foreground hits or P values fails.

**Contract:** SM16; **owner:** R07; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V069"></a>

### V069: Measured gene background for ORA

**Fixture:** Eight measured mapped genes, of which six pass model eligibility, and two unmeasured genes in a resource.

**Oracle:** Exact set of six eligible measured mapped IDs.

**Exact assertion:** Use those six genes as ORA universe and document predictable universe changes under declared filtering.

**Negative case:** Adding the resource genome or unmeasured genes to increase apparent enrichment fails.

**Contract:** SM17; **owner:** R07; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V070"></a>

### V070: Resource/mapping regression fixtures

**Fixture:** One-to-many, duplicate-symbol, projected and zero-overlap cases plus an intentionally corrupted local snapshot.

**Oracle:** Hand mappings, overlap counts and direct SHA-256.

**Exact assertion:** Reproduce mapping/representative/universe identities offline, retaining each loss and corrupted-resource failure.

**Negative case:** Reordering labels or repeated execution cannot alter representative membership; a corrupt snapshot cannot be silently refreshed.

**Contract:** SM14; **owner:** R07; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

## Boundary

No recovered source/audit changes, private-data transfer, invented license, fake backend or unsupported completion claim. A missing prerequisite is NOT_RUN, not a successful acceptance case. Only the Maintainer updates traceability after independent gate execution.
