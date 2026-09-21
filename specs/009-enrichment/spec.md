# Feature Specification: Design-compatible pathways and enrichment

**Phase:** 2. **Packet:** R08. **Status:** 1.2.0-frozen; implementation pending whole-group acceptance and separate authorization.

## Dispatch before enrichment

| Design | Required primary matrix method | Rejected / exploratory |
|---|---|---|
| Independent finite linear | CAMERA | No block ignored. |
| Paired fixed-subject | ROAST with fixed-subject design | CAMERA rejected. |
| Repeated duplicateCorrelation | ROAST with validated block/correlation | CAMERA rejected; no double subject encoding. |
| Native ranks only / alternative primary engine | No automatically equivalent matrix test | fgsea exploratory; separate linear sensitivity explicitly named. |

This is a navigation summary; the normative [SM16 dispatch table](../001-downstream-proteomics/contracts/scientific-methods.md#design-method-dispatch) controls eligibility/nulls.

## Scope

Deliver only FR-071–FR-080 for US3. [Packet index](../001-downstream-proteomics/packet-index.md), [ownership](../001-downstream-proteomics/packet-ownership.json), [data model](../001-downstream-proteomics/data-model.md) and [scientific contract](../001-downstream-proteomics/contracts/scientific-methods.md) define the shared interfaces. No implementation dispatch is authorized yet.

## Requirements

- **FR-071 — Pathway null and method dispatcher:** The system MUST select CAMERA only for eligible independent matrix designs and ROAST for paired/blocked designs; label fgsea exploratory gene-set null and ORA separately.
- **FR-072 — Independent-design CAMERA:** The system MUST match native statistic/P/direction within tolerance and retain correlation/universe diagnostics plus separate central q.
- **FR-073 — Blocked-design ROAST:** The system MUST match reference endpoints with the same RNG; record correct covariance route, nrot and smoke label.
- **FR-074 — Directional and mixed rotation endpoints:** The system MUST export separate self_contained_directional and self_contained_mixed rows/families and recompute each declared central q from its correct raw P.
- **FR-075 — Exploratory fgsea adapter:** The system MUST retain ES/NES/P/native q/central q/log2err/leading-edge/size and warnings with preranked_gene_set hypothesis and exploratory status.
- **FR-076 — Exact background-aware ORA:** The system MUST match exact raw and adjusted values and retain the k=0 set until after BH; empty foreground gives p=1 for every eligible set.
- **FR-077 — Pathway family adjustment:** The system MUST keep method/null/role and declared collection/contrast membership exact; preserve native adjusted values separately from central q.
- **FR-078 — Alternative-engine pathway sensitivity:** The system MUST label fgsea exploratory and linear CAMERA/ROAST results as a separate sensitivity, never same-likelihood primary-engine confirmation.
- **FR-079 — Leading-edge redundancy and diagnostics:** The system MUST display overlap relations while preserving every original hypothesis ID, P/q and universe; no causal mechanism or independent-discovery count is created from clustering.
- **FR-080 — Enrichment golden and failure suite:** The system MUST produce honest complete empty/no-discovery results where appropriate, and explicit failures for missing prerequisites or all eligible numerical tests failing.

## Acceptance Scenarios

<a id="V071"></a>

### V071: Pathway null and method dispatcher

**Fixture:** Independent, paired fixed-subject, duplicate-correlation and native-rank-only plans.

**Oracle:** The SM16 design × method dispatch table.

**Exact assertion:** Select CAMERA only for eligible independent matrix designs and ROAST for paired/blocked designs; label fgsea exploratory gene-set null and ORA separately.

**Negative case:** Any paired/fixed-subject CAMERA exception, ignored block or sample-level fgsea claim fails.

**Contract:** SM16; **owner:** R08; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V072"></a>

### V072: Independent-design CAMERA

**Fixture:** Finite independent gene matrix with two overlapping sets and an identifiable contrast.

**Oracle:** Direct pinned limma::camera using inter.gene.cor=NA and identical gene matrix/design/contrast.

**Exact assertion:** Match native statistic/P/direction within tolerance and retain correlation/universe diagnostics plus separate central q.

**Negative case:** Supplying a subject block causes E_CAMERA_BLOCKED_DESIGN instead of being passed and ignored.

**Contract:** SM16; **owner:** R08; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V073"></a>

### V073: Blocked-design ROAST

**Fixture:** Four-subject paired finite matrix and a separate eligible duplicate-correlation model; fixed seed and smoke nrot=199.

**Oracle:** Direct pinned mroast/roast calls: fixed-subject design alone versus verified block/correlation, midp=FALSE.

**Exact assertion:** Match reference endpoints with the same RNG; record correct covariance route, nrot and smoke label.

**Negative case:** Double-encoded subjects, ignored block or production nrot=199 fails eligibility/resolution checks.

**Contract:** SM16; **owner:** R08; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V074"></a>

### V074: Directional and mixed rotation endpoints

**Fixture:** A set containing opposing effects and another with concordant effects in the same rotation reference run.

**Oracle:** Native directional versus mixed raw columns, retained before central adjustment.

**Exact assertion:** Export separate self_contained_directional and self_contained_mixed rows/families and recompute each declared central q from its correct raw P.

**Negative case:** Swapping mixed/directional P or labeling native per-call FDR as global central q fails.

**Contract:** SM16; **owner:** R08; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V075"></a>

### V075: Exploratory fgsea adapter

**Fixture:** A finite unique-gene rank list with ties, tiny fixed gene sets, explicit eps and fixed RNG/thread settings.

**Oracle:** Separate pinned fgsea call with stable gene-ID tie order and identical settings.

**Exact assertion:** Retain ES/NES/P/native q/central q/log2err/leading-edge/size and warnings with preranked_gene_set hypothesis and exploratory status.

**Negative case:** Nonfinite/duplicate ranks or a TREAT P transformed into a purported native zero-null rank fail E_FGSEA_RANK_INVALID.

**Contract:** SM16; **owner:** R08; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V076"></a>

### V076: Exact background-aware ORA

**Fixture:** Universe N=10, foreground n=2, eligible sets K=2 with k=2 and K=2 with k=0; bounds explicitly 1–10.

**Oracle:** Hypergeometric tail: p1=1/45, p2=1; BH across both gives q1=2/45, q2=1.

**Exact assertion:** Match exact raw and adjusted values and retain the k=0 set until after BH; empty foreground gives p=1 for every eligible set.

**Negative case:** Dropping the zero-hit set before adjustment changes q1 to 1/45 and must fail this gate.

**Contract:** SM17; **owner:** R08; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V077"></a>

### V077: Pathway family adjustment

**Fixture:** Two contrasts, two collections and multiple null types with small explicit raw-P arrays.

**Oracle:** Prespecified family membership expansion and independent BH/BY on each selected array.

**Exact assertion:** Keep method/null/role and declared collection/contrast membership exact; preserve native adjusted values separately from central q.

**Negative case:** Pooling ROAST mixed with CAMERA competitive or omitting a tested zero-hit ORA set fails.

**Contract:** SM12; **owner:** R08; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V078"></a>

### V078: Alternative-engine pathway sensitivity

**Fixture:** A proDA primary model with native zero-null ranks and an explicitly named complete-case limma gene-model sensitivity.

**Oracle:** Actual model IDs, likelihood/null labels and corresponding direct method calls.

**Exact assertion:** Label fgsea exploratory and linear CAMERA/ROAST results as a separate sensitivity, never same-likelihood primary-engine confirmation.

**Negative case:** A CAMERA row carrying engine=proda or the primary proDA model identity fails E_PATHWAY_ENGINE_MISMATCH.

**Contract:** SM16; **owner:** R08; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V079"></a>

### V079: Leading-edge redundancy and diagnostics

**Fixture:** Three sets with shared leading-edge genes and distinct set IDs/P values.

**Oracle:** Hand Jaccard/intersection counts from the exported memberships.

**Exact assertion:** Display overlap relations while preserving every original hypothesis ID, P/q and universe; no causal mechanism or independent-discovery count is created from clustering.

**Negative case:** Merging sets and summing them as independent confirmations or changing their original P values fails.

**Contract:** SM18; **owner:** R08; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V080"></a>

### V080: Enrichment golden and failure suite

**Fixture:** V072–V076 references plus missing snapshot, NA/Inf matrix, no foreground and all-nonfinite eligible test outputs.

**Oracle:** Direct references, schema/state contract and hand empty-set behavior.

**Exact assertion:** Produce honest complete empty/no-discovery results where appropriate, and explicit failures for missing prerequisites or all eligible numerical tests failing.

**Negative case:** All-nonfinite tests cannot become a successful zero-discovery family or a blank healthy report section.

**Contract:** SM18; **owner:** R08; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

## Boundary

No recovered source/audit changes, private-data transfer, invented license, fake backend or unsupported completion claim. A missing prerequisite is NOT_RUN, not a successful acceptance case. Only the Maintainer updates traceability after independent gate execution.
