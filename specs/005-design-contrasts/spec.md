# Feature Specification: Design validation, exact contrasts and blocking

**Phase:** 1. **Packet:** R04. **Status:** 1.2.0-frozen; implementation pending prerequisite acceptance and explicit authorization.

## Scope

Deliver only FR-031–FR-040 for US2. [Packet index](../001-downstream-proteomics/packet-index.md), [ownership](../001-downstream-proteomics/packet-ownership.json), [data model](../001-downstream-proteomics/data-model.md) and [scientific contract](../001-downstream-proteomics/contracts/scientific-methods.md) define the shared interfaces. No implementation dispatch is authorized yet.

## Requirements

- **FR-031 — Safe declarative design grammar:** The system MUST match the frozen coefficient map and numeric design, using the documented escaping/interaction convention and no arbitrary R evaluation.
- **FR-032 — Stable levels and aligned matrices:** The system MUST preserve aligned coefficients and reference meanings across row/column reorderings.
- **FR-033 — Rank and confounding diagnostics:** The system MUST reject the confounded model before fit with E_DESIGN_CONFOUNDED and aliased terms; accept the full-rank crossed design.
- **FR-034 — Independent, paired and repeated plans:** The system MUST preserve subject identity, pass the correct fixed or random representation and record actual consensus correlation for the repeated fit.
- **FR-035 — Numeric contrasts and interaction semantics:** The system MUST match d/t/r and r=d+t; calculate the interaction directly from coefficient weights, not from P-value comparisons.
- **FR-036 — Featurewise estimability and df:** The system MUST export correct featurewise n/df/estimability and named reasons for ordinary available-case failure; retain excluded rows.
- **FR-037 — Exact contrast covariance strategy:** The system MUST match exact SE/covariance within backend tolerance and record the exact path; demonstrate the approximate shortcut differs on at least one chosen fixture.
- **FR-038 — Method-design applicability registry:** The system MUST return each documented eligibility/error type before execution and retain one authoritative primary engine.
- **FR-039 — Frozen analysis plan and matrices:** The system MUST write plan, masks, evaluated design/contrasts, hypotheses and families before any fit; each meaningful change produces a different plan hash.
- **FR-040 — Design adversarial verification:** The system MUST accept only identifiable sufficiently replicated contrasts; exact key alignment makes reordering invariant; singleton inference fails while descriptive QC remains available.

## Acceptance Scenarios

<a id="V031"></a>

### V031: Safe declarative design grammar

**Fixture:** Two groups, a numeric age covariate, one batch factor and group×age interaction; malicious function-call text as a term.

**Oracle:** Explicit design rows from indicators/centered covariates and product columns.

**Exact assertion:** Match the frozen coefficient map and numeric design, using the documented escaping/interaction convention and no arbitrary R evaluation.

**Negative case:** system(), source(), path traversal or undeclared term expressions fail E_DESIGN_TERM before R evaluates them.

**Contract:** SM07; **owner:** R04; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V032"></a>

### V032: Stable levels and aligned matrices

**Fixture:** Independent fixture with shuffled columns/metadata and explicit factor references.

**Oracle:** Keyed alignment plus direct design-matrix construction.

**Exact assertion:** Preserve aligned coefficients and reference meanings across row/column reorderings.

**Negative case:** Missing reference level or a required nonfinite covariate produces a field-specific error, not silent complete-case sample removal.

**Contract:** SM07; **owner:** R04; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V033"></a>

### V033: Rank and confounding diagnostics

**Fixture:** Eight observations with batch exactly equal to treatment; matched nonconfounded crossed version.

**Oracle:** Independent QR/SVD rank and alias computation.

**Exact assertion:** Reject the confounded model before fit with E_DESIGN_CONFOUNDED and aliased terms; accept the full-rank crossed design.

**Negative case:** A pseudoinverse result must not be published as an identifiable treatment effect for the confounded fixture.

**Contract:** SM07; **owner:** R04; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V034"></a>

### V034: Independent, paired and repeated plans

**Fixture:** Four independent pairs under fixed-subject design and a separate repeated model with known subject IDs.

**Oracle:** Direct full-rank fixed-subject fit and pinned limma duplicateCorrelation/lmFit calls.

**Exact assertion:** Preserve subject identity, pass the correct fixed or random representation and record actual consensus correlation for the repeated fit.

**Negative case:** Missing subject IDs or simultaneous fixed and random encoding of the same subject fail; visits do not inflate independent n.

**Contract:** SM08; **owner:** R04; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V035"></a>

### V035: Numeric contrasts and interaction semantics

**Fixture:** C/U/T group means 10/12/11 and a two-stage treatment effect example.

**Oracle:** d=2, t=−1, r=1; interaction equals the difference of treatment-effect contrasts.

**Exact assertion:** Match d/t/r and r=d+t; calculate the interaction directly from coefficient weights, not from P-value comparisons.

**Negative case:** A contrast with reversed weights/mislabeled numerator or one unknown coefficient fails rather than silently changing sign.

**Contract:** SM07; **owner:** R04; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V036"></a>

### V036: Featurewise estimability and df

**Fixture:** Eight-sample nonconfounded design with one sparse feature losing an entire target group and another losing one observation.

**Oracle:** Observed-row QR rank, n−rank and contrast row-space membership.

**Exact assertion:** Export correct featurewise n/df/estimability and named reasons for ordinary available-case failure; retain excluded rows.

**Negative case:** Do not supply zero effects/P=1 for an unestimable feature or borrow a proDA estimate to fill the limma table.

**Contract:** SM08; **owner:** R04; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V037"></a>

### V037: Exact contrast covariance strategy

**Fixture:** A ≤12-observation nonorthogonal covariate design with unequal weights, one feature-specific missing row and a blocked variant.

**Oracle:** Direct featurewise GLS covariance and contrast-as-coefficient refit, independently authored from the adapter.

**Exact assertion:** Match exact SE/covariance within backend tolerance and record the exact path; demonstrate the approximate shortcut differs on at least one chosen fixture.

**Negative case:** Calling contrasts.fit alone on the general fixture fails this gate even when a convenient complete one-way case happened to agree.

**Contract:** SM08; **owner:** R04; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V038"></a>

### V038: Method-design applicability registry

**Fixture:** Valid limma independent/paired requests; DEqMS weights/block/treat and proDA TMT/block/treat requests.

**Oracle:** Frozen adapter/dispatch tables, independent of hit counts and package availability.

**Exact assertion:** Return each documented eligibility/error type before execution and retain one authoritative primary engine.

**Negative case:** Failure of a DEqMS/proDA capability cannot invoke limma under that engine name; an absent adapter is NOT_RUN, not scientifically INAPPLICABLE.

**Contract:** SM10; **owner:** R04; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V039"></a>

### V039: Frozen analysis plan and matrices

**Fixture:** Independent fixture with two otherwise identical configs differing in one contrast weight/resource hash/exclusion.

**Oracle:** Canonical sorted JSON hash plus referenced actual artifacts.

**Exact assertion:** Write plan, masks, evaluated design/contrasts, hypotheses and families before any fit; each meaningful change produces a different plan hash.

**Negative case:** Intercept a fit request with no valid plan or a changed input hash: reject rather than refreeze after seeing its result.

**Contract:** SM07; **owner:** R04; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V040"></a>

### V040: Design adversarial verification

**Fixture:** Unequal group n, singleton group, absent factor level, repeated-unit imbalance and reordered IDs, each in its own tiny fixture.

**Oracle:** Explicit n/rank/df/estimability calculations and independent reference designs.

**Exact assertion:** Accept only identifiable sufficiently replicated contrasts; exact key alignment makes reordering invariant; singleton inference fails while descriptive QC remains available.

**Negative case:** A missing group level cannot silently disappear from the planned contrast/family coverage.

**Contract:** SM07; **owner:** R04; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

## Boundary

No recovered source/audit changes, private-data transfer, invented license, fake backend or unsupported completion claim. A missing prerequisite is NOT_RUN, not a successful acceptance case. Only the Maintainer updates traceability after independent gate execution.
