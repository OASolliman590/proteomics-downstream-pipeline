# Feature Specification: Treatment response, equivalence and independent scores

**Phase:** 2. **Packet:** R09. **Status:** 1.2.0-frozen; implementation pending Phase 1 acceptance and separate authorization.

**No p-value on in-sample selected scores.** See [SM23–SM24](../001-downstream-proteomics/contracts/scientific-methods.md#treatment-response-and-scores); descriptive-score output physically lacks inferential columns.

## Scope

Deliver only FR-081–FR-090 for US4. [Packet index](../001-downstream-proteomics/packet-index.md), [ownership](../001-downstream-proteomics/packet-ownership.json), [data model](../001-downstream-proteomics/data-model.md) and [scientific contract](../001-downstream-proteomics/contracts/scientific-methods.md) define the shared interfaces. No implementation dispatch is authorized yet.

## Requirements

- **FR-081 — Explicit disease-treatment-residual axes:** The system MUST export coherent axes and covariance and verify r=d+t with preserved contrast signs and model ID.
- **FR-082 — Descriptive reversal without double dipping:** The system MUST export only DescriptiveScore with hypothesis_type=descriptive_score and physically no P/q/test-statistic columns; record the overlap reason.
- **FR-083 — Bounded response categories:** The system MUST assign exact boundary classes, NA for an invalid denominator, crossed_control for RI>1 and overshoot only for RI>1.2; RI=3 is never restoration.
- **FR-084 — Ratio uncertainty contract:** The system MUST match bounded or unbounded/disjoint confidence sets and retain interval kind; descriptive-only emits no interval and no ratio P value.
- **FR-085 — Model-correct residual equivalence:** The system MUST both zero-null tests can be nonsignificant, but only the narrow-interval fixture passes unadjusted alpha=.05 TOST; export p_lower/p_upper/max-P and separate family q.
- **FR-086 — Conjunction claim eligibility:** The system MUST retain components, direction-source hash and joint endpoint in FormalRescueResult only when independence and model eligibility pass.
- **FR-087 — Independent fixed-score artifacts:** The system MUST apply exactly the training transforms, preserve missing-feature policy and emit inferential eligibility only for verified disjoint subjects/units.
- **FR-088 — Exact and Monte Carlo randomization:** The system MUST compute exact k/70 using the frozen tie rule; top-four versus bottom-four gives k=2 and 2/70. MC uses the full-space-with-replacement protocol and (k+1)/(B+1), separately labeled with precision.
- **FR-089 — Score/equivalence reporting semantics:** The system MUST keep descriptive movement, overshoot, equivalence, formal conjunction and independent score test distinct in tables/prose; never generate causal rescue language.
- **FR-090 — Circularity/null/overshoot tests:** The system MUST demonstrate shared-control negative covariance and prevent every in-sample-score P route; classify RI=3 as overshoot and preserve the correct exact denominators.

## Acceptance Scenarios

<a id="V081"></a>

### V081: Explicit disease-treatment-residual axes

**Fixture:** Independent means C/U/T=10/12/11 with covariance diag(.25,.25,.25).

**Oracle:** Linear algebra d=2, t=−1, r=1; Cov(d,t)=−.25, Var(d)=Var(t)=Var(r)=.5.

**Exact assertion:** Export coherent axes and covariance and verify r=d+t with preserved contrast signs and model ID.

**Negative case:** Contrasts from incompatible models or residual weights not equal to d+t fail E_AXIS_INCOHERENT.

**Contract:** SM19; **owner:** R09; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V082"></a>

### V082: Descriptive reversal without double dipping

**Fixture:** A score whose selected features/signs use one or more of the tested contrast samples.

**Oracle:** Set intersection of training/selection and validation biological-unit/subject IDs.

**Exact assertion:** Export only DescriptiveScore with hypothesis_type=descriptive_score and physically no P/q/test-statistic columns; record the overlap reason.

**Negative case:** A request to permute already-selected scores cannot produce IndependentScoreTest or a rescue label.

**Contract:** SM23; **owner:** R09; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V083"></a>

### V083: Bounded response categories

**Fixture:** d=1, t=−RI for RI in [−.1,0,.299999,.3,.799999,.8,1,1.2,1.200001,3]; also d=0 and abs(d)<dmin.

**Oracle:** Literal SM20 class inequalities and r=d×(1−RI).

**Exact assertion:** Assign exact boundary classes, NA for an invalid denominator, crossed_control for RI>1 and overshoot only for RI>1.2; RI=3 is never restoration.

**Negative case:** An unlimited full_reversal class or rounding 1.200001 into near_restoration fails.

**Contract:** SM20; **owner:** R09; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V084"></a>

### V084: Ratio uncertainty contract

**Fixture:** Two bivariate normal contrast estimates with known covariance: one stable denominator and one whose uncertainty crosses zero.

**Oracle:** Independent Fieller quadratic inequality solution using the declared covariance/critical value.

**Exact assertion:** Match bounded or unbounded/disjoint confidence sets and retain interval kind; descriptive-only emits no interval and no ratio P value.

**Negative case:** Dividing separate interval endpoints or truncating an unbounded interval to a finite range fails.

**Contract:** SM20; **owner:** R09; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V085"></a>

### V085: Model-correct residual equivalence

**Fixture:** Residual r=0 with df=20, delta=.2, and SE=.05 versus SE=1.

**Oracle:** Direct Student-t CDF one-sided tests and qt(.95,20) for the 90% interval.

**Exact assertion:** Both zero-null tests can be nonsignificant, but only the narrow-interval fixture passes unadjusted alpha=.05 TOST; export p_lower/p_upper/max-P and separate family q.

**Negative case:** Nonsignificance alone, 95%-effect-CI substitution or recycled zero-null q cannot produce equivalence.

**Contract:** SM21; **owner:** R09; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V086"></a>

### V086: Conjunction claim eligibility

**Fixture:** Externally fixed positive disease direction, disease/treatment margins, residual margin and four valid one-sided component P values [.01,.02,.03,.04].

**Oracle:** Intersection-union joint P=max=.04, followed by the declared conjunction-family adjustment.

**Exact assertion:** Retain components, direction-source hash and joint endpoint in FormalRescueResult only when independence and model eligibility pass.

**Negative case:** Direction selected on tested disease samples or intersection of separately BH-filtered lists cannot claim formal_rescue/FDR.

**Contract:** SM22; **owner:** R09; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V087"></a>

### V087: Independent fixed-score artifacts

**Fixture:** A two-feature score manifest with fixed weights/centers/scales and training subjects disjoint from four tested pairs.

**Oracle:** Direct application of the fixed transform and participant-set disjointness check.

**Exact assertion:** Apply exactly the training transforms, preserve missing-feature policy and emit inferential eligibility only for verified disjoint subjects/units.

**Negative case:** A reused subject at another visit, unknown selection cohort or refitted validation transform degrades to descriptive_score with no P columns.

**Contract:** SM23; **owner:** R09; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V088"></a>

### V088: Exact and Monte Carlo randomization

**Fixture:** Eight fixed scores [1,2,3,4,5,6,7,8], every 4-vs-4 allocation and a separate fixed-seed Monte Carlo draw stream.

**Oracle:** Independent itertools combinations and absolute mean differences; exactly 70 allocations, observed split once.

**Exact assertion:** Compute exact k/70 using the frozen tie rule; top-four versus bottom-four gives k=2 and 2/70. MC uses the full-space-with-replacement protocol and (k+1)/(B+1), separately labeled with precision.

**Negative case:** 3/71 for the exhaustive k=2 fixture fails; unrestricted permutation of paired subjects or exclusion of the observed allocation from an unbounded with-replacement MC space fails.

**Contract:** SM24; **owner:** R09; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V089"></a>

### V089: Score/equivalence reporting semantics

**Fixture:** A descriptive near_restoration feature failing TOST, an overshoot feature and an in-sample-selected score.

**Oracle:** Typed input objects and SM20–SM24 vocabulary, with a separate eligible score test for comparison.

**Exact assertion:** Keep descriptive movement, overshoot, equivalence, formal conjunction and independent score test distinct in tables/prose; never generate causal rescue language.

**Negative case:** A near_restoration class cannot be labeled statistically equivalent without its passing eligible EquivalenceResult.

**Contract:** SM25; **owner:** R09; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V090"></a>

### V090: Circularity/null/overshoot tests

**Fixture:** The V081 covariance fixture, V082 overlap fixture and independent 70-allocation arithmetic cases with extreme counts 2 and 6.

**Oracle:** Exact algebra and enumeration; 2/70 and 6/70 are fixture arithmetic, not constants returned by production.

**Exact assertion:** Demonstrate shared-control negative covariance and prevent every in-sample-score P route; classify RI=3 as overshoot and preserve the correct exact denominators.

**Negative case:** Historical 3/71 or 7/71, hardcoded answers for arbitrary input, or a selected-score null permutation labeled confirmatory fails.

**Contract:** SM19; **owner:** R09; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

## Boundary

No recovered source/audit changes, private-data transfer, invented license, fake backend or unsupported completion claim. A missing prerequisite is NOT_RUN, not a successful acceptance case. Only the Maintainer updates traceability after independent gate execution.
