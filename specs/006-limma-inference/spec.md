# Feature Specification: Core limma inference and multiplicity

**Phase:** 1. **Packet:** R05. **Status:** 1.2.0-frozen; implementation pending prerequisite acceptance and explicit authorization.

## Scope

Deliver only FR-041–FR-050 for US2. [Packet index](../001-downstream-proteomics/packet-index.md), [ownership](../001-downstream-proteomics/packet-ownership.json), [data model](../001-downstream-proteomics/data-model.md) and [scientific contract](../001-downstream-proteomics/contracts/scientific-methods.md) define the shared interfaces. No implementation dispatch is authorized yet.

## Requirements

- **FR-041 — Observed-data limma adapter:** The system MUST match effects, moderated SE/df/statistic/P within frozen backend tolerance; retain original NA and no newly imputed primary cells.
- **FR-042 — Exact sparse and weighted contrasts:** The system MUST match exact effects/SE/P and block/weight provenance through export, not only inside a helper.
- **FR-043 — Robust/trend diagnostics:** The system MUST save requested/actual settings and prior/posterior variance diagnostics and match reference values; report robust eBayes as hyperparameter robustness, not robust sample regression.
- **FR-044 — Effect-threshold TREAT endpoint:** The system MUST preserve distinct raw P/q and hypothesis_type for protein_treat versus protein_zero_null; display abs(effect) filtering never changes either P or its family.
- **FR-045 — Backend-correct confidence intervals:** The system MUST match 95% effect intervals and report actual CI method/df; unsupported/unavailable uncertainty stays NA with a reason.
- **FR-046 — Declared multiple-testing families:** The system MUST pool every planned primary zero-null contrast-feature test once, keep secondary separate and report complete planned/eligible/finite counts.
- **FR-047 — Omnibus and interaction outputs:** The system MUST export distinct omnibus hypothesis/family and all planned pairwise/interaction endpoints without unadjusted post-hoc claims.
- **FR-048 — Complete typed differential tables:** The system MUST retain every planned endpoint or explicit row-level exclusion/nonestimability; export hypothesis_type/family_id/engine/estimable/n_obs_by_required_group and mark a numerical-failure family incomplete.
- **FR-049 — Biological-unit influence diagnostics:** The system MUST report whole-subject influence with correct relevant omission counts and descriptive effect changes; do not turn it into held-out validation.
- **FR-050 — Core inference golden suite:** The system MUST pass all reference comparisons; zero discoveries remain COMPLETED when required work finishes; show Cov(d,t)=−.25 and Var(d+t)=.5, not 1; Phase 1 emits no score P-value artifact.

## Acceptance Scenarios

<a id="V041"></a>

### V041: Observed-data limma adapter

**Fixture:** Complete and two-NA 8×12 matrices with declared log2 scale and one primary zero-null model.

**Oracle:** Separate pinned lmFit/eBayes reference using exactly the declared fitting universe/settings.

**Exact assertion:** Match effects, moderated SE/df/statistic/P within frozen backend tolerance; retain original NA and no newly imputed primary cells.

**Negative case:** An adapter that silently fills NA or uses a different fit universe fails even if selected P values look similar.

**Contract:** SM09; **owner:** R05; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V042"></a>

### V042: Exact sparse and weighted contrasts

**Fixture:** The V037 general-design fixtures fitted through the full limma adapter.

**Oracle:** Direct reparameterized limma coefficients plus independent featurewise unscaled covariance.

**Exact assertion:** Match exact effects/SE/P and block/weight provenance through export, not only inside a helper.

**Negative case:** Dropping weights, missing rows or block correlation to obtain a complete-matrix shortcut fails.

**Contract:** SM08; **owner:** R05; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V043"></a>

### V043: Robust/trend diagnostics

**Fixture:** A ≤20-feature variance/abundance trend fixture with a hypervariable feature and explicit trend/robust flags.

**Oracle:** Separate direct eBayes calls with the same flags and package versions.

**Exact assertion:** Save requested/actual settings and prior/posterior variance diagnostics and match reference values; report robust eBayes as hyperparameter robustness, not robust sample regression.

**Negative case:** Missing statmod with robust=true fails dependency preflight; silently changing robust=false is forbidden.

**Contract:** SM09; **owner:** R05; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V044"></a>

### V044: Effect-threshold TREAT endpoint

**Fixture:** Effects below/on/above tau=.5 in an otherwise eligible synthetic fit, with separately declared zero-null and TREAT families.

**Oracle:** Direct limma::treat with lfc=.5 and separate eBayes zero-null call.

**Exact assertion:** Preserve distinct raw P/q and hypothesis_type for protein_treat versus protein_zero_null; display abs(effect) filtering never changes either P or its family.

**Negative case:** A zero-null q table filtered by abs(logFC)≥.5 cannot be exported as ProteinTreatResult.

**Contract:** SM11; **owner:** R05; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V045"></a>

### V045: Backend-correct confidence intervals

**Fixture:** Known coefficient SE, posterior variance and inference df from a direct limma fit.

**Oracle:** effect ± qt(.975,df_inference)×moderated_SE, computed independently.

**Exact assertion:** Match 95% effect intervals and report actual CI method/df; unsupported/unavailable uncertainty stays NA with a reason.

**Negative case:** A normal-quantile or residual-df interval substituted for the declared moderated-t interval fails its numerical oracle.

**Contract:** SM13; **owner:** R05; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V046"></a>

### V046: Declared multiple-testing families

**Fixture:** Raw primary P values [.01,.02,.5,.8] spanning two contrasts; separate secondary P=.001.

**Oracle:** BH primary q=[.04,.04,2/3,.8]; secondary one-test q=.001; BY uses the independent harmonic factor.

**Exact assertion:** Pool every planned primary zero-null contrast-feature test once, keep secondary separate and report complete planned/eligible/finite counts.

**Negative case:** Computing primary BH separately per contrast or adding secondary P into primary fails expected q/family membership.

**Contract:** SM12; **owner:** R05; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V047"></a>

### V047: Omnibus and interaction outputs

**Fixture:** Two-coefficient omnibus and a direct treatment×stage contrast with all primary/auxiliary effects requested.

**Oracle:** Pinned joint F test/reduced-model reference and explicit interaction vector.

**Exact assertion:** Export distinct omnibus hypothesis/family and all planned pairwise/interaction endpoints without unadjusted post-hoc claims.

**Negative case:** One significant and one nonsignificant treatment contrast cannot produce an interaction-significant flag.

**Contract:** SM11; **owner:** R05; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V048"></a>

### V048: Complete typed differential tables

**Fixture:** Eight features×three contrasts with one coverage exclusion and one deliberately failed eligible numerical test.

**Oracle:** Cartesian planned keys and data-model field types.

**Exact assertion:** Retain every planned endpoint or explicit row-level exclusion/nonestimability; export hypothesis_type/family_id/engine/estimable/n_obs_by_required_group and mark a numerical-failure family incomplete.

**Negative case:** A missing eligible row, NA→0 coercion or native-q-as-central-q fails table verification.

**Contract:** SM12; **owner:** R05; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V049"></a>

### V049: Biological-unit influence diagnostics

**Fixture:** Three subjects each with two visits and two technical observations (12 observations), aggregated before fitting.

**Oracle:** Explicit omission of all records of one subject and direct model refit.

**Exact assertion:** Report whole-subject influence with correct relevant omission counts and descriptive effect changes; do not turn it into held-out validation.

**Negative case:** Dropping only one injection/visit or omitting an unrelated group in a relevant-group count fails.

**Contract:** SM13; **owner:** R05; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

<a id="V050"></a>

### V050: Core inference golden suite

**Fixture:** Complete/sparse/paired/blocked direct-reference fixtures; a symmetric equal-group-means null matrix; a C/U/T mean-covariance fixture diag(.25,.25,.25).

**Oracle:** Independent limma references; algebra using d=[−1,1,0], t=[0,−1,1], r=[−1,0,1].

**Exact assertion:** Pass all reference comparisons; zero discoveries remain COMPLETED when required work finishes; show Cov(d,t)=−.25 and Var(d+t)=.5, not 1; Phase 1 emits no score P-value artifact.

**Negative case:** A shared-control opposition claim labeled independent confirmation, a score P column in Phase 1, or zero discoveries mapped to FAILED fails this gate.

**Contract:** SM19; **owner:** R05; **evidence:** pending-after-freeze, none recorded. Numeric comparison uses [frozen tolerances](../001-downstream-proteomics/validation-strategy.md#numeric-tolerances).

## Boundary

No recovered source/audit changes, private-data transfer, invented license, fake backend or unsupported completion claim. A missing prerequisite is NOT_RUN, not a successful acceptance case. Only the Maintainer updates traceability after independent gate execution.
