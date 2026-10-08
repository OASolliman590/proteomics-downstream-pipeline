# ADR 0010: Adaptive policy for post-DE modules and PERMANOVA

Status: accepted — directed by the Maintainer (Omar) on 2026-10-05. Recorded by the operator-authorized Claude route (ADR 0007). Decisions D-59 and D-60 in [decisions.md](../../specs/001-downstream-proteomics/decisions.md).

## Context

The Maintainer's direction: *"the policy choices need to be designed to fit any type of data so pipeline doesn't be a hold up for any analysis. wrap up the pipeline please"*. Several post-DE rules and PERMANOVA rules refused a whole module, or rejected the whole run, for situations where part of the analysis was still valid or a valid alternative existed. Examples: co-abundance on repeated designs, small cohorts, an over-budget biomarker evaluation, adjusted Pearson/Spearman, an aliased pooled phenotype, an invalid set definition among valid ones, and an unbalanced subject layout in PERMANOVA.

## Decision

1. **Never hold up.** Post-DE modules and PERMANOVA adapt to the data and the declaration. They never block the primary differential-abundance results or other modules. An optional module whose declaration is invalid becomes INAPPLICABLE with its typed reason; the run continues.
2. **Refuse only what is invalid, as narrowly as possible.** A refusal is reserved for an output that would be statistically invalid or that cannot be computed. Only that piece is refused (INAPPLICABLE with reason), and the valid alternative runs automatically where one exists. Examples: complete-case analysis instead of phenotype imputation; within-group analyses instead of an aliased pooled estimate; the leakage-safe form of a leaky declaration; the measured-universe null instead of a genome-wide null.
3. **Adaptations are never silent.** Every adaptation records the analysis, the item, what was requested, what was used, and the reason. Records go into the planner decision (`plan.json`), the module's `eligibility.json` and settings, a `W_POST_DE_ADAPTED` (or `W_PERMANOVA_ADAPTED`) warning, and the R14f report table. Labels follow what was actually computed (`correlation_used`, `row_transform`, `exploratory`, the adapted `outer_k`/`B` in the result tables). Data-driven thresholds keep a tier: the recommended level runs normally; below it the analysis runs as exploratory (recorded); below the computable minimum it is refused, with the number and the reason stated.
4. **`required` stays explicit strictness.** A required module or analysis that cannot run as declared still rejects the plan. This covers invalid declarations (grammar, leakage, unknown references, dropped declared content) and computable-minimum refusals. Post-DE modules and PERMANOVA default to optional.
5. **Relation to ADR 0005 (no silent method fallback).** An adaptation is not a silent substitution. The method used is named in the results and the record, the requested method is preserved beside it, and the claim label reflects what was computed.

## Consequences

The scientific-methods contract (SM27, SM28, SM33–SM36, SM39–SM41) and `contracts/post-de.md` are amended accordingly (D-59). The refusal audit in D-60 classifies every refusal code of R13 and R14a–R14f as *validity*, kept and narrowed, or *policy*, converted to adapt-with-record. The acceptance tests that expected the old refusals now assert the adaptation records, and a required variant still asserts the refusal. The V167 generality matrix gained cells for the previously refused situations. Calibration thresholds are test-time gates only and are unchanged.
