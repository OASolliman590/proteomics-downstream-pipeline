# ADR 0006 — Operator-authorized PERMANOVA scope amendment

Status: accepted by the operator (Omar, Maintainer) on 2026-10-01; recorded by the Claude (Opus) route acting for the Maintainer role. This is an explicit scope amendment, not a silent expansion.

## Context

The frozen v1.2.0 kit covers protein-wise inference (limma and qualified adapters), enrichment and response. The operator completed a real analysis in which multivariate group separation (PERMANOVA with PERMDISP) answered a reviewer question that protein-wise tests cannot, and asked that it become a maintained, tested capability. AGENTS.md forbids scope expansion without authorization and forbids renumbering identities.

## Decision

- Add packet **R13** (slice `014-multivariate-permanova`, Phase 2, v0.2-qualified-methods) with **FR-121–FR-130, T121–T130, V121–V130**, appended after the existing ranges. FR-001–FR-120, T001–T120 and V001–V120 keep their identity.
- Bump the scientific methods contract to **v1.2.0** with SM26–SM30; SM01–SM25 are unchanged.
- Dispatch R13 after accepted R05 and R10a, serially, before the R06 ∥ R07 ∥ R09 group (operator order). It is not part of the Phase 1 milestone.
- Shared-interface amendment A-2026-10-01-01: add `permanova` → `proteomics_pipeline.permanova_service` to the fixed Python capability map and `permanova` → `permanova_stage` to the fixed R dispatch map; add the optional `multivariate` configuration block (Phase 1 configurations cannot enable it); add optional R Suggests (vegan, permute and the other packages later packets use).
- The design follows the operator's reference workflow (z-scored Euclidean distance, adonis2 with marginal SS, covariate terms permuted within the primary group, group term from unrestricted permutations, PERMDISP beside every test, mean-R² identity, random equal-size set null for same-data panels). No data, numbers or figures from that private study enter the repository; all fixtures are synthetic.

## Consequences

The kit now has 130 requirements across 13 slices; `scripts/check_spec_kit.py` checks the amended counts and that the first 120 identities are unchanged. Missing vegan/permute is NOT_RUN, never INAPPLICABLE. PERMANOVA results are descriptive multivariate separation, not validation of a panel.

## Verification

V121–V130 (`tests/scientific/test_permanova.py`, `r/proteomicsCore/tests/testthat/test-permanova.R`); receipt in `docs/validation/014-multivariate-permanova/`.
