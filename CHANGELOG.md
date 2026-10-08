# Changelog

No version of the maintained successor has been released. No license has been selected; licensing and publication are operator decisions (AGENTS.md). Statuses below come from executed evidence ([acceptance ledger](docs/validation/acceptance.json)). They are self-verified by the operator-authorized Claude route unless an independent audit is named.

## Unreleased — branch `claude/post-de`

### Added
- Post-differential analysis (spec 015, ADR 0009, Phase 4): sets and concordance (R14a), robustness and sensitivity (R14b), protein–phenotype association (R14c), leakage-safe biomarker evaluation (R14d), co-abundance modules and interaction networks (R14e), and the integrated eligibility/dependency report with the generality matrix (R14f). V131–V167 PASS on synthetic fixtures; independent audit NOT_RUN.
- Post-DE calibration evidence (nested-CV AUC and permutation P under the null, a leaky reference, the connectivity null) in `docs/validation/015-post-de-analysis/R14f/calibration/`.
- Independent review of 2026-10-05 (ACCEPT WITH FIXES) resolved: `docs/validation/review-2026-10-05`.
- Release tooling (R12): hygiene check and versioned successor export (`scripts/maintained/build_release.py`), Maintainer-local legacy regression (`scripts/maintained/legacy_regression.py`), methods and reason-code documentation, and the [Maintainer runbook](docs/validation/013-release/MAINTAINER_RUNBOOK.md). V115–V117 PASS; V111–V114 and V118–V120 NOT_RUN (Maintainer-only gates, operator decisions or final review).

### Changed (Maintainer-directed, 2026-10-05)
- Adaptive post-DE and PERMANOVA policy (ADR 0010, D-59). Modules adapt to the data and never hold up the DE results. Only invalid or non-computable pieces are refused, and the valid alternative runs automatically. Every adaptation is recorded with requested vs used and the reason. Co-abundance supports paired and repeated designs and small cohorts (exploratory). The biomarker compute guard scales work down by default (`compute_policy: adapt`). Adjusted Pearson/Spearman run as partial correlations. Aliased pooled phenotypes run within-group. PERMANOVA handles unbalanced and mixed subject designs. The refusal audit is D-60.
- All maintained text outputs are written with LF on every platform (D-61). The report scatter has padding, a zero line and marked identical points (D-62).

### Fixed
- PR #2 CI (run 37244850636): full-history checkout for V116 with a specific message on shallow clones; an idempotent CRLF conversion in the line-ending test; `.gitattributes` checks text out as LF.
- Re-review follow-ups of 2026-10-03: platform-independent plan hash, PERMANOVA empty-group pairs and exact permutation accounting, V059 setup failures, UTF-8 portability (D-42, D-43).
- A model with exactly one eligible feature no longer fails in the limma moderation diagnostics (A-2026-10-01-19).

### Decisions
- The Maintainer accepted the ADR 0008 contract amendments (SM05, SM27, SM28) on 2026-10-04 (D-50) and the SM27 amendment D-43 on 2026-10-05 (D-58).

## Merged into `main` (PR #1, 1578bc8)

- R02–R11 (intake through reproduction and validation tooling) and R13 PERMANOVA/PERMDISP, with the fixes from the independent audit of 2026-10-02 (ADR 0008). V107 cross-platform CI passed on Ubuntu and Windows (run 37087074374).
- R01 foundation (accepted 2026-09-21).
