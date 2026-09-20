# Research decisions for implementation

Research date: 2026-09-12. Detailed primary-source/API notes are in [scientific-methods-research.md](scientific-methods-research.md), authored from opened official sources and primary methods papers. The scientific contract translates these into project decisions; current online documentation does not establish the historical package versions.

| Decision | Rationale | Alternative and boundary |
|---|---|---|
| Python orchestration + R scientific package | Builds on the recovered limma/R work while separating filesystem/reporting from scientific inference | No new statistical reimplementation or web platform |
| Preserve processed input by default | Avoids double log transformation/normalization and retains provenance uncertainty | Explicit eligible normalization/sensitivities supported |
| limma as core eligible linear backend | Established moderation/contrasts; supports modeled covariates and qualified correlation handling | DEqMS/proDA are assay-qualified explicit alternatives, not hit-count fallbacks |
| Exact contrast refitting for complex missing designs | Standard contrasts.fit may have approximate SE under nonorthogonal missing/weighted designs | Exact one-way shortcut only when its conditions are proven |
| TREAT only for actual effect-threshold question | Zero-null significance plus an observed cutoff does not test a true minimum effect | DEqMS/proDA unsupported threshold combinations fail explicitly |
| Explicit multiplicity families | Each adjusted P must identify its hypothesis universe and primary/sensitivity scope | Per-contrast family is allowed only as an explicit narrower claim |
| CAMERA / ROAST according to covariance model | CAMERA's arguments do not implement duplicateCorrelation blocks; ROAST has a different self-contained null | Preranked fgsea remains optional exploratory analysis |
| Label-independent gene representative | Avoids choosing the most extreme protein statistic per gene | Prespecified aggregation sensitivity is separately modeled and reported |
| Descriptive reversal by default | Shared untreated estimates induce negative covariance; selected-score testing is circular | Independent fixed-score validation and properly defined equivalence supported |
| Frozen resources and per-stage versions | Old environment/session capture cannot reproduce gene sets or loaded libraries | Explicit resource acquisition is separate from offline analysis |
| Spec-of-specs with twelve slices | Broad maturity work needs independently testable bounded implementation and review | No single opaque implementer run is trusted without slice gates |

## Methodology sources

- [limma contrasts.fit](https://raw.githubusercontent.com/bioc/limma/master/man/contrasts.fit.Rd) and [limma manual](https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf): model and contrast applicability; source/package version must be pinned during implementation.
- [DEqMS official vignette](https://bioconductor.org/packages/release/bioc/vignettes/DEqMS/inst/doc/DEqMS-package-vignette.html): count covariates and count-adjusted statistics.
- [proDA introduction](https://const-ae.github.io/proDA/articles/Introduction.html) and [test_diff API](https://const-ae.github.io/proDA/reference/test_diff.html): dropout model assumptions and supported hypotheses.
- [CAMERA implementation](https://raw.githubusercontent.com/bioc/limma/master/R/geneset-camera.R) and [ROAST documentation](https://raw.githubusercontent.com/bioc/limma/master/man/roast.Rd): finite inputs, covariance support, nulls and output semantics.
- [Spec Kit](https://github.github.com/spec-kit/) and [spec-of-specs convention](https://github.github.com/spec-kit/concepts/spec-of-specs.html): constitution, per-feature spec/plan/tasks and bidirectional roadmap links. This project follows the artifact format without claiming CLI integration installation.

## Decisions that implementation must resolve by execution

Choose a compatible current R/Bioconductor environment and exact Python/R package versions, validate installed signatures, build locks and reproduce references. Packet scope and ordering are governed by `packet-index.md`, not by an external CLI probe. Do not label an untested environment or adapter as working. Missing historical resource versions constrain exact legacy reproduction; use a separately versioned maintained analysis and document the difference.

No remaining product preference requires a clarification round for implementation. Unknown experimental metadata cannot be invented; the system must support explicit unknowns and constrained reports. Ownership/licensing/public deployment remain outside this implementation mandate.
