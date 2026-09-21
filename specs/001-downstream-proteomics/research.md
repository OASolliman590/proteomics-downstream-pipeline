# Research/decision record — candidate hardening

Updated 2026-09-20 because the normative contracts changed. The [scientific methods contract v1.1.0](contracts/scientific-methods.md) is authoritative; the older scientific-methods-research.md is retained background at the recovered baseline, not a competing set of executable rules. The immutable audit is not rewritten.

## What changed and why

The candidate separates phase milestones, makes limma the explicit default, prevents selected-score inference by default, freezes enrichment dispatch and forbids silent method fallback. [ADRs 0001–0005](../../docs/adr/0001-phase-milestones.md) record these decisions. The schema now declares source scale, organism, primary design/engine/hypothesis, response/score mode and multiplicity families explicitly; the preimplementation draft migration is in semantic-validation.md. R10a's thin offline HTML retains FR-091–094 rather than weakening the HTML requirement to a Markdown-only placeholder.

The exact-contrast policy uses a reparameterized direct fit as a stable correctness target. Current limma documentation also describes an exact fit-time contrast interface; a later implementation may qualify a mathematically equivalent exact route against independent references. No current documentation page alone establishes that the future pinned installation supports a given argument. The restricted shortcut rule is intentionally narrower than all mathematically possible exact cases.

## Primary-source verification register

These are source documents inspected for this pass, not tested package pins or offline resource snapshots. Contract limits such as .3/.8/1.2 response classes, required rotations and release-calibration thresholds are prespecified product policies; they are not attributed to the sources below as universal biological standards.

| Source | Supported point | Contract use |
|---|---|---|
| [limma package reference manual](https://bioc.r-universe.dev/limma/doc/manual.html) | contrasts.fit can approximate SEs with nonorthogonal missing/weighted designs; direct contrast fitting avoids that approximation. | SM08; actual installed API still requires independent qualification. |
| [limma release reference manual](https://bioc-release.r-universe.dev/limma/doc/manual.html) | CAMERA requires finite expression and exposes estimated versus fixed inter-gene correlation; the documented methods have distinct nulls/settings. | SM16; paired CAMERA is deliberately disallowed by this kit. |
| [DEqMS author/package reference manual](https://bioc.r-universe.dev/DEqMS/doc/manual.html) | outputResult separates count-adjusted sca.t/sca.P.Value from ordinary limma output and is downstream of spectraCounteBayes. | SM10; do not use original P.Value as DEqMS raw P. |
| [DEqMS official vignette mirror](https://bioconductor.statistik.tu-dortmund.de/packages/3.8/bioc/vignettes/DEqMS/inst/doc/DEqMS-package-vignette.html) | Explicit count-linked fitting/output sequence and diagnostics. This is a historical package vignette, not a selected current runtime. | Backend golden-call design, not a version pin. |
| [proDA author introduction](https://const-ae.github.io/proDA/articles/Introduction.html) | LFQ probabilistic dropout modeling without filling missing values. | SM10 input eligibility and uncertainty labeling. |
| [proDA test_diff reference](https://const-ae.github.io/proDA/reference/test_diff.html) | Coefficient Wald and reduced-model tests have native output and df semantics. | Preserve native statistics; no borrowed limma df or invented TREAT. |
| [msigdbr author reference](https://igordot.github.io/msigdbr/reference/msigdbr.html) | Source database species and output target species differ; human genes can be converted to model-organism counterparts. | SM14 provenance and ortholog-projection label, including human→rat. |
| [msigdbr author introduction](https://igordot.github.io/msigdbr/articles/msigdbr-intro.html) | Ortholog mapping carries evidence and source/target distinctions. | Snapshot manifest/mapping-loss requirements. |

Package-name/version strings observed on a documentation site are not a tested environment. No gene-set snapshot is pinned by a URL alone. Runtime production never invokes msigdbr; explicit separate resource preparation may build local snapshots with actual bytes/hash/source/version/terms. Repository licensing remains unselected regardless of the license of an upstream package.

## Unchanged evidence boundary

The known historical circularity, covariance, permutation, overshoot, Negr1, orthology, representative-rank, stochastic-imputation and hardcoded-finalizer findings remain findings from the supplied prompt/recovered audit. This pass does not claim to have rerun historical R or reproduced private numerical results. Synthetic contrast/permutation arithmetic checks only demonstrate contract oracles. Public source, private archive, resource availability and software calibration are separate evidence categories.
