# R16e Claude-route receipt, stages 1 and 2 of 2: QC, differential, biomarker, pathway, network and PERMANOVA figures (ADR 0012)

Stage 1 of packet R16e implements the figures FR-206 to FR-214 (V206 to V214): the four QC figures and the five differential figures. Stage 2 implements FR-215 to FR-222 (V215 to V222): the biomarker, pathway, network and PERMANOVA figures, see the stage-2 section at the end. The retrofit of the existing R10b, R13 and R14 figures (part B) and the `proteomics run` wiring are not in this receipt. Implemented 2026-10-08 by the Claude route (Haiku builder, Opus review), at the operator's direction. This is not work by the routes named in AGENTS.md. The receipt means ready for review, not accepted. Independent audit: **NOT_RUN** (the Opus review follows this receipt).

## Identities

- Contract: `specs/017-visualization/contracts/figures.md` (figure IDs, the three outputs, the Info sheet and the statistics table of D-75, the P-value rule of FR-198 with the D-74 threshold note, the error codes). Decisions D-68, D-74, D-75.
- Working-source manifest (stage 1 tree): 28 files, SHA-256 `96a62553424fcfa224e806088c61551a8e55bcf43b72dca69dce3c96dfc74dcd`. Current manifest of the final stage-2 tree: `evidence/working-source-manifest.txt`, 50 files, SHA-256 `87e971d493fad96067422976b4c3841235186136df05d0da403c6740f337932b` (the value in the V215 to V222 traceability rows). It covers the catalogue, the shared HTML writer, the R renderers, the NAMESPACE and DESCRIPTION, the lock, the reason-code catalogue, the R and Python tests and the synthetic fixtures.
- Amendments of this stage: **A-2026-10-01-37** (NAMESPACE export of `render_catalogue_png`; the R16d heatmap hover rule for the row and column identifiers; ggrepel in R Suggests and the lock; reason codes; traceability; tasks), **A-2026-10-01-38** (the R01 export invariant in `test-foundation.R` is updated to the figure functions). The R16e allowlist itself is unchanged.
- `verified_commit`: **null**. Local commits on `claude/visualization`, not pushed.
- Environment: macOS, Python 3.13.15 (`.venv`), R 4.6.1, ggplot2 4.0.3, ggprism 1.0.7, ragg 1.5.2, ggrepel 0.9.8, `LANG=en_US.UTF-8`.

## What is implemented

- `src/proteomics_pipeline/figures/catalogue.py`: one builder per figure. Each figure is registered from its synthetic source table with the R16a registry; the source copy is hashed, and all three outputs are produced from that copy only.
  - PNG: `proteomicsCore::render_catalogue_png`, drawn in R with the R16c Prism theme at 300 dpi (1417 x 1181 px at 120 mm x 100 mm).
  - .pzfx: the R16b writer, with the figure's data tables, the D-75 statistics table where the figure carries inferential statistics (volcano, MA plot, per-protein dot plot), and the Info sheet (statistics `none` otherwise).
  - HTML: the R16d writer, with the trace data built from the same source rows.
  - Refusals: `E_FIGURE_SOURCE_IMPUTED`, `E_FIGURE_SOURCE_MISSING`, `E_FIGURE_SET_MISMATCH`, `E_VENN_K` (the UpSet still runs).
- `r/proteomicsCore/R/figures_catalogue.R`: the PNG renderers.
  - Missingness: a binary observed / missing map (observed light grey, missing black), features ordered by missing fraction, the missing percentage per sample in the column labels.
  - PCA: axis titles "PC1 (xx.x%)" and "PC2 (xx.x%)", the variance from the same prcomp as the scores.
  - MA plot: facet spacing of 8 mm.
  - Clustered DEP heatmap: a group colour bar above the columns in the stable group colours (Control, Acute); the cell colours are computed in R (blue low, white 0, orange high).
  - UpSet: exclusive intersection bars with their counts, a dot matrix whose rows carry the set names, and set-size bars with their totals.
  - Venn (k <= 3): unit circles with centres 1.2 apart; every exclusive region carries its count (the centroid of the grid cells with that exact membership); each circle carries its set name and total.
  - Volcano labels repelled with ggrepel (at most 10, D-74). Brackets are drawn with the FR-198 labels. The clustered heatmap uses Euclidean distance on row z-scores with complete linkage.
- `tests/fixtures/figures/catalogue/`: synthetic source tables from `make_fixtures.py`. The PCA scores and the variance percentages are the output of `stats::prcomp` on the SM06 display matrix. The DEP heatmap tables carry the group of each sample.

## Gates executed

| Command | Exit | Result | Log |
|---|---|---|---|
| `.venv/bin/python -m pytest tests/integration/test_figure_catalogue.py -q -rs -p no:cacheprovider` | 0 | 21 passed | `evidence/python-tests-catalogue.log` |
| testthat `test-figures-catalogue.R` (from a temporary directory) | 0 | 142 expectations, 0 failed | `evidence/r-test-figures-catalogue.log` |
| full R testthat suite, `testthat::test_dir(r/proteomicsCore/tests/testthat)` | 0 | 919 expectations, 0 failed, 0 skipped, 18 files | `evidence/full-r-testthat.log` |
| `R CMD INSTALL --library="$R_LIBS_USER" r/proteomicsCore` | 0 | installed | `evidence/r-install.log` |
| `.venv/bin/python scripts/check_spec_kit.py` | 0 | PASS | `evidence/spec-kit-check.log` |
| `reproduction.check_r_lock('renv.lock')` | 0 | qualified | `evidence/renv-lock-check.log` |
| `.venv/bin/python -m pytest tests/regression tests/unit -q -p no:cacheprovider` | 0 | 85 passed | `evidence/regression-unit.log` |
| full Python suite, `.venv/bin/python -m pytest tests -q -p no:cacheprovider` (background, run on the final tree) | 0 | 515 passed, 0 failed, 0 skipped (`-rs`) | `evidence/full-python-suite.log` |

## Acceptance cases

| Case | Status | Positive oracle (computed independently of the builder) | Negative case |
|---|---|---|---|
| V206 | PASS | per-sample type-7 quantiles of the observed cells, and the observed mask, equal in the HTML box data and in the .pzfx columns (missing cells empty) | an intensity source with no observed mask fails `E_FIGURE_SOURCE_IMPUTED` |
| V207 | PASS | the drawn cell state (1 observed, 0 missing) equals the source mask for every cell; the observed fraction per feature and per sample equals the oracle; features are ordered by missing fraction; the column labels carry the missing percentage; the .pzfx state table equals the mask | a map that draws the missing cells as observed no longer equals the mask (count of mismatches equals the missing count) |
| V208 | PASS | every correlation cell equals the pure-Python Pearson on the cells observed in both samples, within 1e-10; NA below 3 shared cells | a correlation on an imputed matrix (no observed mask) fails `E_FIGURE_SOURCE_IMPUTED` |
| V209 | PASS | scores equal `stats::prcomp` on the median-filled, centred matrix, sign-aligned, within 1e-10; the axis titles equal "PC1 (xx.x%)" and "PC2 (xx.x%)" with the variance from the same prcomp; every ellipse point has Mahalanobis distance squared equal to qchisq(0.95, 2) = 5.9915 within 1e-8 | the group of two (variant) gets no ellipse and the omission is labelled; an ellipse drawn for it without the label fails |
| V210 | PASS | classes from the family q and the sign of the effect equal the HTML classes for all 200 features; the labelled features are exactly the declared subset (10) | a label for a feature absent from the source (`F999`) fails `E_FIGURE_SOURCE_MISSING` |
| V211 | PASS | A = (mean numerator + mean denominator) / 2 and M = mean numerator - mean denominator equal the HTML values within 1e-10, for 2 contrasts of 200 features; M follows the declared numerator | a reversed sign of M does not equal the oracle |
| V212 | PASS | individual values, group means and SDs per protein equal the data; each bracket label equals the FR-198 exact label of its family P (written in the test) | a bracket whose P is not in the family table fails `E_FIGURE_SOURCE_MISSING` |
| V213 | PASS | row order equals `stats::hclust` (Euclidean, complete) on the row z-scores; the z-scores equal the sample-SD formula within 1e-10; only the declared DEP set appears | a non-DEP feature in the source fails `E_FIGURE_SET_MISMATCH` |
| V214 | PASS | exclusive region counts equal the brute-force enumeration of the 2^k - 1 patterns for k = 3 and k = 5; the text labels read from `ggplot_build` of the Venn draw every exclusive count and each set total; the UpSet labels read from `ggplot_build` carry the set names, the intersection counts and the set totals | a Venn for k = 5 fails `E_VENN_K`; the UpSet for k = 5 still produces its three outputs |

## Notes and open points

- **Statistics table.** Per D-75 the volcano and the MA plot carry the statistics table with one row per feature (for the MA, per contrast and feature), and the dot plot carries one row per protein comparison. The QC figures, the clustered heatmap and the UpSet/Venn have no inferential statistics and record `Statistics: none`.
- **Heatmap hover (A-2026-10-01-37).** The R16d check required a `customdata` identifier, which a heatmap does not carry. The check now accepts a heatmap whose hover names the row and column (`%{x}` and `%{y}`), and refuses any other heatmap.
- **Dot-plot hover.** The per-point hover names the protein, the group and the replicate; the P values are shown by the brackets (the FR-198 label), not by the point hover.
- **Venn geometry.** The layout is a standard three-circle Venn with unit radius and centres 1.2 apart. The region label positions are centroids of the grid cells, so an exclusive region with no grid area would be refused (`E_FIGURE_OUTPUT_MISSING`); none occurs in the fixtures.
- **UpSet ordering.** Intersections are ordered by count, then by pattern. The set-size bars are drawn in their own panel with the set order of the declared sets.
- **Clustering through R.** The row order of the clustered heatmap is computed by `stats::hclust` through the R runner from the builder; the test recomputes it with R independently.
- **Reason code `E_FIGURE_OUTPUT_MISSING`.** It is also used for a failed clustering call and for an empty Venn region, which are looser uses of the code than the contract's ("one output is not produced"). The Maintainer may prefer dedicated codes.
- **Font status.** The renderer returns the font status (Arial or the recorded substitution); the builder does not yet write the substitution into the registry (stage 2 wiring).
- **Prism-openability.** The .pzfx files of this stage are not opened in GraphPad Prism. That check stays a manual Maintainer check (NOT_RUN).
- **Visual review.** The renders were reviewed by eye for the volcano, UpSet, Venn, missingness, PCA, MA, and clustered heatmap. No HTML of the new figures was opened in a browser in this stage.
- No private data, no study output. All fixtures are synthetic.

## Acceptance

All gates of stage 1 exited 0, and every case, including its negative case, passed. The R16e stage-1 figures are recorded as ready for review (PASS by the Claude route; independent audit NOT_RUN). Stage 2 (the biomarker, pathway, network and PERMANOVA figures) is recorded in the stage-2 section below; the retrofit and the workflow wiring are not started.

## Full-suite results

- Full R testthat suite on the final tree: 919 expectations, 0 failed, 0 skipped (`evidence/full-r-testthat.log`).
- Full Python suite on the final stage-1 tree: **515 passed, 0 failed, 0 skipped** (run with `-rs`, so any skip would be named; none). An earlier run showed 3 failures: two were caused by the reviewer reinstalling proteomicsCore during the run (they pass in isolation), one was a real missing `fileEncoding = "UTF-8"` in `figures_catalogue.R`, fixed before this run.

## Stage 2 (V215 to V222): biomarker, pathway, network and PERMANOVA figures

Two commits on `claude/visualization`, not pushed: `ab4f5aa` (group 1: biomarker figures V215 to V217) and `8780e13` (group 2: pathway, network and PERMANOVA figures V218 to V222). Claude route (Haiku builder, Opus review), at the operator's direction. The receipt means ready for review, not accepted. Independent audit: **NOT_RUN**.

### What is implemented

- `src/proteomics_pipeline/figures/catalogue.py`: one builder per stage-2 figure, each registered from its synthetic source table and producing the three outputs from the registered copy only (PNG via the R renderer, `.pzfx` with the data tables, the D-75 statistics table where the figure carries inferential statistics, and the Info sheet, HTML via the R16d writer).
- `r/proteomicsCore/R/figures_catalogue.R`: one renderer per stage-2 figure in `render_catalogue_png`. The network layout is `catalogue_network_layout`, internal (not exported), deterministic (modules on a circle of radius 0.6, members on a circle of radius 0.25, modules sorted, no random seed), computed in R and called from the builder for the HTML. The NAMESPACE is unchanged.
- Claim labels are the frozen labels of `post_de_biomarker.R`. The ROC and the CV AUC require `cross_validated_nested`; the panel plot accepts a frozen label; a training-fold ROC (`in_sample`) is refused.
- New code: `E_FIGURE_PERMUTATION_P` (V217). Reused: `E_BIOMARKER_CLAIM`, `E_FIGURE_SOURCE_MISSING`, `E_FIGURE_SET_MISMATCH`. The reason-code catalogue was regenerated with `scripts/maintained/build_release.py codes`.
- Fixtures (synthetic, no study data), generated by `tests/fixtures/figures/catalogue/make_fixtures.py`: pROC `roc` and `ci.se` (seed 20261008, B = 2000), `stats::t.test` (Welch), `stats::phyper` with BH, `stats::cmdscale`, `vegan::betadisper` and `vegan::permutest` (seed 20261003, B = 199), and the R13 PERMANOVA fixture factory (`tests/fixtures/permanova/make_fixtures.py`, loaded by file path). The generator is deterministic; re-running it reproduces every source byte for byte.

### Acceptance cases

| Case | Status | Positive oracle (computed independently of the builder) | Negative case |
|---|---|---|---|
| V215 | PASS | curve points (the set of thresholds) and AUC (pairwise Mann-Whitney) within 1e-10; the drawn band equals `pROC::ci.se` with the recorded seed within 1e-10; the claim label is in the title, the source and the Info sheet | a ROC from training-fold scores (`in_sample`) fails `E_BIOMARKER_CLAIM` |
| V216 | PASS | estimates and intervals equal `stats::t.test` (Welch) within 1e-10; the statistics table P equals the oracle | a panel plot without a claim label fails `E_BIOMARKER_CLAIM` |
| V217 | PASS | observed = mean of the nested-CV repeats; every null value and every repeat drawn; the displayed P equals (k + 1) / (B + 1) = 0.050 with B = 19 shown | a displayed P of 0.10 fails `E_FIGURE_PERMUTATION_P` |
| V218 | PASS | gene ratio = k / n and -log10 q equal within 1e-10; only sets with a foreground hit are plotted; the statistics table P and q equal | SET_E (no foreground hit) requested as enriched fails `E_FIGURE_SET_MISMATCH` |
| V219 | PASS | every cell equals the group mean of the sign-aligned PC1 of the scaled module (`prcomp`, computed in the test) within 1e-10; only declared modules appear | an undeclared module (M4) fails `E_FIGURE_SET_MISMATCH` |
| V220 | PASS | nodes and edges equal the induced subnetwork of the interaction snapshot; layout coordinates are not compared | edge N01-N08 (absent from the snapshot) fails `E_FIGURE_SOURCE_MISSING` |
| V221 | PASS | coordinates equal `stats::cmdscale` within 1e-8 after sign alignment, and the percent of positive eigenvalues within 1e-8; each drawn centroid equals the mean of its own group within 1e-10 | a centroid taken from the full sample differs from the drawn centroid (asserted) |
| V222 | PASS | distances equal `vegan::betadisper` within 1e-10; the PERMDISP P equals `vegan::permutest` with the recorded seed; B and seed appear in the title and the statistics table | a P without its seed fails `E_FIGURE_SOURCE_MISSING` |

### Gates

Before each commit (the gates named in the brief):

| Commit | Catalogue integration tests | testthat `test-figures-catalogue.R` | Regression, unit and text encoding | Spec kit |
|---|---|---|---|---|
| `ab4f5aa` (V215 to V217) | 27 passed | 166 expectations, 0 failed | 91 passed | PASS |
| `8780e13` (V218 to V222) | 37 passed | 188 expectations, 0 failed | 91 passed | PASS |

On the final tree (after commit `8780e13`, which is the content of the reviewed manifest):

- Full Python suite, `.venv/bin/python -m pytest tests -q -rs -p no:cacheprovider`, run once in the background: **531 passed, 0 failed, 0 skipped** in 3714.55 s (`evidence/stage2-full-python-suite.log`). Stage 1 had 515; the 16 new tests are the V215 to V222 cases.
- Full R testthat suite, `testthat::test_dir('r/proteomicsCore/tests/testthat', ...)`: **966 expectations, 0 failed, 0 skipped**, 18 files (`evidence/stage2-full-r-testthat.log`). Stage 1 had 919.
- renv lock check: qualified (`evidence/stage2-renv-lock-check.log`). The lock is unchanged in stage 2: pROC, vegan and ggrepel were already in it.
- The acceptance ledger: 188 PASS and 34 NOT_RUN (stage 1: 180 and 42). `gate-results.json` lists all commands of both stages with their logs.

### Amendments of this stage

- `packet-ownership.json` maintainer amendments 38 (commit `ab4f5aa`) and 39 (commit `8780e13`): the reason-code catalogue, traceability status and task checkboxes for V215 to V222, the acceptance ledger. The R16e write allowlist itself is unchanged.
- Working-source manifest: the V215 to V222 rows point to the final manifest (87e971d4...). The stage-1 rows (V206 to V214) keep their stage-1 manifest hash; see the doubts below.

### Doubts for review

1. The panel plot accepts any frozen claim label, including `in_sample`. A panel evaluated in sample is a performance estimate, so it is not refused, but a reviewer may want a stricter rule.
2. The Prism table kinds: the contract lists column, XY and grouped tables. The network graph writes a column table (node degree) and an XY table (the edge list as node indices), and the ROC writes its band as two XY tables. The statistics table is omitted for the ROC and the network, which have no inferential statistics; the ROC AUC and its interval are in the Info notes.
3. `E_FIGURE_PERMUTATION_P` is a new code that the contract does not list. The Maintainer may prefer an existing code.
4. The V206 to V214 rows keep the stage-1 manifest hash, which predates the stage-2 changes to the shared catalogue code. The Maintainer should decide whether those rows are re-verified.
5. The commit trailers name Claude Opus 5.5, as the operator's brief asks. The builder for this work is the Haiku route, as the stage-1 receipt also records.
6. Not run: GraphPad Prism opening of any `.pzfx` (a Maintainer check). The `pzfx` R package reads all eight stage-2 `.pzfx` files without error (manual check, not a test). The HTML files were not opened in a browser; the static scan and the hover checks run in the writer.
7. Visual review by eye of the eight stage-2 renders. The captions of the ROC, the CV AUC and the network were shortened after review, and the pROC `progress` argument (deprecated) was removed.
8. The R suite prints warnings that also appear in the stage-1 log (a ggplot2 S3 guide deprecation in V203 and the `review minor 2` messages of the post-DE biomarker tests). They are not new.
