# Feature Specification: Visualization outputs (Prism files, Prism-style PNG, interactive HTML)

**Phase:** 6 (v1.3-visualization). **Packets:** R16a–R16e (dispatch units). **Status:** Maintainer scope decision of 2026-10-08 ([ADR 0012](../../docs/adr/0012-visualization-outputs.md)); drafted, not yet implemented. It adds FR-195–FR-222, T195–T222 and V195–V222 and does not renumber any existing identity.

## Scope

Deliver FR-195–FR-222 for one journey:
- **US12:** the analyst receives every figure as three outputs: a GraphPad Prism `.pzfx` file the supervisor can restyle, a finished Prism-style PNG at 300 dpi for print, and a standalone interactive HTML file for offline review. All figure areas are in the first build: QC and PCA; differential; biomarker and ROC; pathways and networks; PERMANOVA PCoA and PERMDISP.

Figures are generated only from saved figure-source tables; no figure is computed from a source that the report does not also write. The figure-level rules are SM25 and the scientific rules of the figure's analysis (SM26–SM30, SM31, SM34, SM36, SM39, SM40). Configuration, output layout and figure IDs are in [contracts/figures.md](contracts/figures.md).

This slice is additive. The existing R10b, R13 and R14 figure outputs keep their current files and names. The new catalogue lives in its own directory (see the contract).

### Normative principles (for every packet in this slice)

1. **Three outputs per figure.** Every catalogue figure produces its `.pzfx`, its 300 dpi PNG and its HTML. An output that cannot be produced is reported as NOT_RUN with its reason, and the other outputs of that figure still run. A missing output is never hidden.
2. **Source-first.** Each output is checked against the figure's source table. No output is drawn from a different table, and no value is rounded before it is written.
3. **No network.** Each HTML file is standalone. plotly.js is vendored, pinned and hash-checked. The HTML makes no request to any URL.
4. **Declared style.** Colours, font, widths, dpi and annotation style come only from `report.figures` (and the contract defaults). Group colours are stable across every figure of a run.
5. **Manual checks are separate.** Opening a `.pzfx` in GraphPad Prism is a Maintainer check, recorded as NOT_RUN in automated tests.

## Requirements

### R16a — Style system and figure registry

- **FR-195 — `report.figures` configuration:** The system MUST accept an optional `report.figures` block with `formats` (default `["png", "pzfx", "html"]`), `dpi` (300; any other value for the PNG output is refused), `journal_width_mm` (85, 120 or 180; default 120, D-74), `palette`, `group_colours`, `font` and `p_annotation`, and MUST refuse unknown keys and unknown formats. Existing `report.figure_formats` keeps its meaning and is never merged silently into `report.figures`.
- **FR-196 — Colourblind-safe palette with stable group colours:** The system MUST colour groups from the Okabe–Ito palette in declared group order by default, or from an explicit `group_colours` mapping. The same group MUST have the same colour in every figure of a run, independent of figure order.
- **FR-197 — Figure registry from source tables:** The system MUST register each catalogue figure ID with its source table path and SHA-256 and MUST generate a figure only from its registered source table. A figure request without its source table MUST be refused with no output.
- **FR-198 — P-value annotation:** The system MUST annotate comparisons either as stars or as exact values, as declared. The annotation MUST be computed from the P value in the source table by the contract rule, and no annotation may appear without its P value in the source.

### R16b — Prism `.pzfx` writer

- **FR-199 — Column tables for group comparisons:** The system MUST write each group comparison as a Prism column table with one column per group containing the individual replicate values exactly as in the source table, so that Prism computes means and SD on those values.
- **FR-200 — XY tables:** The system MUST write the volcano (effect against −log10 P) and ROC (1 − specificity against sensitivity) data as Prism XY tables, one point per source row, with values exactly as in the source table.
- **FR-201 — Grouped tables for multi-protein panels:** The system MUST write multi-protein panels as Prism grouped tables with one row per replicate and one data set per protein and group, with values exactly as in the source table.

### R16c — Static Prism-style PNG

- **FR-202 — Prism-style theme:** The system MUST render each catalogue figure as a PNG with ggplot2 and the ggprism package using the Prism-style theme: white background, thick black axes, outward ticks, no gridlines, Arial (or a recorded substitution), and mean ± SD drawn over the individual points.
- **FR-203 — PNG dimensions and data fidelity:** The system MUST export the PNG at the declared dpi, with pixel width equal to the width in inches multiplied by dpi, rounded down to an integer (the raster device truncates), and its plotted data MUST equal the source table as reported by `ggplot_build`.

### R16d — Interactive Plotly HTML

- **FR-204 — Standalone HTML with exact data and hover fields:** The system MUST write one standalone HTML file per figure whose embedded trace data equal the source table and whose hover text includes the identifying fields of the figure (group, feature or set identifier, and the P or q shown).
- **FR-205 — Offline, pinned and licensed plotly.js:** The system MUST embed plotly.js from the vendored copy only after checking its SHA-256 against the pinned value. The HTML MUST make no network request and contain no external resource reference. The MIT licence notice for plotly.js MUST be included with the vendored file and in the HTML.

### R16e — Figure catalogue

- **FR-206 — QC: intensity distributions:** The system MUST plot per-sample log2 intensity distributions from observed cells only, with the observed mask in the source table.
- **FR-207 — QC: missingness heatmap:** The system MUST plot a missingness heatmap of observed fraction per feature and sample, with missing cells shown as missing and never as zero.
- **FR-208 — QC: sample correlation heatmap:** The system MUST plot pairwise sample correlations computed on genuinely observed cells with the coverage rule of SM06, and MUST NOT correlate an imputed matrix.
- **FR-209 — QC: PCA with ellipses:** The system MUST plot PCA scores with a 95 % ellipse per group. A group with fewer than three units MUST NOT receive an ellipse, and the omission MUST be labelled.
- **FR-210 — Differential: labelled volcano:** The system MUST plot the zero-null family as a volcano with the declared q cutoff, label only features present in the source table, and classify each point by its family membership.
- **FR-211 — Differential: MA plot:** The system MUST plot average log2 abundance against log2 fold change for each contrast, with the sign of the change following the declared numerator group.
- **FR-212 — Differential: per-protein dot plots with P brackets:** The system MUST plot individual values per group for declared proteins with mean ± SD, and MUST draw P brackets whose labels follow FR-198 from the family P values.
- **FR-213 — Differential: clustered DEP heatmap:** The system MUST plot row-scaled values of the declared DEP set with hierarchical clustering by the declared distance, and MUST NOT include features outside that set.
- **FR-214 — Differential: UpSet and Venn:** The system MUST draw an UpSet plot for any number of sets and a Venn only for three or fewer sets, with region counts equal to the exclusive-region partition of [R14a](../015-post-de-analysis/spec.md).
- **FR-215 — Biomarker: ROC with CI band:** The system MUST plot the ROC curve of the declared claim label with its bootstrap confidence band. The claim label MUST appear in the figure and its source table.
- **FR-216 — Biomarker: panel forest plot:** The system MUST plot panel coefficients or effects with their confidence intervals, with the claim label of the panel evaluation.
- **FR-217 — Biomarker: CV AUC against permutation null:** The system MUST plot the observed cross-validated AUC against the permutation null distribution and show the permutation P as (k + 1) / (B + 1) with B shown.
- **FR-218 — Pathways: enrichment dot plot:** The system MUST plot enrichment sets with their gene ratio and −log10 q for eligible sets with at least one foreground hit, and MUST NOT show sets without foreground hits as enriched.
- **FR-219 — Networks: module heatmap:** The system MUST plot module eigengenes per group or phenotype level for the declared modules only.
- **FR-220 — Networks: network graph:** The system MUST plot the declared subnetwork's nodes and edges exactly as in the interaction snapshot. Layout coordinates are presentation only and are excluded from data equality.
- **FR-221 — PERMANOVA: PCoA with centroids:** The system MUST plot principal coordinates of the declared distance with group centroids computed from each group's members, and MUST NOT draw a centroid from the full sample.
- **FR-222 — PERMANOVA: PERMDISP plot:** The system MUST plot distances to group centroids per group with the PERMDISP permutation P, showing the number of permutations and the seed used.

## Acceptance Scenarios

<a id="V195"></a>

### V195: `report.figures` configuration

**Fixture:** (a) A configuration without `report.figures`. (b) A configuration with explicit defaults. (c) A configuration with `dpi: 150`. (d) A configuration with `formats: ["jpg"]`. (e) A configuration with an unknown key.

**Oracle:** The default table in the contract, written in the test.

**Exact assertion:** (a) and (b) resolve to formats png, pzfx and html, dpi 300, width 120 mm, palette Okabe–Ito, font Arial and exact P annotation, and the resolved values are written to the run settings.

**Negative case:** (c) fails E_FIGURE_DPI. (d) fails E_FIGURE_FORMAT. (e) fails with the unknown-key error. No figure is produced in any negative case.

**Contract:** SM25; **owner:** R16a.

<a id="V196"></a>

### V196: Palette and stable group colours

**Fixture:** Three groups in declared order, and three different figure types generated in two different orders.

**Oracle:** The Okabe–Ito hex values written in the test and the declared group order.

**Exact assertion:** Group colours equal the oracle in every figure, and the colour of each group is identical across the two figure orders.

**Negative case:** A colour assignment that depends on figure order fails E_FIGURE_COLOUR_UNSTABLE.

**Contract:** SM25; **owner:** R16a.

<a id="V197"></a>

### V197: Figure registry

**Fixture:** A run with three registered catalogue figures and one request whose source table is absent.

**Oracle:** The registered figure IDs from the contract and the SHA-256 of each source table, computed in the test.

**Exact assertion:** Every registered figure lists its source path and the SHA-256 that matches the oracle, and every output is generated from its registered source.

**Negative case:** The request without a source table fails E_FIGURE_SOURCE_MISSING and produces no output.

**Contract:** SM25; **owner:** R16a.

<a id="V198"></a>

### V198: P-value annotation

**Fixture:** P values 0.2, 0.04, 0.0004 and 1e-7 in one table.

**Oracle:** The annotation rule table in the contract, written in the test.

**Exact assertion:** Stars mode gives ns, *, ***, and **** respectively under the contract thresholds. Exact mode gives p = 0.20, p = 0.040, p = 0.00040 and p < 0.0001 under the contract rounding (two significant digits with trailing zeros kept, never rounded onto or across a threshold, floor at 0.0001; review fix d2b8414).

**Negative case:** Stars shown when the declared mode is exact fails.

**Contract:** SM25; **owner:** R16a.

<a id="V199"></a>

### V199: Column tables for group comparisons

**Fixture:** Two groups of six replicate values for three proteins.

**Oracle:** The values read from the source table and parsed from the written `.pzfx` XML with the standard library XML parser, written in the test.

**Exact assertion:** Each column holds the source values exactly, with round-trip numeric precision (D-12), and the table type is Column. When the CRAN pzfx package is installed, `pzfx::read_pzfx` reads the same values back; otherwise that read-back is NOT_RUN. Opening the file in GraphPad Prism is a manual Maintainer check, NOT_RUN in automated tests.

**Negative case:** A single changed value, or rounding to three decimals, fails the comparison.

**Contract:** SM25; **owner:** R16b.

<a id="V200"></a>

### V200: XY tables

**Fixture:** A volcano source of 200 features and a ROC source of 30 points.

**Oracle:** Source rows, and the XML parsed independently in the test.

**Exact assertion:** Every X and Y value equals its source value, one point per source row, and the table type is XY. The `pzfx::read_pzfx` read-back is performed when pzfx is installed, otherwise NOT_RUN.

**Negative case:** Writing the values rounded to three decimals fails.

**Contract:** SM25; **owner:** R16b.

<a id="V201"></a>

### V201: Grouped tables for multi-protein panels

**Fixture:** Four proteins, two groups and six replicates per group.

**Oracle:** The source table and the independently parsed XML, written in the test.

**Exact assertion:** Rows are replicates, data sets are proteins by group, and every value equals the source. When the CRAN pzfx package is installed, `pzfx::read_pzfx` reads the same values back. When it is not installed, that read-back is NOT_RUN.

**Negative case:** A missing protein data set fails the comparison.

**Contract:** SM25; **owner:** R16b.

<a id="V202"></a>

### V202: Prism-style theme

**Fixture:** A grouped dot plot built with the Prism-style theme.

**Oracle:** The theme values declared in the style contract, written in the test.

**Exact assertion:** Background is white. Major and minor gridlines are blank. Axis lines have the declared thickness. Ticks point outward with the declared length. Text uses Arial, or a recorded substitution.

**Negative case:** A theme with gridlines fails E_FIGURE_STYLE.

**Contract:** SM25; **owner:** R16c.

<a id="V203"></a>

### V203: PNG dimensions and data fidelity

**Fixture:** One figure at 85 mm, 120 mm and 180 mm.

**Oracle:** Pixel width = width in inches × 300, rounded down (1003 px at 85 mm, 1417 px at 120 mm and 2125 px at 180 mm), and `ggplot_build` layer data taken from the source table.

**Exact assertion:** The PNG header dimensions equal the oracle, and each plotted layer equals its source rows.

**Negative case:** An export at 72 dpi fails E_FIGURE_DPI.

**Contract:** SM25; **owner:** R16c.

<a id="V204"></a>

### V204: Embedded trace data and hover fields

**Fixture:** A volcano source with group, feature identifier, effect, P and q.

**Oracle:** The source table and the trace arrays parsed from the HTML JSON, written in the test.

**Exact assertion:** The x, y and text arrays equal the source. The hover template includes the group, feature identifier and P or q fields.

**Negative case:** A hover template without the feature identifier fails E_PLOTLY_HOVER.

**Contract:** SM25; **owner:** R16d.

<a id="V205"></a>

### V205: Offline HTML and pinned plotly.js

**Fixture:** The HTML output of a figure, and the vendored plotly.js with its licence notice and pin record.

**Oracle:** SHA-256 of the vendored file computed in the test, compared to the pin record. A static scan with the rules in the contract.

**Exact assertion:** The hash equals the pin. The HTML has no `script src`, `link href`, `fetch`, `XMLHttpRequest` or dynamic `import` to a URL, and every `http(s)://` string is in the namespace allowlist. The MIT notice is present.

**Negative case:** A CDN script tag fails E_PLOTLY_NETWORK. A modified vendored file fails E_PLOTLY_HASH.

**Contract:** SM25; **owner:** R16d.

<a id="V206"></a>

### V206: QC intensity distributions

**Fixture:** Twelve samples in which two have a planted shift, with missing cells.

**Oracle:** Per-sample quantiles of observed cells, computed in the test.

**Exact assertion:** All three outputs exist. The quantiles in each output equal the oracle, and the observed mask matches.

**Negative case:** A source with imputed cells and no observed mask fails E_FIGURE_SOURCE_IMPUTED.

**Contract:** SM25, SM06; **owner:** R16e.

<a id="V207"></a>

### V207: QC missingness heatmap

**Fixture:** A 30 × 12 matrix with a planted missing block.

**Oracle:** Observed fraction per feature and per sample, computed in the test.

**Exact assertion:** All three outputs exist. Each cell equals the oracle and missing cells are rendered as missing.

**Negative case:** Missing cells drawn as zero fail.

**Contract:** SM25, SM06; **owner:** R16e.

<a id="V208"></a>

### V208: QC sample correlation heatmap

**Fixture:** Ten samples in two blocks, with missing cells.

**Oracle:** Pairwise correlations on observed cells with the SM06 coverage rule, computed in the test.

**Exact assertion:** All three outputs exist and every cell equals the oracle within 1e-10.

**Negative case:** A correlation computed on an imputed matrix fails E_FIGURE_SOURCE_IMPUTED.

**Contract:** SM25, SM06; **owner:** R16e.

<a id="V209"></a>

### V209: QC PCA with ellipses

**Fixture:** Two groups of eight samples, and a variant with a group of two.

**Oracle:** PCA scores sign-aligned to the reference, and 95 % ellipses from the group covariance, computed in the test.

**Exact assertion:** All three outputs exist. Scores agree with the oracle within 1e-10. The group of two has no ellipse and the omission is labelled.

**Negative case:** An ellipse drawn for the group of two without a label fails.

**Contract:** SM25; **owner:** R16e.

<a id="V210"></a>

### V210: Labelled volcano

**Fixture:** 200 features with planted up and down effects, and a declared q cutoff.

**Oracle:** Point class labels from the family q values, and the planted feature list, computed in the test.

**Exact assertion:** All three outputs exist. Point classes equal the oracle. Labelled features are exactly the declared subset present in the source.

**Negative case:** A label for a feature absent from the source table fails E_FIGURE_SOURCE_MISSING.

**Contract:** SM12, SM25; **owner:** R16e.

<a id="V211"></a>

### V211: MA plot

**Fixture:** Two contrasts over 200 features with a known fold change.

**Oracle:** A = (mean of numerator + mean of denominator) / 2 and M = mean of numerator − mean of denominator, computed in the test.

**Exact assertion:** All three outputs exist and A and M equal the oracle within 1e-10. The sign of M follows the declared numerator.

**Negative case:** A reversed sign of M fails.

**Contract:** SM25; **owner:** R16e.

<a id="V212"></a>

### V212: Per-protein dot plots with P brackets

**Fixture:** Three proteins in three groups of five replicates.

**Oracle:** Group means, SD and the family P values with the FR-198 label rule, computed in the test.

**Exact assertion:** All three outputs exist. Means, SD and every individual value equal the oracle, and each bracket label equals its FR-198 annotation.

**Negative case:** A bracket showing a P value that is not in the family table fails.

**Contract:** SM25, SM12; **owner:** R16e.

<a id="V213"></a>

### V213: Clustered DEP heatmap

**Fixture:** Forty planted DEPs in two clusters over twelve samples, plus 40 non-DEP features.

**Oracle:** Row order from `stats::hclust` with the declared distance and linkage, computed in the test.

**Exact assertion:** All three outputs exist. Row order and row-scaled values equal the oracle. Only the declared DEP set appears.

**Negative case:** A non-DEP feature in the heatmap fails E_FIGURE_SET_MISMATCH.

**Contract:** SM25, SM12; **owner:** R16e.

<a id="V214"></a>

### V214: UpSet and Venn

**Fixture:** Three and five overlapping synthetic sets, as in the R14a fixtures.

**Oracle:** Exclusive region counts by brute-force enumeration of the 2^k − 1 regions, written in the test.

**Exact assertion:** The UpSet figure has all three outputs for k = 3 and k = 5. The Venn variant is produced for k = 3 only, and its region counts equal the oracle and sum to the union.

**Negative case:** A Venn requested for k = 5 fails E_VENN_K and the UpSet still runs.

**Contract:** SM31, SM25; **owner:** R16e.

<a id="V215"></a>

### V215: ROC with CI band

**Fixture:** Two classes of 15 units with planted separation and out-of-fold scores.

**Oracle:** ROC points and AUC from the scores, and the stratified bootstrap band with the recorded seed, computed in the test.

**Exact assertion:** All three outputs exist. Curve points and AUC equal the oracle, the band equals the oracle, and the claim label is printed in the figure and in its source.

**Negative case:** A ROC from training-fold scores without the cross-validated label fails E_BIOMARKER_CLAIM.

**Contract:** SM34, SM25; **owner:** R16e.

<a id="V216"></a>

### V216: Panel forest plot

**Fixture:** A panel of five features with coefficients and intervals.

**Oracle:** Estimates and intervals from the model, computed in the test.

**Exact assertion:** All three outputs exist, and every estimate and interval equals the oracle. The claim label is shown.

**Negative case:** A panel plot without a claim label fails E_BIOMARKER_CLAIM.

**Contract:** SM34, SM25; **owner:** R16e.

<a id="V217"></a>

### V217: CV AUC against permutation null

**Fixture:** Nested-CV AUC repeats and a permutation null with B = 19 values from the V156 design.

**Oracle:** Permutation P = (k + 1) / (B + 1) computed in the test.

**Exact assertion:** All three outputs exist. The observed AUC and every null value equal the source. The displayed P equals the oracle and B is displayed.

**Negative case:** A displayed P that differs from the oracle fails.

**Contract:** SM36, SM25; **owner:** R16e.

<a id="V218"></a>

### V218: Enrichment dot plot

**Fixture:** The V176 ORA rows, including one set with no foreground hit.

**Oracle:** Gene ratio = k / n and −log10 q computed in the test.

**Exact assertion:** All three outputs exist. Gene ratio and −log10 q equal the oracle, and only eligible sets with at least one hit are plotted as enriched.

**Negative case:** A set with no foreground hit plotted as enriched fails.

**Contract:** SM18, SM25; **owner:** R16e.

<a id="V219"></a>

### V219: Module heatmap

**Fixture:** Three planted modules over two groups.

**Oracle:** Eigengenes as the first principal component of each scaled module, sign-aligned, computed in the test.

**Exact assertion:** All three outputs exist and each cell equals the oracle within 1e-10. Only declared modules appear.

**Negative case:** An undeclared module shown fails E_FIGURE_SET_MISMATCH.

**Contract:** SM39, SM25; **owner:** R16e.

<a id="V220"></a>

### V220: Network graph

**Fixture:** A 30-node interaction snapshot with a planted cluster and a declared subset.

**Oracle:** The induced subnetwork edge list from the snapshot, computed in the test.

**Exact assertion:** All three outputs exist. Nodes and edges equal the oracle. Layout coordinates are excluded from the equality.

**Negative case:** An edge absent from the snapshot fails E_FIGURE_SOURCE_MISSING.

**Contract:** SM40, SM25; **owner:** R16e.

<a id="V221"></a>

### V221: PCoA with centroids

**Fixture:** The distance matrix of the PERMANOVA synthetic fixture.

**Oracle:** Principal coordinates from `stats::cmdscale` sign-aligned to the reference, and centroids as group means, computed in the test.

**Exact assertion:** All three outputs exist. Coordinates agree with the oracle within 1e-8, and each centroid equals the mean of its own group.

**Negative case:** A centroid computed from the full sample fails.

**Contract:** SM26, SM25; **owner:** R16e.

<a id="V222"></a>

### V222: PERMDISP plot

**Fixture:** The PERMANOVA dispersion fixture with two groups.

**Oracle:** Distances to group centroids and the PERMDISP permutation P with the recorded seed, computed in the test.

**Exact assertion:** All three outputs exist. Distances and P equal the oracle, and the number of permutations and seed are displayed.

**Negative case:** A P displayed without B or seed fails.

**Contract:** SM29, SM25; **owner:** R16e.
