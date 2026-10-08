# Figure configuration, output layout and figure IDs (R16a–R16e)

This contract fixes the configuration grammar, the three outputs of every figure, the output layout and the catalogue figure IDs. Figure values come only from saved figure-source tables. Unknown keys are refused. Defaults below are prespecified product policies, not claims about the data.

## Configuration: `report.figures` (optional block of analysis.schema.json)

| Key | Allowed values | Default | Meaning |
|---|---|---|---|
| `formats` | subset of `png`, `pzfx`, `html`, non-empty | `["png", "pzfx", "html"]` | Outputs produced for each catalogue figure. |
| `dpi` | `300` | `300` | PNG resolution. Any other value fails E_FIGURE_DPI. |
| `journal_width_mm` | `85`, `180` | `85` | Figure width preset. |
| `palette` | `okabe_ito`, or an explicit list of hex colours with one entry per group | `okabe_ito` | Group palette. |
| `group_colours` | optional map from group identifier to hex colour | none | Overrides the palette. Stable across all figures of a run. |
| `font` | `Arial` | `Arial` | If unavailable, the substitution is recorded in the registry. |
| `p_annotation` | `stars`, `exact` | `exact` | P-value annotation style (FR-198). |

### Relation to `report.figure_formats`

`report.figure_formats` (R10b) keeps its meaning and drives the existing PDF, SVG and PNG figures exactly as before. `report.figures` drives only the catalogue in this contract. If both are declared, both run, with separate outputs and separate directories. The two keys are never merged silently. An invalid `report.figures` block makes only the catalogue INAPPLICABLE with its typed reason. The existing report is unaffected, in line with [ADR 0010](../../../docs/adr/0010-adaptive-post-de-policy.md).

## Style defaults

| Property | Default |
|---|---|
| Background | white |
| Gridlines | none (major and minor) |
| Axis lines | black, thick: 1.5 pt |
| Ticks | outward, length 4 pt |
| Base font | Arial, 9 pt |
| Points in dot plots | individual values, mean as a horizontal bar, error bar ± SD |
| Okabe–Ito palette | `#E69F00`, `#56B4E9`, `#009E73`, `#F0E442`, `#0072B2`, `#D55E00`, `#CC79A7`, `#000000` in that order |

Group order is the declared group order of the frozen plan. A `group_colours` map overrides it.

## P-value annotation rule (FR-198)

| P value | `stars` | `exact` |
|---|---|---|
| P ≥ 0.05 | `ns` | `p = ` two significant digits in plain decimal, for example `p = 0.20` |
| 0.01 ≤ P < 0.05 | `*` | same rule, for example `p = 0.04` |
| 0.001 ≤ P < 0.01 | `**` | same rule |
| 0.0001 ≤ P < 0.001 | `***` | same rule, for example `p = 0.0004` |
| P < 0.0001 | `****` | `p < 0.0001` |

The annotation is computed from the P value in the source table. The rule applies to the family P value named in the figure (for example the zero-null family).

## Output layout

```text
<run>/report/figure_catalogue/
  <figure_id>.png          Prism-style PNG at the declared dpi (R16c)
  <figure_id>.pzfx         GraphPad Prism file (R16b)
  <figure_id>.html        standalone interactive HTML (R16d)
  sources/<figure_id>.tsv  figure-source table, the only input to all three outputs
  registry.json            per figure: state, reason, source path and SHA-256, output paths and SHA-256, font substitution if any
```

`<run>/report/figures/` and `<run>/report/figure_sources/` (R10b) are not written by the catalogue. Each output has a SHA-256 in `registry.json`. A missing output is recorded as NOT_RUN with its reason (E_FIGURE_OUTPUT_MISSING). The other outputs of that figure still run.

## Vendored plotly.js (R16d, FR-205)

```text
vendor/plotly/
  plotly.min.js        vendored file, version pinned
  LICENSE              the MIT licence notice of plotly.js, unmodified
  PINNED.json          {package, version, sha256, licence, vendored_by_packet}
```

`PINNED.json` is written when the Maintainer approves the vendored version, and its SHA-256 is recorded then. Until approval, R16d is NOT_RUN. The HTML embeds the vendored file inline after the hash check. It has no external resource reference. The static scan allows only the XML namespace strings `http://www.w3.org/2000/svg`, `http://www.w3.org/1999/xhtml` and `http://www.w3.org/1999/xlink`, and rejects every other URL, `<script src>`, `<link href>`, `fetch(`, `XMLHttpRequest` and dynamic `import(` to a URL (E_PLOTLY_NETWORK).

## Prism `.pzfx` table rules (R16b)

| Figure kind | Prism table | Rows and columns |
|---|---|---|
| Group comparisons (dot plots, per-protein values) | Column | One column per group. Each row is a replicate value, exactly as in the source. |
| Volcano, ROC, XY scatter | XY | One point per source row. X and Y are exactly the source values. |
| Multi-protein panels | Grouped | One row per replicate. One data set per protein and group. |

The writer self-check compares every value it writes with the source table (E_PZFX_VALUE_MISMATCH). Values are written with round-trip numeric precision (D-12). Opening the file in GraphPad Prism is a Maintainer check, recorded as NOT_RUN in automated tests.

## Figure IDs

| Area | Figure ID | Packet | Source family or analysis | Claim or method label |
|---|---|---|---|---|
| QC | `qc_intensity_distributions` | R16e | observed log2 intensity (SM06) | descriptive |
| QC | `qc_missingness_heatmap` | R16e | observed mask (SM04) | descriptive |
| QC | `qc_sample_correlation_heatmap` | R16e | SM06 coverage rule | descriptive |
| QC | `qc_pca_ellipses` | R16e | PCA on observed cells | descriptive |
| Differential | `diff_volcano_labelled` | R16e | zero-null family (SM12) | family q |
| Differential | `diff_ma_plot` | R16e | contrast means | descriptive |
| Differential | `diff_protein_dotplot_brackets` | R16e | family P (FR-198) | family P |
| Differential | `diff_dep_heatmap_clustered` | R16e | declared DEP set | descriptive |
| Differential | `diff_upset_venn` | R16e | SM31 exclusive regions (Venn for k ≤ 3) | descriptive |
| Biomarker | `bio_roc_ci_band` | R16e | SM34 | claim label |
| Biomarker | `bio_panel_forest` | R16e | SM34 | claim label |
| Biomarker | `bio_cv_auc_vs_permutation` | R16e | SM36 | claim label |
| Pathways | `path_enrichment_dotplot` | R16e | SM18 | family q |
| Networks | `net_module_heatmap` | R16e | SM39 | module_level |
| Networks | `net_network_graph` | R16e | SM40 | descriptive |
| PERMANOVA | `perm_pcoa_centroids` | R16e | SM26 | descriptive |
| PERMANOVA | `perm_permdisp` | R16e | SM29 | family P |

Seventeen figure IDs. Each has one V case (V206–V222) that checks all three outputs against the same source table.

## Error codes introduced by this slice

| Code | Packet | Behaviour |
|---|---|---|
| `E_FIGURE_FORMAT` | R16a | Unknown format. No figure is produced. |
| `E_FIGURE_DPI` | R16a, R16c | A dpi other than 300 for PNG. |
| `E_FIGURE_PALETTE` | R16a | Invalid palette or group colour list. |
| `E_FIGURE_COLOUR_UNSTABLE` | R16a | Group colours depend on figure order. |
| `E_FIGURE_SOURCE_MISSING` | R16a | A figure has no registered source table, or a feature is absent from its source. |
| `E_FIGURE_SOURCE_IMPUTED` | R16e | A source includes imputed cells without the observed mask. |
| `E_FIGURE_SET_MISMATCH` | R16e | A figure shows a feature, set or module outside its declared set. |
| `E_FIGURE_STYLE` | R16c | A theme property differs from the style defaults. |
| `E_FIGURE_OUTPUT_MISSING` | R16a–R16e | One output is not produced. Recorded as NOT_RUN; the other outputs run. |
| `E_PZFX_VALUE_MISMATCH` | R16b | The writer self-check finds a value that differs from the source. |
| `E_PLOTLY_HOVER` | R16d | Hover text lacks an identifying field. |
| `E_PLOTLY_NETWORK` | R16d | An external resource reference or network call in the HTML. |
| `E_PLOTLY_HASH` | R16d | The vendored plotly.js hash differs from the pin. |

Existing codes are reused: `E_VENN_K` (R14a) and `E_BIOMARKER_CLAIM` (R14d).

## Dependencies

R Suggests: ggplot2 and ggprism (R16c), pzfx (optional read-back, R16b; a missing package is NOT_RUN). The Python writers need only the standard library for `.pzfx` and HTML. A missing dependency produces NOT_RUN for that output, never a silent change of format.
