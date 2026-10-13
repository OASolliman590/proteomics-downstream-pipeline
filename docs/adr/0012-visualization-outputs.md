# ADR 0012 — Visualization outputs scope amendment

Status: scope decided by the Maintainer on 2026-10-08. Recorded by the Claude route acting under [ADR 0007](0007-claude-route-and-integration-amendments.md). The specification drafted under this ADR is pending the Maintainer's review; no packet is dispatched by this ADR. Decisions D-68 and D-72 in [decisions.md](../../specs/001-downstream-proteomics/decisions.md).

## Context

- The supervisor restyles figures in GraphPad Prism, so every figure needs a Prism-editable data file, not only a picture.
- Figures must be readable in print at journal widths and also explorable offline by a reviewer, without a network connection.
- The existing report produces figures from saved figure-source tables ([R10b](../../specs/001-downstream-proteomics/packet-index.md), FR-097), in PDF, SVG and base-graphics PNG. Those outputs stay as they are. They do not provide the Prism file, the Prism-style 300 dpi PNG or the interactive HTML.

## Decision

- Add slice **017-visualization**, Phase 6 (v1.3-visualization), with dispatch units **R16a–R16e** and **FR-195–FR-222, T195–T222, V195–V222**, appended after FR-194. No identity is renumbered or retired.
- Every figure in the new catalogue has three outputs:
  1. a GraphPad Prism `.pzfx` file (column tables for group comparisons, XY tables for volcano and ROC data, grouped tables for multi-protein panels);
  2. a finished Prism-style static PNG at 300 dpi (white background, thick black axes, outward ticks, no gridlines, Arial, mean ± SD over individual points), produced with R ggplot2 and the ggprism package;
  3. a standalone interactive HTML file per figure, using Plotly. plotly.js is vendored offline, version-pinned, hash-checked and accompanied by its MIT licence notice. The file makes no network requests.
- Figure areas, all in the first build: QC and PCA; differential; biomarker and ROC; pathways and networks; PERMANOVA PCoA and PERMDISP.
- Every figure is generated from its saved figure-source table, and each output is checked against that table.
- Phase 6 may be built before Phase 5 because the two phases have no dependency edge and disjoint write allowlists. This ADR does not authorize concurrent dispatch.

## Consequences

- Shared-interface amendment before R16a, made by the Maintainer: an optional `report.figures` block in `analysis.schema.json` (formats defaulting to `["png", "pzfx", "html"]`, dpi 300, journal width presets, palette, font and p-value annotation). The existing `report.figure_formats` key remains valid and keeps its current meaning; the two are never merged silently (see [the contract](../../specs/017-visualization/contracts/figures.md)).
- The existing R10b, R13 and R14 figure producers are accepted packets with disjoint allowlists. They are not retrofitted by this ADR. Retrofitting them to the three-output rule would need a separate Maintainer amendment.
- R Suggests entries: ggplot2, ggprism, pzfx (optional, read-back oracle) and plotly-related vendored files. A missing optional package gives NOT_RUN, never a silent fallback to another output format.
- Automated tests cannot open a file in Prism. Prism-openability is a manual Maintainer check and is NOT_RUN in automated tests.

## Verification

V195–V222. Automated oracles: independent XML parsing of each `.pzfx`; PNG header dimensions equal width in inches × 300 (rounded); `ggplot_build` data equal to the source tables; embedded Plotly trace data equal to the source tables; a static scan of each HTML file for external loads; and the vendored plotly.js hash equal to the pinned value.
