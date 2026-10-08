# R16c Claude-route receipt: static Prism-style PNG (ggplot2 + ggprism, 300 dpi, ADR 0012)

Packet R16c (FR-202, FR-203 / T202, T203 / V202, V203) was implemented on 2026-10-08 by the Claude route (Haiku builder, Opus review), at the operator's direction. This is not work by the routes named in AGENTS.md. The receipt means ready for review, not accepted. Independent audit: **NOT_RUN** (the Opus review follows this receipt).

## Identities

- Contract: `specs/017-visualization/contracts/figures.md` ("Style defaults", the Dependencies section, the error codes). Decisions: D-68, D-74 (120 mm, Okabe–Ito, rounding-down PNG width). Style pattern: the approved preview `analysis/20_FIGURE_STYLE_PREVIEW/make_preview.R`, read for its theme code only; none of its data, paths or outputs were used.
- Working-source manifest: `evidence/working-source-manifest.txt`, 9 files, SHA-256 `67c41f35bafa50ac2e7624e385290834058f590a3318130cc73dcd3e5239401a` (after the review fix). It covers the R renderer, its testthat file, the Python acceptance test, the synthetic fixture, the package NAMESPACE and the release scanner (A-2026-10-01-34, -35), and the shared files of A-2026-10-01-32 (DESCRIPTION, renv.lock, reason-code catalogue).
- Shared-interface amendment made first: **A-2026-10-01-32** (`packet-ownership.json`): ggplot2, ggprism, ragg and systemfonts in R Suggests; ggprism 1.0.7, ragg 1.5.2 and digest 0.6.39 in `renv.lock`; the reason-code catalogue (`E_FIGURE_STYLE`); traceability; task ticks T202–T203.
- `verified_commit`: **null**. The R16c commit is local on `claude/visualization` and is not pushed.
- Environment: macOS, Python 3.13.15 (`.venv`), R 4.6.1, ggplot2 4.0.3, ggprism 1.0.7, ragg 1.5.2, systemfonts 1.3.2, `LANG=en_US.UTF-8`. Arial is installed on this machine, so the substitution path is not exercised here; the tests accept both outcomes and record the substitution when it occurs.

## Implementation

- `r/proteomicsCore/R/figures_prism_png.R`
  - `prism_figure_theme()`: ggprism `theme_prism` plus the contract values (white background, no gridlines, 1.5 pt black axis lines, 4 pt outward ticks, bold axis titles and tick labels, 8 pt title and axis-title gaps, 4 mm legend spacing, margins of 10, 14, 10 and 10 pt).
  - `figure_font_status()`: records the Arial substitution instead of failing.
  - `prism_dot_plot()`: the approved dot-plot pattern (points, mean bar, mean ± SD error bars) with stable group colours passed in.
  - `render_prism_png()`: `ragg::agg_png` at 300 dpi, a white background, the pixel size read from the IHDR. Refusals: `E_FIGURE_DPI` (dpi other than 300), `E_FIGURE_STYLE` (gridlines), `E_CONFIG_SCHEMA` (width preset).
  - `layer_data_strings()`: the ggplot_build layer data as 17-significant-digit strings.
- The pixel width is `floor(width_mm / 25.4 * 300)`. ragg truncates, so 85, 120 and 180 mm give 1003, 1417 and 2125 px, and the height of 100 mm gives 1181 px.

## Gates executed

| Command | Exit | Result | Log |
|---|---|---|---|
| `.venv/bin/python -m pytest tests/scientific/test_prism_png.py -q -rs -p no:cacheprovider` | 0 | 6 passed | `evidence/python-tests-prism-png.log` |
| `Rscript --vanilla -e "testthat::test_file('r/proteomicsCore/tests/testthat/test-figures-prism-png.R', …)"` | 0 | 19 expectations, 0 failed, 0 skipped | `evidence/r-test-figures-prism-png.log` |
| `R CMD INSTALL --library="$R_LIBS_USER" r/proteomicsCore` | 0 | installed (DESCRIPTION changed) | `evidence/r-install.log` |
| `.venv/bin/python scripts/check_spec_kit.py` | 0 | PASS | `evidence/spec-kit-check.log` |
| `.venv/bin/python -m pytest tests/regression tests/unit -q -p no:cacheprovider` | 0 | 85 passed | `evidence/regression-unit.log` |

The renv lock check (`reproduction.check_r_lock('renv.lock')`) is qualified, with no mismatches (`evidence/renv-lock-check.log`).

## Acceptance cases

| Case | Status | Positive oracle | Negative case |
|---|---|---|---|
| V202 | PASS | the style values written in the test: white background, gridlines blank, axis line 1.5 pt (linewidth 1.5 / (72.27/25.4) mm), ticks 4 pt and outward, font Arial or a recorded substitution. Read back from the ggplot2 theme in R. | a theme with gridlines (`theme_grey`) fails `E_FIGURE_STYLE` (Python and testthat). |
| V203 | PASS | the PNG IHDR width equals `floor(width_mm / 25.4 * 300)` (1003, 1417, 2125 px), the height equals `floor(100 / 25.4 * 300)` = 1181 px, and the ggplot_build layer data equal the 15 source rows (y exactly, x within 0.18 of the group position, jittered) of `tests/fixtures/figures/png/dotplot_source.tsv`. | an export at 72 dpi fails `E_FIGURE_DPI` and writes no file. |

## Notes and open points

- **Prism-style check depends on the theme.** `check_prism_style()` reads the complete theme through `utils::getFromNamespace("plot_theme", "ggplot2")`, because ggplot2 4.0 does not export `plot_theme()`. If that internal changes, the gridline check must be updated.
- **Axis-line unit.** The contract says "1.5 pt". The renderer sets `linewidth = 1.5 / ggplot2::.pt` mm, and the oracle uses the same conversion. The Maintainer should confirm the unit on the approved preview.
- **Arial substitution is not exercised on this machine.** On Linux or Windows CI, the test accepts `substitution` with the theme family set to the substitute, and the status is returned to the caller for the registry. The registry write belongs to R16e.
- No private data was used. The fixture is synthetic.

## Review fix (commit after b7ae9ea; the R16c commit is 2c55bbb)

A reviewer rendered the dot plot and compared it with the approved preview. Three fixes, each with tests:

1. **Jitter.** The points were stacked on the group line. They are now jittered with `position_jitter(width = 0.18, height = 0, seed = 2026)`, so the layout is deterministic. V203 compares y exactly and requires |x - group position| <= 0.18. A negative case refuses the stacked layout.
2. **Axes.** The y scale uses `ggprism::guide_prism_offset_minor()` (the offset axis with minor ticks, as in the preview). The theme keeps the contract's 1.5 pt axis line, and `base_line_size = 0.8` as in the preview. The V202 test asserts the guide class. The x axis is discrete, so the offset guide is on the y axis only; the preview does the same. A visual check of the rendered PNG was made.
3. **Exports and tests.** The public functions (`prism_figure_theme`, `prism_dot_plot`, `render_prism_png`, `png_size`, `figure_font_status`, `check_prism_style`, `layer_data_strings`, `theme_report`) are exported in `NAMESPACE` (A-2026-10-01-34). The tests call the installed package, not the source file.

Other changes:
- **Stray `Rplots.pdf`.** `grid::convertUnit()` in `theme_report` needs a graphics device, so R created `Rplots.pdf`. The unit is now read directly (`as.numeric` and `grid::unitType`). The Python tests run R in a temporary directory, and the testthat run is executed from a temporary directory too. No file is left in the repository.
- **Release scanner (A-2026-10-01-35).** The release secret rule flagged a JavaScript expression in the vendored plotly.js (`access_token="+(...`). The scanner now exempts only `vendor/plotly/plotly.min.js` while its SHA-256 equals the pin. This narrows a security control and needs Maintainer review. A modified copy is still scanned (checked).

## Acceptance
All gates exited 0, and every case, including its negative case, passed. R16c is recorded as ready for review (PASS by the Claude route; independent audit NOT_RUN).
