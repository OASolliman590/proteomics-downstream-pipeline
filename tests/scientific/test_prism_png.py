"""R16c acceptance: V202 (Prism-style theme) and V203 (PNG dimensions and data fidelity), FR-202 and FR-203.

The renderer is r/proteomicsCore/R/figures_prism_png.R, run through the project's R runner. Oracles are written here:
- the style values of contracts/figures.md (white background, no gridlines, 1.5 pt axis lines, 4 pt outward ticks,
  Arial or a recorded substitution);
- the PNG pixel size, read from the IHDR chunk in this file, and floor(width_mm / 25.4 * 300) for the width;
- the source rows of the synthetic fixture tests/fixtures/figures/png/dotplot_source.tsv, read with the csv module.
If Rscript or ggplot2, ggprism or ragg is missing, the tests are skipped and reported as NOT_RUN.
"""
from __future__ import annotations

import csv
import json
import math
import shutil
import struct
from pathlib import Path

import pytest

from proteomics_pipeline.runtime import run_r_code

ROOT = Path(__file__).resolve().parents[2]
R_FILE = ROOT / "r" / "proteomicsCore" / "R" / "figures_prism_png.R"
SOURCE = ROOT / "tests" / "fixtures" / "figures" / "png" / "dotplot_source.tsv"
GROUPS = ["Control", "Acute", "Chronic"]           # declared group order (the first three Okabe-Ito colours)


def _r_available() -> bool:
    if shutil.which("Rscript") is None:
        return False
    probe = run_r_code('cat(all(vapply(c("ggplot2", "ggprism", "ragg", "systemfonts", "jsonlite"), requireNamespace, logical(1), quietly = TRUE)))')
    return probe.returncode == 0 and probe.stdout.strip().endswith("TRUE")


pytestmark = pytest.mark.skipif(not _r_available(), reason="NOT_RUN: Rscript with ggplot2, ggprism, ragg, systemfonts and jsonlite is unavailable")

DRIVER_RENDER = """
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(proteomicsCore))
out_dir <- args[1]
src <- read.delim(args[2], stringsAsFactors = FALSE)
src$value <- as.numeric(src$value)
cols <- c(Control = "#0072B2", Acute = "#E69F00", Chronic = "#D55E00")
p <- prism_dot_plot(src, cols, y_label = "log2 abundance", title = "Synthetic dot plot")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
sizes <- list()
for (w in c(85, 120, 180)) {
  r <- render_prism_png(p, file.path(out_dir, sprintf("dotplot_%d.png", w)), width_mm = w)
  sizes[[as.character(w)]] <- list(width_px = r$width_px, height_px = r$height_px, font = r$font$used)
}
layer <- layer_data_strings(p, 1L, c("x", "y"))
cat(jsonlite::toJSON(list(sizes = sizes, layer_x = layer$x, layer_y = layer$y), auto_unbox = TRUE, null = "null"))
"""

DRIVER_THEME = """
suppressPackageStartupMessages(library(proteomicsCore))
src <- data.frame(group = c("A", "B"), value = c(1, 2))
p <- prism_dot_plot(src, c(A = "#0072B2", B = "#E69F00"))
cat(jsonlite::toJSON(list(theme = theme_report(p), font = figure_font_status("Arial")), auto_unbox = TRUE, null = "null", digits = NA))
"""

DRIVER_GRIDLINES = """
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(proteomicsCore))
p <- ggplot2::ggplot(data.frame(g = c("A", "B"), v = c(1, 2)), ggplot2::aes(g, v)) + ggplot2::geom_point() + ggplot2::theme_grey()
render_prism_png(p, file.path(args[1], "grid.png"), width_mm = 120)
"""

DRIVER_DPI = """
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(proteomicsCore))
src <- data.frame(group = c("A", "B"), value = c(1, 2))
p <- prism_dot_plot(src, c(A = "#0072B2", B = "#E69F00"))
render_prism_png(p, file.path(args[1], "dpi72.png"), width_mm = 120, dpi = 72)
"""


def _png_size(path: Path) -> tuple[int, int]:
    """Width and height from the IHDR chunk: the 8-byte signature, the chunk length and type, then two big-endian ints."""
    data = path.read_bytes()
    assert data[:8] == b"\x89PNG\r\n\x1a\n", "not a PNG file"
    assert data[12:16] == b"IHDR"
    return struct.unpack(">II", data[16:24])


def _render(tmp_path: Path) -> dict:
    result = run_r_code(DRIVER_RENDER, [tmp_path / "png", SOURCE], cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout.strip().splitlines()[-1])


def _source_rows() -> list[dict]:
    with SOURCE.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


# ---------------------------------------------------------------- V202: Prism-style theme

def test_v202_theme_values_equal_the_style_contract(tmp_path):
    result = run_r_code(DRIVER_THEME, [], cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    report = json.loads(result.stdout.strip().splitlines()[-1])
    theme, font = report["theme"], report["font"]
    # Oracle: contracts/figures.md "Style defaults" (1.5 pt axis line in ggplot points, 4 pt outward ticks, white).
    assert theme["background"].upper() == "#FFFFFF"
    assert theme["grid_major_blank"] is True and theme["grid_minor_blank"] is True
    assert theme["axis_line_pt"] == pytest.approx(1.5, abs=1e-9)
    assert theme["axis_line_linewidth_mm"] == pytest.approx(1.5 / (72.27 / 25.4), abs=1e-9)
    assert theme["tick_length"] == pytest.approx(4.0, abs=1e-9)
    assert theme["tick_unit"] == "points"          # grid names the pt unit "points"
    assert "prism_offset_minor" in theme["y_guide"]   # ggprism offset axis with minor ticks
    assert theme["tick_direction"] == "outside"
    # Arial, or a recorded substitution when Arial is not installed (never a silent change).
    if font["available"]:
        assert theme["font_family"] == "Arial"
        assert font["substitution"] is None
    else:
        assert font["substitution"]
        assert theme["font_family"] == font["used"]


def test_v202_negative_a_theme_with_gridlines_fails_e_figure_style(tmp_path):
    result = run_r_code(DRIVER_GRIDLINES, [tmp_path], cwd=tmp_path)
    assert result.returncode != 0
    assert "E_FIGURE_STYLE" in result.stderr


# ---------------------------------------------------------------- V203: PNG dimensions and data fidelity

def _floor_pixels(width_mm: float, dpi: int = 300) -> int:
    return math.floor(width_mm / 25.4 * dpi)


def test_v203_png_header_dimensions_equal_the_oracle_at_three_widths(tmp_path):
    # Oracle: pixel width = floor(width in inches x 300): 1003, 1417 and 2125 px for 85, 120 and 180 mm.
    assert {w: _floor_pixels(w) for w in (85, 120, 180)} == {85: 1003, 120: 1417, 180: 2125}
    rendered = _render(tmp_path)
    for width in (85, 120, 180):
        width_px, height_px = _png_size(tmp_path / "png" / f"dotplot_{width}.png")
        assert width_px == _floor_pixels(width)
        assert rendered["sizes"][str(width)]["width_px"] == width_px
        assert height_px == _floor_pixels(100)                  # the default 100 mm height
        assert height_px == rendered["sizes"][str(width)]["height_px"]


def test_v203_plotted_layer_equals_the_source_rows(tmp_path):
    rendered = _render(tmp_path)
    rows = _source_rows()
    # Oracle: each source row is one point. y is the value, exactly. x is the group position plus a deterministic jitter
    # of at most 0.18 of a category width (contracts/figures.md, dot-plot jitter).
    positions = [float(GROUPS.index(r["group"]) + 1) for r in rows]
    expected_y = [float(r["value"]) for r in rows]
    got_x = [float(v) for v in rendered["layer_x"]]
    got_y = [float(v) for v in rendered["layer_y"]]
    assert got_y == expected_y
    assert len(got_x) == len(positions)
    assert all(abs(x - p) <= 0.18 + 1e-9 for x, p in zip(got_x, positions))
    assert any(abs(x - p) > 1e-6 for x, p in zip(got_x, positions)), "points must be jittered, not stacked"


def test_v203_negative_unjittered_points_stack_on_the_group_line(tmp_path):
    rendered = _render(tmp_path)
    rows = _source_rows()
    stacked = [float(GROUPS.index(r["group"]) + 1) for r in rows]
    # the stacked layout is what the fidelity rule refuses: the rendered x values differ from the group positions
    assert [float(v) for v in rendered["layer_x"]] != stacked


def test_v203_negative_export_at_72_dpi_fails_e_figure_dpi(tmp_path):
    result = run_r_code(DRIVER_DPI, [tmp_path], cwd=tmp_path)
    assert result.returncode != 0
    assert "E_FIGURE_DPI" in result.stderr
    assert not (tmp_path / "dpi72.png").exists()
