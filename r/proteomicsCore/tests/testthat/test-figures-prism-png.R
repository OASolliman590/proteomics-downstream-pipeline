# R16c at the R level (V202, V203; FR-202, FR-203). Oracles are written here: the style values of
# contracts/figures.md, floor(width_mm / 25.4 * 300) for the pixel width, and the synthetic source rows of
# tests/fixtures/figures/png/dotplot_source.tsv. The Python acceptance test is tests/scientific/test_prism_png.py.

ns <- asNamespace("proteomicsCore")
fn <- function(name) get(name, envir = ns)
fixture <- normalizePath(file.path(testthat::test_path(), "..", "..", "..", "..", "tests", "fixtures", "figures", "png", "dotplot_source.tsv"), mustWork = TRUE)

testthat::test_that("V202 the Prism-style theme has white background, no gridlines, 1.5 pt axis lines and 4 pt outward ticks", {
  testthat::skip_if_not_installed("ggprism")
  src <- utils::read.delim(fixture, stringsAsFactors = FALSE)
  p <- fn("prism_dot_plot")(src, c(Control = "#0072B2", Acute = "#E69F00", Chronic = "#D55E00"))
  t <- fn("theme_report")(p)
  testthat::expect_identical(toupper(t$background), "#FFFFFF")
  testthat::expect_true(t$grid_major_blank && t$grid_minor_blank)
  testthat::expect_equal(t$axis_line_pt, 1.5, tolerance = 1e-9)
  testthat::expect_equal(t$tick_length, 4, tolerance = 1e-9)
  testthat::expect_identical(t$tick_direction, "outside")
  status <- fn("figure_font_status")("Arial")
  testthat::expect_true(identical(t$font_family, status$used))
})

testthat::test_that("V202 negative: a theme with gridlines fails E_FIGURE_STYLE", {
  testthat::skip_if_not_installed("ggprism")
  p <- ggplot2::ggplot(data.frame(g = c("A", "B"), v = c(1, 2)), ggplot2::aes(g, v)) + ggplot2::geom_point() + ggplot2::theme_grey()
  testthat::expect_error(fn("check_prism_style")(p), "E_FIGURE_STYLE")
})

testthat::test_that("V203 PNG pixel width is floor(width_mm / 25.4 * 300) at 85, 120 and 180 mm, and the layer equals the source", {
  testthat::skip_if_not_installed("ragg")
  src <- utils::read.delim(fixture, stringsAsFactors = FALSE)
  p <- fn("prism_dot_plot")(src, c(Control = "#0072B2", Acute = "#E69F00", Chronic = "#D55E00"))
  out <- tempfile("prism-png-")
  dir.create(out)
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  for (w in c(85, 120, 180)) {
    r <- fn("render_prism_png")(p, file.path(out, sprintf("d%d.png", w)), width_mm = w)
    testthat::expect_identical(r$width_px, as.integer(floor(w / 25.4 * 300)))
    testthat::expect_identical(r$height_px, as.integer(floor(100 / 25.4 * 300)))
  }
  layer <- fn("layer_data_strings")(p, 1L, c("x", "y"))
  testthat::expect_identical(as.numeric(layer$y), as.numeric(src$value))
  testthat::expect_identical(as.numeric(layer$x), as.numeric(match(src$group, c("Control", "Acute", "Chronic"))))
})

testthat::test_that("V203 negative: a dpi of 72 fails E_FIGURE_DPI and writes nothing", {
  testthat::skip_if_not_installed("ragg")
  p <- ggplot2::ggplot(data.frame(g = c("A", "B"), v = c(1, 2)), ggplot2::aes(g, v)) + ggplot2::geom_point()
  out <- tempfile("prism-dpi-")
  dir.create(out)
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  testthat::expect_error(fn("render_prism_png")(p, file.path(out, "x.png"), width_mm = 120, dpi = 72), "E_FIGURE_DPI")
  testthat::expect_identical(length(list.files(out)), 0L)
})
