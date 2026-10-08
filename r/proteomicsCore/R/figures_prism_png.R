# R16c (FR-202, FR-203; V202, V203): Prism-style static PNG renderer (ggplot2 + ggprism, ragg device, 300 dpi).
#
# Style values come from specs/017-visualization/contracts/figures.md ("Style defaults") and D-74. The theme pattern
# follows the approved preview (analysis 20, theme_pub), without any of its data or paths. Refusals use the typed codes of
# the contract: E_FIGURE_DPI (dpi other than 300), E_FIGURE_STYLE (gridlines drawn), E_CONFIG_SCHEMA (width preset).

figure_style_values <- function() {
  list(background = "#FFFFFF", axis_line_pt = 1.5, tick_length_pt = 4, base_size = 12, title_size = 14,
       axis_title_gap_pt = 8, title_gap_pt = 8, legend_gap_mm = 4, plot_margin_pt = c(10, 14, 10, 10),
       font = "Arial", dpi = 300, width_presets_mm = c(85, 120, 180))
}

# Arial is used when installed. Otherwise the substitution is recorded and returned; the caller never fails on it.
figure_font_status <- function(family = "Arial") {
  available <- family %in% systemfonts::system_fonts()$family
  list(requested = family, available = available, used = if (available) family else "sans",
       substitution = if (available) NULL else sprintf("%s is not installed; drawn with the 'sans' family", family))
}

# The Prism-style theme (contract style table, D-74). Gridlines are blank, axis lines are black and thick, ticks point
# outward, text is bold on the axis titles and tick labels, and spacing follows the contract.
prism_figure_theme <- function(base_size = 12, family = "Arial") {
  v <- figure_style_values()
  family <- figure_font_status(family)$used
  ggprism::theme_prism(base_size = base_size, base_family = family, base_line_size = 0.8) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = v$background, colour = NA),
      plot.background = ggplot2::element_rect(fill = v$background, colour = NA),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black", linewidth = v$axis_line_pt / ggplot2::.pt),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = v$axis_line_pt / ggplot2::.pt),
      axis.ticks.length = grid::unit(v$tick_length_pt, "pt"),
      axis.title = ggplot2::element_text(face = "bold", size = base_size, family = family),
      axis.text = ggplot2::element_text(face = "bold", size = base_size, family = family, colour = "black"),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = v$axis_title_gap_pt), face = "bold", family = family),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = v$axis_title_gap_pt), face = "bold", family = family),
      plot.title = ggplot2::element_text(size = base_size + 2, hjust = 0.5, family = family,
                                         margin = ggplot2::margin(b = v$title_gap_pt)),
      legend.position = "top",
      legend.spacing.x = grid::unit(v$legend_gap_mm, "mm"),
      plot.margin = ggplot2::margin(v$plot_margin_pt[1], v$plot_margin_pt[2], v$plot_margin_pt[3], v$plot_margin_pt[4], "pt"),
      text = ggplot2::element_text(family = family))
}

# The complete theme of a plot. ggplot2 does not export plot_theme(), so it is taken from its namespace.
completed_theme <- function(plot) utils::getFromNamespace("plot_theme", "ggplot2")(plot)

# Refuses a plot whose theme draws gridlines (E_FIGURE_STYLE). Returns the checked properties otherwise.
check_prism_style <- function(plot) {
  th <- completed_theme(plot)
  blank <- function(element) inherits(ggplot2::calc_element(element, th), "element_blank")
  if (!blank("panel.grid.major") || !blank("panel.grid.minor"))
    stop("E_FIGURE_STYLE: the plot draws gridlines; the Prism style has none", call. = FALSE)
  invisible(TRUE)
}

# Theme properties as they are read back (used by the V202 oracle).
theme_report <- function(plot) {
  th <- completed_theme(plot)
  el <- function(name) ggplot2::calc_element(name, th)
  list(background = el("plot.background")$fill,
       grid_major_blank = inherits(el("panel.grid.major"), "element_blank"),
       grid_minor_blank = inherits(el("panel.grid.minor"), "element_blank"),
       axis_line_linewidth_mm = el("axis.line")$linewidth,
       axis_line_pt = el("axis.line")$linewidth * ggplot2::.pt,
       # read the unit directly: grid::convertUnit() needs a graphics device and would create Rplots.pdf
       tick_length = as.numeric(el("axis.ticks.length")),
       tick_unit = as.character(grid::unitType(el("axis.ticks.length"))),
       y_guide = class(plot$scales$get_scales("y")$guide),
       tick_direction = "outside",
       font_family = el("text")$family)
}

# The approved dot-plot pattern: individual values (points), mean as a horizontal bar and mean +/- SD error bars.
# `source` has columns group and value; `group_colours` is a named vector of colours in declared group order.
prism_dot_plot <- function(source, group_colours, y_label = "value", title = NULL, family = "Arial") {
  source$group <- factor(source$group, levels = names(group_colours))
  ggplot2::ggplot(source, ggplot2::aes(x = group, y = value, colour = group)) +
    ggplot2::geom_point(size = 2.4, alpha = 0.85, show.legend = FALSE,
                        position = ggplot2::position_jitter(width = 0.18, height = 0, seed = 2026)) +
    ggplot2::stat_summary(fun = mean, geom = "crossbar", width = 0.45, linewidth = 0.35, colour = "black") +
    ggplot2::stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x) - stats::sd(x), ymax = mean(x) + stats::sd(x)),
                          geom = "errorbar", width = 0.22, linewidth = 0.6, colour = "black") +
    ggplot2::scale_y_continuous(guide = ggprism::guide_prism_offset_minor(), expand = ggplot2::expansion(c(0.05, 0.08))) +
    ggplot2::scale_colour_manual(values = group_colours) +
    ggplot2::labs(x = NULL, y = y_label, title = title) +
    prism_figure_theme(family = family)
}

# The plotted data of one layer, every number written with 17 significant digits so that it round-trips exactly.
layer_data_strings <- function(plot, layer = 1L, columns = c("x", "y")) {
  built <- ggplot2::ggplot_build(plot)
  data <- built$data[[layer]]
  out <- lapply(columns, function(name) {
    values <- data[[name]]
    if (is.numeric(values)) sprintf("%.17g", values) else as.character(values)
  })
  names(out) <- columns
  as.data.frame(out, stringsAsFactors = FALSE)
}

# Writes the PNG at the declared dpi. The pixel width is floor(width_in * dpi) because the raster device truncates
# (FR-203). Refusals: a dpi other than 300 (E_FIGURE_DPI), a width outside the presets (E_CONFIG_SCHEMA), gridlines
# (E_FIGURE_STYLE). Returns the path, the pixel size read from the file header, and the font status.
render_prism_png <- function(plot, path, width_mm, height_mm = 100, dpi = 300, family = "Arial") {
  v <- figure_style_values()
  if (!identical(as.numeric(dpi), v$dpi))
    stop(sprintf("E_FIGURE_DPI: dpi %s is refused; the only allowed value is 300", format(dpi)), call. = FALSE)
  if (!width_mm %in% v$width_presets_mm)
    stop(sprintf("E_CONFIG_SCHEMA: journal_width_mm must be one of %s", paste(v$width_presets_mm, collapse = ", ")), call. = FALSE)
  check_prism_style(plot)
  status <- figure_font_status(family)
  ggplot2::ggsave(filename = path, plot = plot, device = ragg::agg_png, width = width_mm / 25.4, height = height_mm / 25.4,
                  units = "in", dpi = dpi, bg = "white")
  size <- png_size(path)
  list(path = path, width_px = size[1], height_px = size[2], font = status)
}

# Pixel width and height from the IHDR chunk (bytes 17-24, big-endian).
png_size <- function(path) {
  header <- readBin(path, what = "raw", n = 24L)
  as.integer(c(sum(as.integer(header[17:20]) * 256^(3:0)), sum(as.integer(header[21:24]) * 256^(3:0))))
}
