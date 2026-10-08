# R16e (FR-206 to FR-214; V206-V214): static Prism-style PNG renderers for the QC and differential figures.
#
# Each renderer reads only its figure-source table (the registered copy written by the R16a registry), draws with the
# R16c Prism theme (prism_figure_theme, render_prism_png) and refuses with a typed code where the contract requires it:
# E_FIGURE_SOURCE_IMPUTED (no observed mask), E_FIGURE_SOURCE_MISSING (a label or bracket without its source),
# E_FIGURE_SET_MISMATCH (a feature outside the declared set), E_VENN_K (Venn for k > 3), E_FIGURE_STYLE (gridlines).
# Source-table columns are documented in tests/fixtures/figures/catalogue/make_fixtures.py.

catalogue_colours <- function() {
  c(up = "#D55E00", down = "#0072B2", ns = "grey70", Control = "#0072B2", Acute = "#E69F00", Chronic = "#D55E00")
}

catalogue_refuse <- function(code, message) stop(sprintf("%s: %s", code, message), call. = FALSE)

# Coverage-aware Pearson correlations (the rule of R03 qc_correlations): cells observed in both samples only, and at
# least `minimum_shared` shared cells, else NA.
catalogue_correlations <- function(values, observed, minimum_shared = 3L) {
  samples <- colnames(values)
  out <- matrix(NA_real_, length(samples), length(samples), dimnames = list(samples, samples))
  for (i in seq_along(samples)) for (j in seq_along(samples)) {
    shared <- observed[, i] & observed[, j]
    if (sum(shared) >= minimum_shared) out[i, j] <- stats::cor(values[shared, i], values[shared, j])
  }
  out
}

catalogue_wide <- function(source, row_col, col_col, value_col) {
  m <- tapply(source[[value_col]], list(source[[row_col]], source[[col_col]]), function(x) x[1])
  matrix(as.numeric(m), nrow(m), ncol(m), dimnames = dimnames(m))
}

# 95 per cent ellipse of a group from its covariance (chi-square with 2 df), as a closed polygon.
catalogue_ellipse <- function(x, y, level = 0.95, n = 100) {
  center <- c(mean(x), mean(y))
  cov_m <- stats::cov(cbind(x, y))
  radius <- sqrt(stats::qchisq(level, df = 2))
  eig <- eigen(cov_m, symmetric = TRUE)
  theta <- seq(0, 2 * pi, length.out = n)
  circle <- rbind(cos(theta), sin(theta))
  points <- center + radius * (eig$vectors %*% diag(sqrt(pmax(eig$values, 0))) %*% circle)
  data.frame(x = points[1, ], y = points[2, ])
}

catalogue_qc_intensity <- function(source) {
  if (!"observed" %in% names(source)) catalogue_refuse("E_FIGURE_SOURCE_IMPUTED", "the source has no observed mask; imputed cells cannot be plotted as observed")
  obs <- source[source$observed, ]
  ggplot2::ggplot(obs, ggplot2::aes(x = sample, y = value)) +
    ggplot2::geom_boxplot(outlier.size = 0.6, linewidth = 0.4, fill = "white", colour = "black") +
    ggplot2::labs(x = NULL, y = "log2 intensity (observed cells)", title = "Intensity distributions") +
    prism_figure_theme()
}

# Binary observed / missing map: observed cells light grey, missing cells black (never an intensity colour scale).
# Features are ordered by their missing fraction; the per-sample missing percentage is in the column labels.
catalogue_missingness <- function(source) {
  if (!"observed" %in% names(source)) catalogue_refuse("E_FIGURE_SOURCE_IMPUTED", "the source has no observed mask")
  feature_missing <- tapply(!source$observed, source$feature, mean)
  features <- names(sort(feature_missing, decreasing = TRUE))
  samples <- unique(source$sample)
  sample_pct <- tapply(!source$observed, source$sample, mean)[samples] * 100
  source$state <- ifelse(source$observed, "Observed", "Missing")
  source$feature <- factor(source$feature, levels = rev(features))
  source$sample <- factor(source$sample, levels = samples,
                          labels = sprintf("%s\n%.0f%%", samples, sample_pct))
  ggplot2::ggplot(source, ggplot2::aes(x = sample, y = feature, fill = state)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = c(Observed = "#D9D9D9", Missing = "black"), breaks = c("Observed", "Missing"),
                               name = "Observed / Missing") +
    ggplot2::labs(x = "sample (missing %)", y = "feature (ordered by missing fraction)", title = "Missingness") +
    prism_figure_theme() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5), axis.text.x = ggplot2::element_text(size = 7, face = "plain"))
}

catalogue_qc_correlation <- function(source) {
  if (!"observed" %in% names(source)) catalogue_refuse("E_FIGURE_SOURCE_IMPUTED", "the correlation needs the observed mask; an imputed matrix is never correlated")
  values <- catalogue_wide(source, "feature", "sample", "value")
  observed <- !is.na(values)
  cor_m <- catalogue_correlations(values, observed)
  long <- data.frame(a = rep(rownames(cor_m), ncol(cor_m)), b = rep(colnames(cor_m), each = nrow(cor_m)), r = as.vector(cor_m))
  long$a <- factor(long$a, levels = rownames(cor_m))
  long$b <- factor(long$b, levels = rev(colnames(cor_m)))
  ggplot2::ggplot(long, ggplot2::aes(a, b, fill = r)) + ggplot2::geom_tile(colour = "white") +
    ggplot2::scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0, limits = c(-1, 1), na.value = "grey90", name = "Pearson r") +
    ggplot2::labs(x = NULL, y = NULL, title = "Sample correlations (observed cells)") + prism_figure_theme()
}

# PCA scores with the variance explained on each axis ("PC1 (xx.x%)", from the same prcomp as the scores).
catalogue_qc_pca <- function(source, variant = FALSE) {
  groups <- unique(source$group)
  small <- names(which(table(source$group) < 3))
  pct1 <- source$PC1_var_pct[1]
  pct2 <- source$PC2_var_pct[1]
  ellipses <- do.call(rbind, lapply(setdiff(groups, small), function(g) {
    d <- source[source$group == g, ]
    cbind(catalogue_ellipse(d$PC1, d$PC2), group = g)
  }))
  note <- if (length(small)) sprintf("no ellipse for group(s) %s (fewer than three samples)", paste(small, collapse = ", ")) else NULL
  p <- ggplot2::ggplot(source, ggplot2::aes(PC1, PC2, colour = group)) +
    ggplot2::geom_point(size = 2.4, alpha = 0.85) +
    ggplot2::scale_colour_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#009E73")[seq_along(groups)]) +
    ggplot2::labs(x = sprintf("PC1 (%.1f%%)", pct1), y = sprintf("PC2 (%.1f%%)", pct2),
                  title = "PCA (95% ellipse per group)", caption = note) + prism_figure_theme()
  if (!is.null(ellipses) && nrow(ellipses)) p <- p + ggplot2::geom_path(data = ellipses, ggplot2::aes(x, y, colour = group), linewidth = 0.5, inherit.aes = FALSE)
  p
}

catalogue_volcano <- function(source, labels, q_cutoff = 0.05) {
  if (!is.null(labels)) {
    absent <- setdiff(labels, source$feature)
    if (length(absent)) catalogue_refuse("E_FIGURE_SOURCE_MISSING", sprintf("label(s) %s are not in the source table", paste(absent, collapse = ", ")))
    if (length(labels) > 10) catalogue_refuse("E_FIGURE_SET_MISMATCH", "at most 10 labels per volcano (D-74)")
  }
  source$class <- ifelse(source$q < q_cutoff & source$effect > 0, "up", ifelse(source$q < q_cutoff & source$effect < 0, "down", "ns"))
  source$neglog10P <- -log10(source$P)
  p <- ggplot2::ggplot(source, ggplot2::aes(effect, neglog10P, colour = class)) +
    ggplot2::geom_point(size = 2, alpha = 0.85) +
    ggplot2::scale_colour_manual(values = catalogue_colours()[c("up", "down", "ns")], breaks = c("up", "down", "ns")) +
    ggplot2::labs(x = "log2 fold change", y = "-log10 P", title = "Volcano (family q below cutoff)") + prism_figure_theme()
  if (!is.null(labels)) {
    lab <- source[source$feature %in% labels, ]
    p <- p + ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = feature), colour = "black", size = 3.6, box.padding = 0.6,
                                      point.padding = 0.3, min.segment.length = 0, max.overlaps = 30, seed = 1, show.legend = FALSE)
  }
  p
}

catalogue_ma <- function(source) {
  source$A <- (source$mean_numerator + source$mean_denominator) / 2
  source$M <- source$mean_numerator - source$mean_denominator
  ggplot2::ggplot(source, ggplot2::aes(A, M)) + ggplot2::geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
    ggplot2::geom_point(size = 1.6, alpha = 0.8, colour = "#0072B2") + ggplot2::facet_wrap(~contrast) +
    ggplot2::labs(x = "A = mean of the two groups (log2)", y = "M = numerator - denominator (log2)", title = "MA plot") +
    prism_figure_theme() +
    ggplot2::theme(panel.spacing = grid::unit(8, "mm"), axis.text.x = ggplot2::element_text(size = 9))
}

catalogue_dotplot_brackets <- function(values, brackets) {
  values$group <- factor(values$group, levels = unique(values$group))
  p <- ggplot2::ggplot(values, ggplot2::aes(group, value, colour = group)) +
    ggplot2::geom_point(size = 2.4, alpha = 0.85, show.legend = FALSE, position = ggplot2::position_jitter(width = 0.18, height = 0, seed = 2026)) +
    ggplot2::stat_summary(fun = mean, geom = "crossbar", width = 0.45, linewidth = 0.35, colour = "black") +
    ggplot2::stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x) - stats::sd(x), ymax = mean(x) + stats::sd(x)),
                          geom = "errorbar", width = 0.22, linewidth = 0.6, colour = "black") +
    ggplot2::scale_colour_manual(values = unname(catalogue_colours()[levels(values$group)])) +
    ggplot2::labs(x = NULL, y = "log2 abundance", title = NULL) + prism_figure_theme()
  p <- p + ggplot2::facet_wrap(~protein, scales = "free_y")
  if (!is.null(brackets) && nrow(brackets)) {
    # brackets are drawn from the family table: x positions are the group positions, the label is the FR-198 annotation
    levels_g <- levels(values$group)
    tip <- 0.02 * diff(range(values$value))
    br <- brackets
    br$x1 <- match(br$group1, levels_g); br$x2 <- match(br$group2, levels_g)
    if (any(is.na(br$x1) | is.na(br$x2))) catalogue_refuse("E_FIGURE_SOURCE_MISSING", "a bracket names a group that is not in the figure")
    segments <- rbind(
      data.frame(protein = br$protein, x = br$x1, xend = br$x2, y = br$y.position, yend = br$y.position),
      data.frame(protein = br$protein, x = br$x1, xend = br$x1, y = br$y.position - tip, yend = br$y.position),
      data.frame(protein = br$protein, x = br$x2, xend = br$x2, y = br$y.position - tip, yend = br$y.position))
    labels <- data.frame(protein = br$protein, x = (br$x1 + br$x2) / 2, y = br$y.position + tip, label = br$label)
    p <- p + ggplot2::geom_segment(data = segments, ggplot2::aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE, linewidth = 0.5) +
      ggplot2::geom_text(data = labels, ggplot2::aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 3.2, family = "Arial", vjust = 0)
  }
  p + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8, angle = 30, hjust = 1))
}

# Clustered DEP heatmap: row z-scores (sample SD), complete linkage on Euclidean distance (stats::hclust). A group colour bar
# sits above the columns, in the stable group colours. Cell colours are computed here (blue low, white 0, orange high) so that
# the bar and the heatmap can share one fill scale.
catalogue_dep_heatmap <- function(source) {
  if (any(!source$declared)) catalogue_refuse("E_FIGURE_SET_MISMATCH", "the heatmap may only show the declared DEP set")
  values <- catalogue_wide(source, "feature", "sample", "value")
  scaled <- t(scale(t(values)))
  hc <- stats::hclust(stats::dist(scaled, method = "euclidean"), method = "complete")
  order <- rownames(scaled)[hc$order]
  samples <- colnames(scaled)
  groups <- unique(source[, c("sample", "group")])
  group_of <- stats::setNames(groups$group, groups$sample)
  limit <- max(abs(scaled), na.rm = TRUE)
  pal <- scales::div_gradient_pal("#0072B2", "white", "#D55E00")
  cell_col <- pal(scales::rescale(as.vector(scaled), from = c(-limit, limit), to = c(0, 1)))
  heat <- data.frame(panel = "row z-scores", feature = rep(rownames(scaled), ncol(scaled)),
                     sample = rep(samples, each = nrow(scaled)), fill = cell_col, stringsAsFactors = FALSE)
  heat$feature <- factor(heat$feature, levels = rev(order))
  bar <- data.frame(panel = "group", feature = "group", sample = samples,
                    fill = unname(catalogue_colours()[group_of[samples]]), stringsAsFactors = FALSE)
  all <- rbind(heat, bar[, names(heat)])
  all$sample <- factor(all$sample, levels = samples)
  all$panel <- factor(all$panel, levels = c("group", "row z-scores"))
  ggplot2::ggplot(all, ggplot2::aes(sample, feature, fill = fill)) + ggplot2::geom_tile() +
    ggplot2::scale_fill_identity() +
    ggplot2::facet_grid(panel ~ ., scales = "free_y", space = "free_y") +
    ggplot2::labs(x = NULL, y = NULL, title = "Clustered DEP heatmap",
                  caption = "row z-scores: blue low, orange high; bar above the columns: group") +
    prism_figure_theme() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5), strip.text = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank())
}

# Region counts of the exclusive partition by brute-force enumeration of the 2^k - 1 patterns (SM31).
catalogue_regions <- function(membership) {
  k <- ncol(membership)
  patterns <- expand.grid(rep(list(c(0L, 1L)), k))[-1, , drop = FALSE]
  counts <- apply(patterns, 1, function(p) sum(apply(membership, 1, function(row) all(row == p))))
  data.frame(pattern = apply(patterns, 1, paste, collapse = ""), count = as.integer(counts), stringsAsFactors = FALSE)
}

# Standard circle layout for k <= 3 sets: unit radius, centres 1.2 apart (substantial overlap).
catalogue_venn_layout <- function(k) {
  switch(as.character(k),
    "2" = data.frame(x = c(-0.6, 0.6), y = c(0, 0)),
    "3" = data.frame(x = c(-0.6, 0.6, 0), y = c(-0.35, -0.35, 0.7)),
    catalogue_refuse("E_VENN_K", sprintf("a circle Venn needs k <= 3 sets; %d were given (UpSet is still produced)", k)))
}

catalogue_circle <- function(cx, cy, r = 1, n = 120) {
  theta <- seq(0, 2 * pi, length.out = n)
  data.frame(x = cx + r * cos(theta), y = cy + r * sin(theta))
}

# Venn of the exclusive regions: every region carries its count (the centroid of the grid cells with that exact membership),
# and every circle carries its set name with the set total outside it.
catalogue_venn <- function(source, set_names) {
  k <- length(set_names)
  layout <- catalogue_venn_layout(k)
  membership <- as.matrix(source[, set_names, drop = FALSE])
  regions <- catalogue_regions(membership)
  grid <- expand.grid(x = seq(-2.3, 2.3, by = 0.02), y = seq(-2.0, 2.4, by = 0.02))
  inside <- vapply(seq_len(k), function(i) (grid$x - layout$x[i])^2 + (grid$y - layout$y[i])^2 <= 1, logical(nrow(grid)))
  grid_pattern <- apply(inside * 1L, 1, paste, collapse = "")
  centre <- do.call(rbind, lapply(regions$pattern, function(pt) {
    cells <- grid_pattern == pt
    if (!any(cells)) catalogue_refuse("E_FIGURE_OUTPUT_MISSING", sprintf("region %s has no area in the layout", pt))
    data.frame(pattern = pt, x = mean(grid$x[cells]), y = mean(grid$y[cells]), stringsAsFactors = FALSE)
  }))
  labels <- merge(regions, centre, by = "pattern")
  circles <- do.call(rbind, lapply(seq_len(k), function(i) cbind(set = set_names[i], catalogue_circle(layout$x[i], layout$y[i]))))
  outward <- sqrt(layout$x^2 + layout$y^2)
  totals <- data.frame(label = sprintf("%s (n = %d)", set_names, as.integer(colSums(membership))),
                       x = layout$x + 1.35 * layout$x / outward,
                       y = layout$y + 1.35 * layout$y / outward, stringsAsFactors = FALSE)
  ggplot2::ggplot() +
    ggplot2::geom_polygon(data = circles, ggplot2::aes(x, y, group = set, fill = set), alpha = 0.3, colour = "black", linewidth = 0.4) +
    ggplot2::geom_text(data = labels, ggplot2::aes(x, y, label = count), size = 4) +
    ggplot2::geom_text(data = totals, ggplot2::aes(x, y, label = label), size = 3.8) +
    ggplot2::scale_fill_manual(values = c("#0072B2", "#E69F00", "#D55E00")[seq_len(k)], guide = "none") +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::labs(title = "Venn (exclusive regions)", x = NULL, y = NULL) + prism_figure_theme() +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), axis.line = ggplot2::element_blank())
}

# UpSet: exclusive intersection bars with their counts; a dot matrix whose rows are labelled with the set names; and the set
# size bars (the size of each set).
catalogue_upset <- function(source, set_names) {
  membership <- as.matrix(source[, set_names, drop = FALSE])
  k <- length(set_names)
  regions <- catalogue_regions(membership)
  regions <- regions[regions$count > 0, , drop = FALSE]          # zero-count intersections are not drawn
  regions <- regions[order(-regions$count, regions$pattern), ]
  regions$idx <- seq_len(nrow(regions))
  bars <- data.frame(panel = "intersections", x = regions$idx, y = regions$count, label = as.character(regions$count), stringsAsFactors = FALSE)
  dots <- expand.grid(x = regions$idx, y = seq_len(k))
  dots$on <- vapply(seq_len(nrow(dots)), function(i) substr(regions$pattern[dots$x[i]], dots$y[i], dots$y[i]) == "1", logical(1))
  dots$panel <- "matrix"
  set_labels <- data.frame(panel = "matrix", x = 0.3, y = seq_len(k), label = set_names, stringsAsFactors = FALSE)
  sizes <- data.frame(panel = "set size", x = as.numeric(colSums(membership)), y = seq_len(k),
                      label = sprintf("%s (%d)", set_names, as.integer(colSums(membership))), set = set_names, stringsAsFactors = FALSE)
  ggplot2::ggplot() +
    ggplot2::geom_col(data = bars, ggplot2::aes(x, y), fill = "black", width = 0.6) +
    ggplot2::geom_text(data = bars, ggplot2::aes(x, y, label = label), vjust = -0.4, size = 2.8) +
    ggplot2::geom_point(data = dots, ggplot2::aes(x, y, colour = on), size = 2) +
    ggplot2::geom_text(data = set_labels, ggplot2::aes(x, y, label = label), hjust = 1, size = 3.2) +
    ggplot2::geom_segment(data = sizes, ggplot2::aes(x = 0, xend = x, y = y, yend = y), linewidth = 4, colour = "black") +
    ggplot2::geom_text(data = sizes, ggplot2::aes(x, y, label = label), hjust = -0.1, size = 2.8) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey80"), guide = "none") +
    ggplot2::facet_wrap(~panel, ncol = 1, scales = "free") +
    ggplot2::labs(x = NULL, y = NULL, title = "UpSet (exclusive intersections)") + prism_figure_theme() +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), axis.line.x = ggplot2::element_blank()) +
    ggplot2::coord_cartesian(clip = "off")
}

catalogue_upset_venn <- function(source, set_names, venn = FALSE) {
  if (venn) return(catalogue_venn(source, set_names))
  catalogue_upset(source, set_names)
}

# The text labels drawn by a figure, from ggplot_build (every layer's label column), in layer order.
catalogue_text_labels <- function(plot) {
  built <- ggplot2::ggplot_build(plot)
  unlist(lapply(built$data, function(d) if ("label" %in% names(d)) as.character(d$label) else character(0)), use.names = FALSE)
}

# Renders one catalogue figure to PNG (300 dpi, Prism style). `options` is a named list (labels, brackets, q_cutoff,
# set_names, venn, variant).
render_catalogue_png <- function(figure_id, source_path, out_path, options = list(), width_mm = 120, height_mm = 100) {
  source <- utils::read.delim(source_path, stringsAsFactors = FALSE, na.strings = "NA", fileEncoding = "UTF-8")
  plot <- switch(figure_id,
    qc_intensity_distributions = catalogue_qc_intensity(source),
    qc_missingness_heatmap = catalogue_missingness(source),
    qc_sample_correlation_heatmap = catalogue_qc_correlation(source),
    qc_pca_ellipses = catalogue_qc_pca(source, isTRUE(options$variant)),
    diff_volcano_labelled = catalogue_volcano(source, options$labels, if (is.null(options$q_cutoff)) 0.05 else options$q_cutoff),
    diff_ma_plot = catalogue_ma(source),
    diff_protein_dotplot_brackets = catalogue_dotplot_brackets(source, if (is.null(options$brackets)) NULL else as.data.frame(options$brackets, stringsAsFactors = FALSE)),
    diff_dep_heatmap_clustered = catalogue_dep_heatmap(source),
    diff_upset_venn = catalogue_upset_venn(source, options$set_names, isTRUE(options$venn)),
    catalogue_refuse("E_FIGURE_SOURCE_MISSING", sprintf("no renderer for figure %s", figure_id)))
  render_prism_png(plot, out_path, width_mm = width_mm, height_mm = height_mm)
}
