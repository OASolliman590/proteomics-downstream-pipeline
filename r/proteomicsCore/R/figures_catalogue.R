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

catalogue_qc_missingness <- function(source) {
  if (!"observed" %in% names(source)) catalogue_refuse("E_FIGURE_SOURCE_IMPUTED", "the source has no observed mask")
  source$shown <- ifelse(source$observed, source$value, NA_real_)
  source$feature <- factor(source$feature, levels = rev(unique(source$feature)))
  source$sample <- factor(source$sample, levels = unique(source$sample))
  ggplot2::ggplot(source, ggplot2::aes(x = sample, y = feature, fill = shown)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient(low = "#F0E442", high = "#0072B2", na.value = "grey20", name = "log2 intensity (grey = missing)") +
    ggplot2::labs(x = NULL, y = NULL, title = "Missingness") +
    prism_figure_theme() + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6))
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

catalogue_qc_pca <- function(source, variant = FALSE) {
  groups <- unique(source$group)
  small <- names(which(table(source$group) < 3))
  ellipses <- do.call(rbind, lapply(setdiff(groups, small), function(g) {
    d <- source[source$group == g, ]
    cbind(catalogue_ellipse(d$PC1, d$PC2), group = g)
  }))
  note <- if (length(small)) sprintf("no ellipse for group(s) %s (fewer than three samples)", paste(small, collapse = ", ")) else NULL
  p <- ggplot2::ggplot(source, ggplot2::aes(PC1, PC2, colour = group)) +
    ggplot2::geom_point(size = 2.4, alpha = 0.85) +
    ggplot2::scale_colour_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#009E73")[seq_along(groups)]) +
    ggplot2::labs(x = "PC1", y = "PC2", title = "PCA (95% ellipse per group)", caption = note) + prism_figure_theme()
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
    ggplot2::labs(x = "A = mean of the two groups (log2)", y = "M = numerator - denominator (log2)", title = "MA plot") + prism_figure_theme()
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

catalogue_dep_heatmap <- function(source) {
  if (any(!source$declared)) catalogue_refuse("E_FIGURE_SET_MISMATCH", "the heatmap may only show the declared DEP set")
  values <- catalogue_wide(source, "feature", "sample", "value")
  scaled <- t(scale(t(values)))
  hc <- stats::hclust(stats::dist(scaled, method = "euclidean"), method = "complete")
  order <- rownames(scaled)[hc$order]
  long <- data.frame(feature = rep(rownames(scaled), ncol(scaled)), sample = rep(colnames(scaled), each = nrow(scaled)), z = as.vector(scaled))
  long$feature <- factor(long$feature, levels = rev(order))
  ggplot2::ggplot(long, ggplot2::aes(sample, feature, fill = z)) + ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0, name = "row z") +
    ggplot2::labs(x = NULL, y = NULL, title = "Clustered DEP heatmap") + prism_figure_theme() + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5))
}

# Region counts of the exclusive partition by brute-force enumeration of the 2^k - 1 patterns (SM31).
catalogue_regions <- function(membership) {
  k <- ncol(membership)
  patterns <- expand.grid(rep(list(c(0L, 1L)), k))[-1, , drop = FALSE]
  counts <- apply(patterns, 1, function(p) sum(apply(membership, 1, function(row) all(row == p))))
  data.frame(pattern = apply(patterns, 1, paste, collapse = ""), count = as.integer(counts), stringsAsFactors = FALSE)
}

catalogue_upset_venn <- function(source, set_names, venn = FALSE) {
  membership <- as.matrix(source[, set_names, drop = FALSE])
  k <- length(set_names)
  regions <- catalogue_regions(membership)
  regions <- regions[order(-regions$count), ]
  if (venn) {
    if (k > 3) catalogue_refuse("E_VENN_K", sprintf("a circle Venn needs k <= 3 sets; %d were given (UpSet is still produced)", k))
    centres <- data.frame(set = set_names, x = c(0, 1, 0.5)[seq_len(k)], y = c(0, 0, 0.87)[seq_len(k)])
    circles <- do.call(rbind, lapply(seq_len(k), function(i) cbind(set = set_names[i], catalogue_circle(centres$x[i], centres$y[i]))))
    return(ggplot2::ggplot(circles, ggplot2::aes(x, y, group = set, fill = set)) + ggplot2::geom_polygon(alpha = 0.3, colour = "black", linewidth = 0.4) +
      ggplot2::coord_equal() + ggplot2::labs(title = "Venn (exclusive regions)", x = NULL, y = NULL) + prism_figure_theme() +
      ggplot2::scale_fill_manual(values = c("#0072B2", "#E69F00", "#D55E00")[seq_len(k)]) + ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), axis.line = ggplot2::element_blank()))
  }
  regions$idx <- seq_len(nrow(regions))
  bars <- regions; bars$panel <- "intersections"; bars$y <- bars$count
  dots <- expand.grid(idx = regions$idx, set = seq_len(k))
  dots$on <- as.integer(substr(regions$pattern[dots$idx], dots$set, dots$set)) == 1L
  dots$panel <- "matrix"; dots$y <- dots$set; dots$setname <- set_names[dots$set]
  ggplot2::ggplot() +
    ggplot2::geom_col(data = bars, ggplot2::aes(factor(idx), y), fill = "black", width = 0.6) +
    ggplot2::geom_point(data = dots, ggplot2::aes(factor(idx), y, colour = on), size = 2) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey80"), guide = "none") +
    ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_y") +
    ggplot2::labs(x = NULL, y = NULL, title = "UpSet (exclusive intersections)") + prism_figure_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
}

catalogue_circle <- function(cx, cy, r = 0.55, n = 80) {
  theta <- seq(0, 2 * pi, length.out = n)
  data.frame(x = cx + r * cos(theta), y = cy + r * sin(theta))
}

# Renders one catalogue figure to PNG (300 dpi, Prism style). `options` is a named list (labels, brackets, q_cutoff,
# set_names, venn, variant).
render_catalogue_png <- function(figure_id, source_path, out_path, options = list(), width_mm = 120, height_mm = 100) {
  source <- utils::read.delim(source_path, stringsAsFactors = FALSE, na.strings = "NA")
  plot <- switch(figure_id,
    qc_intensity_distributions = catalogue_qc_intensity(source),
    qc_missingness_heatmap = catalogue_qc_missingness(source),
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
