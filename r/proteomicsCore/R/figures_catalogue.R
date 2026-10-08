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

# ---------------------------------------------------------------- stage 2 (V215-V222)

catalogue_cv_claim <- "cross_validated_nested"
catalogue_panel_claims <- c("in_sample", "cross_validated_nested", "fixed_panel_cv", "independently_validated")
catalogue_palette <- c("#0072B2", "#E69F00", "#D55E00", "#009E73", "#CC79A7", "#56B4E9", "#F0E442", "#000000")

# The one claim label of a figure, which must be a frozen label (E_BIOMARKER_CLAIM otherwise).
catalogue_claim <- function(source, allowed) {
  labels <- unique(source$claim_label)
  if (length(labels) != 1L || !labels %in% allowed) {
    catalogue_refuse("E_BIOMARKER_CLAIM", sprintf("the claim label is %s; the figure needs one of %s",
                                                  paste(sQuote(labels), collapse = ", "), paste(allowed, collapse = ", ")))
  }
  labels
}

# ROC with the stratified bootstrap band (V215). The curve and the band are the figure source (pROC, seed recorded).
catalogue_roc_ci_band <- function(source) {
  claim <- catalogue_claim(source, catalogue_cv_claim)
  if (any(source$score_source != "out_of_fold")) {
    catalogue_refuse("E_BIOMARKER_CLAIM", "the ROC is drawn from training-fold scores; only out-of-fold scores carry the cross-validated claim")
  }
  curve <- source[source$component == "curve", ]
  band <- source[source$component == "band", ]
  first <- source[1, ]
  ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = band, ggplot2::aes(x = fpr, ymin = band_lower, ymax = band_upper), fill = "#D55E00", alpha = 0.18) +
    ggplot2::geom_abline(linetype = "dashed", linewidth = 0.4, colour = "grey50") +
    ggplot2::geom_step(data = curve, ggplot2::aes(x = fpr, y = sensitivity), colour = "#D55E00", linewidth = 0.9, direction = "vh") +
    ggplot2::annotate("text", x = 0.42, y = 0.1, hjust = 0, size = 3.2, family = "Arial",
                      label = sprintf("AUC = %.3f\n95%% CI %.3f-%.3f", first$auc, first$auc_lower, first$auc_upper)) +
    ggplot2::coord_equal() +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
    ggplot2::labs(x = "1 - specificity", y = "sensitivity", title = "ROC with 95% CI band",
                  caption = sprintf("claim: %s\nbootstrap band: seed %d, B = %d", claim, as.integer(first$seed), as.integer(first$B))) +
    prism_figure_theme()
}

# Panel forest plot (V216): Case minus Control effects with their intervals, the claim label in the caption.
catalogue_forest <- function(source) {
  claim <- catalogue_claim(source, catalogue_panel_claims)
  source$feature <- factor(source$feature, levels = rev(source$feature))
  ggplot2::ggplot(source, ggplot2::aes(estimate, feature)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey50") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_lower, xmax = ci_upper), orientation = "y", width = 0.2, linewidth = 0.6) +
    ggplot2::geom_point(size = 2.6, colour = "#0072B2") +
    ggplot2::labs(x = "log2 difference (Case - Control), 95% CI", y = NULL, title = "Panel forest plot",
                  caption = sprintf("claim: %s", claim)) +
    prism_figure_theme()
}

# CV AUC against the permutation null (V217). The displayed P must equal (k + 1) / (B + 1) (E_FIGURE_PERMUTATION_P).
catalogue_cv_permutation <- function(source) {
  claim <- catalogue_claim(source, catalogue_cv_claim)
  cv <- source[source$series == "cv_repeat", ]
  nulls <- source[source$series == "null", ]
  B <- as.integer(source$B[1])
  if (nrow(nulls) != B) catalogue_refuse("E_FIGURE_SOURCE_MISSING", sprintf("the permutation null has %d values for B = %d", nrow(nulls), B))
  observed <- mean(cv$value)
  k <- as.integer(sum(nulls$value >= observed))
  p <- (k + 1) / (B + 1)
  if (any(abs(source$P_display - p) > 1e-12)) {
    catalogue_refuse("E_FIGURE_PERMUTATION_P", sprintf("the displayed permutation P differs from (k + 1) / (B + 1) = %g with k = %d and B = %d", p, k, B))
  }
  nulls$x <- "permutation null"
  cv$x <- "nested-CV repeats"
  ggplot2::ggplot() +
    ggplot2::geom_point(data = nulls, ggplot2::aes(x, value), colour = "grey60", size = 2, alpha = 0.85,
                        position = ggplot2::position_jitter(width = 0.12, height = 0, seed = 2026)) +
    ggplot2::geom_point(data = cv, ggplot2::aes(x, value), colour = "#0072B2", size = 2.4,
                        position = ggplot2::position_jitter(width = 0.12, height = 0, seed = 2026)) +
    ggplot2::geom_point(data = data.frame(x = "nested-CV repeats", value = observed), ggplot2::aes(x, value),
                        shape = 23, fill = "black", colour = "black", size = 3.5) +
    ggplot2::labs(x = NULL, y = "AUC", title = "CV AUC against permutation null",
                  caption = sprintf("permutation P = %.3f = (k + 1) / (B + 1), B = %d\nclaim: %s", p, B, claim)) +
    prism_figure_theme()
}

# Enrichment dot plot (V218): gene ratio k / n against -log10 q. Sets without a foreground hit are never plotted as enriched.
catalogue_enrichment <- function(source, show_sets = NULL) {
  shown <- if (is.null(show_sets)) source$set_id[source$k >= 1] else as.character(show_sets)
  unknown <- setdiff(shown, source$set_id)
  if (length(unknown)) catalogue_refuse("E_FIGURE_SOURCE_MISSING", sprintf("set(s) %s are not in the ORA source", paste(unknown, collapse = ", ")))
  plotted <- source[match(shown, source$set_id), ]
  empty <- plotted$set_id[plotted$k < 1]
  if (length(empty)) catalogue_refuse("E_FIGURE_SET_MISMATCH", sprintf("set(s) %s have no foreground hit; they are not enriched and are not plotted", paste(empty, collapse = ", ")))
  plotted$gene_ratio <- plotted$k / plotted$n
  plotted$neglog10_q <- -log10(plotted$q_value)
  plotted$significant <- plotted$q_value < 0.05
  # most significant set first in the rows; factor levels run bottom to top, so the most significant set is drawn at the top
  plotted <- plotted[order(plotted$q_value, plotted$set_id), ]
  plotted$set_id <- factor(plotted$set_id, levels = rev(plotted$set_id))
  ggplot2::ggplot(plotted, ggplot2::aes(gene_ratio, set_id, size = k, fill = neglog10_q)) +
    ggplot2::geom_point(shape = 21, colour = "grey40", stroke = 0.3) +
    ggplot2::geom_point(data = plotted[plotted$significant, ], ggplot2::aes(gene_ratio, set_id, size = k), shape = 21,
                        fill = NA, colour = "black", stroke = 1, show.legend = FALSE, inherit.aes = FALSE) +
    ggplot2::scale_fill_gradient(low = "grey85", high = "#D55E00", name = "-log10 q") +
    ggplot2::scale_size_continuous(range = c(2, 6), name = "genes in set (k)") +
    ggplot2::labs(x = "gene ratio (k / n)", y = NULL, title = "Enrichment dot plot (family q)",
                  caption = "outline: q < 0.05; top: most significant") +
    prism_figure_theme() +
    ggplot2::theme(legend.position = "right", legend.box = "vertical", legend.title = ggplot2::element_text(size = 9, face = "plain", hjust = 0))
}

# Module eigengene heatmap (V219): group means of the declared modules only (E_FIGURE_SET_MISMATCH otherwise).
catalogue_module_heatmap <- function(source, modules) {
  if (!length(modules)) catalogue_refuse("E_FIGURE_SET_MISMATCH", "no declared module set was given; no module is plotted")
  outside <- unique(source$module[!source$declared | !source$module %in% modules])
  if (length(outside)) catalogue_refuse("E_FIGURE_SET_MISMATCH", sprintf("module(s) %s are outside the declared set", paste(outside, collapse = ", ")))
  groups <- unique(source$group)
  source$module <- factor(source$module, levels = rev(modules))
  source$group <- factor(source$group, levels = groups)
  limit <- max(abs(source$eigengene))
  ggplot2::ggplot(source, ggplot2::aes(group, module, fill = eigengene)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", eigengene)), size = 3.2, family = "Arial") +
    ggplot2::scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0, limits = c(-limit, limit),
                                  name = "eigengene") +
    ggplot2::labs(x = NULL, y = NULL, title = "Module eigengenes (group means)") +
    prism_figure_theme()
}

# Deterministic layout of the network (presentation only): each module on a circle of radius 0.6 around the centre, its
# members on a circle of radius 0.25; modules in sorted order and members in source order. No random seed is used.
catalogue_network_layout <- function(nodes) {
  nodes <- as.data.frame(nodes, stringsAsFactors = FALSE)
  modules <- sort(unique(nodes$module))
  m <- length(modules)
  coords <- do.call(rbind, lapply(seq_along(modules), function(i) {
    members <- nodes[nodes$module == modules[i], , drop = FALSE]
    n <- nrow(members)
    cx <- if (m == 1L) 0 else 0.6 * cos(2 * pi * (i - 1) / m)
    cy <- if (m == 1L) 0 else 0.6 * sin(2 * pi * (i - 1) / m)
    r <- if (n == 1L) 0 else 0.25
    angle <- 2 * pi * (seq_len(n) - 1) / n + pi / 2
    data.frame(id = members$id, x = cx + r * cos(angle), y = cy + r * sin(angle), stringsAsFactors = FALSE)
  }))
  coords[match(nodes$id, coords$id), c("id", "x", "y")]
}

# Network graph (V220): the declared subnetwork's edges and nodes as in the snapshot; edges and nodes drawn in ggplot2.
catalogue_network_graph <- function(source) {
  nodes <- source[source$kind == "node", c("id", "module"), drop = FALSE]
  edges <- source[source$kind == "edge", , drop = FALSE]
  xy <- catalogue_network_layout(nodes)
  edges$x <- xy$x[match(edges$from, xy$id)]
  edges$y <- xy$y[match(edges$from, xy$id)]
  edges$xend <- xy$x[match(edges$to, xy$id)]
  edges$yend <- xy$y[match(edges$to, xy$id)]
  if (anyNA(c(edges$x, edges$xend))) catalogue_refuse("E_FIGURE_SOURCE_MISSING", "an edge has an endpoint that is not a node of the figure")
  xy$module <- nodes$module[match(xy$id, nodes$id)]
  modules <- sort(unique(xy$module))
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = edges, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), colour = "grey60", linewidth = 0.4) +
    ggplot2::geom_point(data = xy, ggplot2::aes(x, y, colour = module), size = 3) +
    ggplot2::geom_text(data = xy, ggplot2::aes(x, y, label = id), size = 2.4, vjust = -1.1, family = "Arial") +
    ggplot2::scale_colour_manual(values = stats::setNames(catalogue_palette[seq_along(modules)], modules)) +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::labs(x = NULL, y = NULL, title = "Declared subnetwork",
                  caption = "edges as in the interaction snapshot") +
    prism_figure_theme() +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), axis.line = ggplot2::element_blank())
}

# PCoA with group centroids (V221): centroids are the means of each group's own members.
catalogue_pcoa <- function(source) {
  groups <- unique(source$group)
  centroids <- stats::aggregate(cbind(axis1, axis2) ~ group, data = source, FUN = mean)
  names(centroids) <- c("group", "cx", "cy")
  linked <- merge(source, centroids, by = "group")
  palette <- stats::setNames(catalogue_palette[seq_along(groups)], groups)
  ggplot2::ggplot(source, ggplot2::aes(axis1, axis2, colour = group)) +
    ggplot2::geom_segment(data = linked, ggplot2::aes(x = axis1, y = axis2, xend = cx, yend = cy), colour = "grey75",
                          linewidth = 0.3, inherit.aes = FALSE) +
    ggplot2::geom_point(size = 2.2, alpha = 0.85) +
    ggplot2::geom_point(data = centroids, ggplot2::aes(cx, cy, fill = group), shape = 23, size = 4, colour = "black", inherit.aes = FALSE) +
    ggplot2::scale_colour_manual(values = palette, name = "group") +
    ggplot2::scale_fill_manual(values = palette, name = "group") +
    ggplot2::labs(x = sprintf("PCoA 1 (%.1f%%)", source$percent_axis1[1]), y = sprintf("PCoA 2 (%.1f%%)", source$percent_axis2[1]),
                  title = "PCoA with group centroids", caption = "diamonds: group centroids (means of each group's members)") +
    prism_figure_theme() +
    ggplot2::theme(legend.title = ggplot2::element_text(size = 9, face = "plain", hjust = 0))
}

# PERMDISP plot (V222): distances to the group centroids, with the PERMDISP P, B and seed. A P without B or seed is refused.
catalogue_permdisp <- function(source) {
  if (anyNA(source$permutations) || anyNA(source$seed)) {
    catalogue_refuse("E_FIGURE_SOURCE_MISSING", "the PERMDISP P is shown without its number of permutations B and its seed")
  }
  B <- as.integer(source$permutations[1])
  seed <- as.integer(source$seed[1])
  groups <- unique(source$group)
  source$group <- factor(source$group, levels = groups)
  palette <- stats::setNames(catalogue_palette[seq_along(groups)], groups)
  ggplot2::ggplot(source, ggplot2::aes(group, distance_to_centroid, colour = group)) +
    ggplot2::geom_boxplot(colour = "black", fill = NA, outlier.shape = NA, linewidth = 0.4) +
    ggplot2::geom_point(size = 1.8, alpha = 0.85, position = ggplot2::position_jitter(width = 0.15, height = 0, seed = 2026)) +
    ggplot2::scale_colour_manual(values = palette) +
    ggplot2::labs(x = NULL, y = "distance to group centroid", title = "PERMDISP: distances to group centroids",
                  caption = sprintf("PERMDISP P = %.3f; B = %d permutations; seed %d", source$permdisp_p[1], B, seed)) +
    prism_figure_theme() +
    ggplot2::theme(legend.position = "none")
}

# Renders one catalogue figure to PNG (300 dpi, Prism style). `options` is a named list (labels, brackets, q_cutoff,
# set_names, venn, variant; stage 2: show_sets, modules).
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
    bio_roc_ci_band = catalogue_roc_ci_band(source),
    bio_panel_forest = catalogue_forest(source),
    bio_cv_auc_vs_permutation = catalogue_cv_permutation(source),
    path_enrichment_dotplot = catalogue_enrichment(source, options$show_sets),
    net_module_heatmap = catalogue_module_heatmap(source, options$modules),
    net_network_graph = catalogue_network_graph(source),
    perm_pcoa_centroids = catalogue_pcoa(source),
    perm_permdisp = catalogue_permdisp(source),
    catalogue_refuse("E_FIGURE_SOURCE_MISSING", sprintf("no renderer for figure %s", figure_id)))
  render_prism_png(plot, out_path, width_mm = width_mm, height_mm = height_mm)
}

# ---------------------------------------------------------------- retrofit of existing figures (D-72)
#
# The registered figure-source tables of the R10b, R14 producers are drawn here with the R16c theme. The Python retrofit rule
# (figures/retrofit.py) chooses the kind and the column names; nothing is recomputed except the transformations the rule names
# (absolute value, group means of the eigengenes). Refusals use the typed codes of the contract.

retro_flag <- function(values, flag_values) {
  if (is.null(flag_values) || !length(flag_values)) return(rep("other", length(values)))
  ifelse(as.character(values) %in% as.character(flag_values), "flag", "other")
}

retro_plot_scatter <- function(source, o) {
  source$.x <- as.numeric(source[[o$x]])
  source$.y <- as.numeric(source[[o$y]])
  if (isTRUE(o$abs_y)) source$.y <- abs(source$.y)
  source <- source[is.finite(source$.x) & is.finite(source$.y), , drop = FALSE]
  source$.flag <- if (!is.null(o$flag_col)) retro_flag(source[[o$flag_col]], o$flag_values) else rep("other", nrow(source))
  p <- ggplot2::ggplot(source, ggplot2::aes(.x, .y, colour = .flag)) +
    ggplot2::geom_point(size = 2, alpha = 0.85) +
    ggplot2::scale_colour_manual(values = c(flag = "#D55E00", other = "grey60"), guide = "none") +
    ggplot2::labs(x = o$xlab, y = o$ylab, title = o$title, caption = o$caption) + prism_figure_theme()
  if (isTRUE(o$diagonal)) p <- p + ggplot2::geom_abline(linetype = "dashed", linewidth = 0.4, colour = "grey50")
  p
}

retro_plot_curve_band <- function(source, o) {
  band <- utils::read.delim(o$band_path, stringsAsFactors = FALSE, na.strings = "NA", fileEncoding = "UTF-8")
  source$.x <- as.numeric(source[[o$x]])
  source$.y <- as.numeric(source[[o$y]])
  band$.x <- as.numeric(band[[o$band_x]])
  band$.lo <- as.numeric(band[[o$band_lo]])
  band$.hi <- as.numeric(band[[o$band_hi]])
  ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = band, ggplot2::aes(.x, ymin = .lo, ymax = .hi), fill = "#D55E00", alpha = 0.18) +
    ggplot2::geom_abline(linetype = "dashed", linewidth = 0.4, colour = "grey50") +
    ggplot2::geom_step(data = source, ggplot2::aes(.x, .y), colour = "#D55E00", linewidth = 0.9, direction = "vh") +
    ggplot2::coord_equal() +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
    ggplot2::labs(x = o$xlab, y = o$ylab, title = o$title, caption = o$caption) + prism_figure_theme()
}

retro_plot_histogram <- function(source, o) {
  source$.x <- as.numeric(source[[o$x]])
  source <- source[is.finite(source$.x), , drop = FALSE]
  p <- ggplot2::ggplot(source, ggplot2::aes(.x)) +
    ggplot2::geom_histogram(bins = 30, fill = "grey75", colour = "white") +
    ggplot2::labs(x = o$xlab, y = "count", title = o$title, caption = o$caption) + prism_figure_theme()
  if (!is.null(o$ref)) p <- p + ggplot2::geom_vline(xintercept = as.numeric(o$ref), colour = "#D55E00", linewidth = 0.8)
  p
}

retro_plot_strip <- function(source, o) {
  source$.x <- as.numeric(source[[o$x]])
  source <- source[is.finite(source$.x), , drop = FALSE]
  source$.g <- "repeats"
  ggplot2::ggplot(source, ggplot2::aes(.g, .x)) +
    ggplot2::geom_point(colour = "#0072B2", size = 2.4, alpha = 0.85, position = ggplot2::position_jitter(width = 0.08, height = 0, seed = 2026)) +
    ggplot2::labs(x = NULL, y = o$xlab, title = o$title, caption = o$caption) + prism_figure_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
}

# Horizontal bars with the label on the axis, in the source order (the source is already ordered by the producer).
retro_plot_bar <- function(source, o, descending = TRUE) {
  source$.v <- as.numeric(source[[o$value]])
  source$.l <- as.character(source[[o$label]])
  source <- source[is.finite(source$.v), , drop = FALSE]
  source$.l <- factor(source$.l, levels = rev(source$.l))
  ggplot2::ggplot(source, ggplot2::aes(.v, .l)) +
    ggplot2::geom_col(fill = "#4c78a8", width = 0.6) +
    ggplot2::labs(x = o$xlab, y = NULL, title = o$title, caption = o$caption) + prism_figure_theme()
}

# Regions with their counts (UpSet-like bars), labelled by the defining sets, largest first.
retro_plot_regions <- function(source, o) {
  source$.v <- as.numeric(source[[o$value]])
  source$.l <- as.character(source[[o$label]])
  source <- source[order(-source$.v, source$.l), , drop = FALSE]
  source$.l <- factor(source$.l, levels = rev(source$.l))
  ggplot2::ggplot(source, ggplot2::aes(.v, .l)) +
    ggplot2::geom_col(fill = "black", width = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = .v), hjust = -0.2, size = 2.8) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(x = o$xlab, y = NULL, title = o$title, caption = o$caption) + prism_figure_theme()
}

# Tile heatmap of a long table: x and y are identifiers, fill the value (NA cells in grey).
retro_plot_heat <- function(long, o) {
  long$.x <- factor(as.character(long$.x), levels = unique(as.character(long$.x)))
  long$.y <- factor(as.character(long$.y), levels = rev(unique(as.character(long$.y))))
  limit <- max(abs(long$.v), na.rm = TRUE)
  if (!is.finite(limit) || limit == 0) limit <- 1
  ggplot2::ggplot(long, ggplot2::aes(.x, .y, fill = .v)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0, limits = c(-limit, limit),
                                  na.value = "grey90", name = o$fill_title) +
    ggplot2::labs(x = o$xlab, y = o$ylab, title = o$title, caption = o$caption) + prism_figure_theme() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6))
}

retro_plot_association <- function(source, o) {
  long <- data.frame(.x = as.character(source[[o$x]]), .y = as.character(source[[o$y]]), .v = as.numeric(source[[o$fill]]), stringsAsFactors = FALSE)
  retro_plot_heat(long, o)
}

# Group means of the eigengenes (one column per module): the modules are the value columns of the source.
retro_plot_eigen <- function(source, o) {
  modules <- setdiff(names(source), c("observation_id", "group", "claim_label"))
  groups <- unique(source$group)
  long <- do.call(rbind, lapply(modules, function(m) data.frame(.x = groups, .y = m,
    .v = vapply(groups, function(g) mean(as.numeric(source[[m]])[source$group == g]), numeric(1)), stringsAsFactors = FALSE)))
  retro_plot_heat(long, o)
}

# Node-link drawing with the layout of the source (x, y), edges read from the companion table of the same stage.
retro_plot_network <- function(source, o) {
  edges <- utils::read.delim(o$edges_path, stringsAsFactors = FALSE, na.strings = "NA", fileEncoding = "UTF-8")
  source$.x <- as.numeric(source[[o$x]])
  source$.y <- as.numeric(source[[o$y]])
  ids <- as.character(source[[o$label]])
  edges$.x1 <- source$.x[match(edges[[o$edge_a]], ids)]
  edges$.y1 <- source$.y[match(edges[[o$edge_a]], ids)]
  edges$.x2 <- source$.x[match(edges[[o$edge_b]], ids)]
  edges$.y2 <- source$.y[match(edges[[o$edge_b]], ids)]
  edges <- edges[is.finite(edges$.x1) & is.finite(edges$.x2), , drop = FALSE]
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = edges, ggplot2::aes(x = .x1, y = .y1, xend = .x2, yend = .y2), colour = "grey60", linewidth = 0.4) +
    ggplot2::geom_point(data = source, ggplot2::aes(.x, .y), colour = "#0072B2", size = 3) +
    ggplot2::geom_text(data = source, ggplot2::aes(.x, .y, label = .data[[o$label]]), size = 2.4, vjust = -1.1, family = "Arial") +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::labs(x = NULL, y = NULL, title = o$title, caption = o$caption) + prism_figure_theme() +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), axis.line = ggplot2::element_blank())
}

# Renders one retrofitted figure from its registered source copy. `kind` and `options` come from the Python retrofit rule.
render_retrofit_png <- function(kind, source_path, out_path, options = list(), width_mm = 120, height_mm = 100) {
  source <- utils::read.delim(source_path, stringsAsFactors = FALSE, na.strings = "NA", fileEncoding = "UTF-8")
  plot <- switch(kind,
    scatter = retro_plot_scatter(source, options),
    curve_band = retro_plot_curve_band(source, options),
    histogram = retro_plot_histogram(source, options),
    strip = retro_plot_strip(source, options),
    bar = retro_plot_bar(source, options),
    regions = retro_plot_regions(source, options),
    heat_table = retro_plot_association(source, options),
    eigen_means = retro_plot_eigen(source, options),
    network = retro_plot_network(source, options),
    catalogue_refuse("E_FIGURE_SOURCE_MISSING", sprintf("no retrofit renderer for kind %s", kind)))
  render_prism_png(plot, out_path, width_mm = width_mm, height_mm = height_mm)
}
