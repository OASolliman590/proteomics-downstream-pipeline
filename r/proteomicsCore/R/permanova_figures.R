# PERMANOVA figures (packet R13, FR-129).  Base graphics only; every figure is
# drawn from a published source table and written as PNG (plus requested PDF;
# SVG only when svglite is installed, otherwise recorded as NOT_RUN).

.pm_palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000")

.pm_devices <- function(out, stem, formats, width, height, draw) {
  files <- list(); records <- list()
  dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
  for (fmt in unique(c("png", formats))) {
    relative <- file.path("figures", paste0(stem, ".", fmt)); path <- file.path(out, relative)
    if (identical(fmt, "svg") && !requireNamespace("svglite", quietly = TRUE)) { records[[length(records) + 1L]] <- list(file = relative, state = "NOT_RUN", reason_code = "E_DEPENDENCY_UNAVAILABLE"); next }
    switch(fmt,
           png = grDevices::png(path, width = width, height = height, units = "in", res = 150),
           pdf = grDevices::pdf(path, width = width, height = height, useDingbats = FALSE),
           svg = svglite::svglite(path, width = width, height = height))
    ok <- tryCatch({ draw(); TRUE }, finally = grDevices::dev.off())
    files[[length(files) + 1L]] <- list(relative_path = relative, artifact_id = paste0("figure_", gsub("[^A-Za-z0-9_.-]", "_", stem), "_", fmt))
    records[[length(records) + 1L]] <- list(file = relative, state = "COMPLETED")
  }
  list(files = files, records = records)
}

.pm_label_p <- function(display) ifelse(grepl("^<", display), paste("p", display), paste("p =", display))

permanova_figures <- function(out, tests, ordination, axes, centroids, nulls, feature_r2, random_null, context, formats) {
  formats <- intersect(formats, c("png", "pdf", "svg"))
  all_files <- list(); all_records <- list()
  add <- function(x) { all_files <<- c(all_files, x$files); all_records <<- c(all_records, x$records) }
  for (set_id in unique(ordination$feature_set_id)) {
    o <- ordination[ordination$feature_set_id == set_id, ]; ax <- axes[axes$feature_set_id == set_id, ]; ce <- centroids[centroids$feature_set_id == set_id, ]
    g <- tests[tests$feature_set_id == set_id & tests$analysis == "global", ][1, ]
    groups <- unique(o$group); cols <- stats::setNames(.pm_palette[seq_along(groups)], groups)
    add(.pm_devices(out, paste0("ordination_", set_id), formats, 6.5, 5.8, function() {
      graphics::par(mar = c(6, 4.5, 4, 1))
      graphics::plot(o$axis1, o$axis2, type = "n", xlab = sprintf("PCoA 1 (%.0f%%)", ax$percent_positive_eigenvalue[1]), ylab = sprintf("PCoA 2 (%.0f%%)", ax$percent_positive_eigenvalue[2]),
                     main = sprintf("%s: PERMANOVA R2 = %.2f, %s\nPERMDISP %s", set_id, g$r2, .pm_label_p(g$p_display), .pm_label_p(g$permdisp_p_display)), cex.main = 0.9)
      graphics::segments(o$axis1, o$axis2, o$centroid_axis1, o$centroid_axis2, col = grDevices::adjustcolor(cols[o$group], 0.4))
      graphics::points(o$axis1, o$axis2, pch = 19, col = cols[o$group])
      graphics::points(ce$axis1, ce$axis2, pch = 23, bg = cols[ce$group], col = "black", cex = 2)
      graphics::legend("topright", legend = groups, col = cols, pch = 19, bty = "n", cex = 0.8)
      graphics::mtext("Diamonds: group centroids; lines join each unit to its centroid", side = 1, line = 4.6, cex = 0.7)
    }))
    add(.pm_devices(out, paste0("permdisp_", set_id), formats, 5.5, 4.5, function() {
      graphics::par(mar = c(4, 4.5, 3, 1))
      graphics::boxplot(distance_to_centroid ~ factor(group, levels = groups), data = o, col = grDevices::adjustcolor(cols, 0.35), outline = FALSE,
                        xlab = "", ylab = "Distance to group centroid", main = sprintf("%s: PERMDISP %s", set_id, .pm_label_p(g$permdisp_p_display)), cex.main = 0.9)
      graphics::stripchart(distance_to_centroid ~ factor(group, levels = groups), data = o, vertical = TRUE, method = "jitter", pch = 19, col = cols, add = TRUE)
    }))
    fr <- feature_r2[feature_r2$feature_set_id == set_id, ]
    fr <- fr[order(fr$univariate_group_r2), ]
    show <- utils::tail(fr, 40L)
    add(.pm_devices(out, paste0("feature_r2_", set_id), formats, 6, max(3, 0.18 * nrow(show) + 1.5), function() {
      graphics::par(mar = c(4, 8, 3, 1))
      graphics::barplot(100 * show$univariate_group_r2, names.arg = show$feature_id, horiz = TRUE, las = 1, col = "#7f7fbf", border = NA, cex.names = 0.6,
                        xlab = "Variance explained by group (%)", main = sprintf("%s: per-feature group R2 (top %d)", set_id, nrow(show)), cex.main = 0.9)
      graphics::abline(v = 100 * mean(fr$univariate_group_r2), lty = 2)
    }))
  }
  for (test_id in unique(paste(nulls$feature_set_id, nulls$test_id, sep = "|"))) {
    n <- nulls[paste(nulls$feature_set_id, nulls$test_id, sep = "|") == test_id, ]
    add(.pm_devices(out, paste0("null_", gsub("[^A-Za-z0-9_-]", "_", test_id)), formats, 5.5, 4, function() {
      graphics::par(mar = c(4, 4, 3, 1))
      rng <- range(c(n$pseudo_f, n$observed_pseudo_f[1]), finite = TRUE)
      graphics::hist(n$pseudo_f, breaks = 50, col = "grey75", border = "white", xlim = rng, xlab = "Pseudo-F under permutation", main = sprintf("%s (%d permutations)", sub("\\|", ": ", test_id), nrow(n)), cex.main = 0.85)
      graphics::abline(v = n$observed_pseudo_f[1], col = "#b2182b", lwd = 2)
    }))
  }
  add(.pm_devices(out, "r2_summary", formats, 7, 4.5, function() {
    t <- tests[tests$analysis %in% c("global", "pairwise", "group_adjusted", "covariate"), ]
    graphics::par(mar = c(8, 4.5, 3, 1))
    labels <- paste(t$feature_set_id, t$comparison, sep = ": ")
    mids <- graphics::barplot(100 * t$r2, names.arg = labels, las = 2, cex.names = 0.6, col = "#4c78a8", border = NA, ylab = "Variation explained (R2, %)",
                              main = "PERMANOVA R2 by test", ylim = c(0, max(100 * t$r2, 1) * 1.25), cex.main = 0.9)
    graphics::text(mids, 100 * t$r2, labels = ifelse(is.na(t$p_adjusted_display), .pm_label_p(t$p_display), paste("adj.", .pm_label_p(t$p_adjusted_display))), pos = 3, cex = 0.55)
  }))
  if (!is.null(random_null)) for (set_id in unique(random_null$feature_set_id)) {
    r <- random_null[random_null$feature_set_id == set_id, ]; ctx <- context[context$feature_set_id == set_id, ][1, ]
    add(.pm_devices(out, paste0("random_null_", set_id), formats, 5.5, 4, function() {
      graphics::par(mar = c(4, 4, 3, 1))
      rng <- range(c(r$r2, ctx$set_r2, ctx$best_possible_r2), finite = TRUE)
      graphics::hist(100 * r$r2, breaks = 50, col = "grey75", border = "white", xlim = 100 * rng, xlab = "Group R2 of random equal-size sets (%)",
                     main = sprintf("%s: in-sample R2 versus %d random sets", set_id, nrow(r)), cex.main = 0.85)
      graphics::abline(v = 100 * ctx$set_r2, col = "#b2182b", lwd = 2)
      if (is.finite(ctx$best_possible_r2)) graphics::abline(v = 100 * ctx$best_possible_r2, lty = 2)
    }))
  }
  list(files = all_files, records = all_records)
}
