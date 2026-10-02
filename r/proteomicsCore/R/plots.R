# Publication figures from saved figure-source tables (packet R10b, FR-097).
# Every plotted coordinate comes from the source TSV written by the report;
# no rounding, no extra points, no significance stars.

plot_figure_source <- function(source, kind, target, format) {
  d <- .pc_read_tsv(source)
  num <- function(x) suppressWarnings(as.numeric(ifelse(x == "NA", NA, x)))
  device <- switch(format, pdf = function() grDevices::pdf(target, width = 6, height = 4.5, useDingbats = FALSE),
                   png = function() grDevices::png(target, width = 6, height = 4.5, units = "in", res = 200),
                   svg = function() { if (!requireNamespace("svglite", quietly = TRUE)) stop("E_DEPENDENCY_UNAVAILABLE: svglite", call. = FALSE); svglite::svglite(target, width = 6, height = 4.5) },
                   stop("E_FIGURE_FORMAT: unsupported figure format", call. = FALSE))
  device(); on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(4.5, 4.5, 3, 1))
  if (identical(kind, "volcano")) {
    x <- num(d$effect_log2); y <- num(d$neg_log10_p); hit <- d$q_at_or_below_cutoff == "true"
    graphics::plot(x, y, pch = 19, col = ifelse(hit, "#b2182b", "#7f7f7fB0"), xlab = "Effect (log2 difference)", ylab = "-log10 P",
                   main = sprintf("%d tested endpoints; red: family q at or below cutoff", nrow(d)), cex.main = 0.85)
    graphics::abline(v = 0, lty = 3)
  } else stop("E_FIGURE_KIND: unknown figure kind", call. = FALSE)
  invisible(target)
}
