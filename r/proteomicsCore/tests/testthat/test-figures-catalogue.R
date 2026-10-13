# R16e at the R level (V206-V214): the R functions of the catalogue, checked against base R. The Python acceptance test is
# tests/integration/test_figure_catalogue.py; the oracles here are base R: cor(), prcomp(), hclust(), brute-force regions,
# the chi-square contour and the typed refusals of the renderers.

ns <- asNamespace("proteomicsCore")
for (nm in ls(ns, pattern = "^catalogue_")) assign(nm, get(nm, envir = ns))

fixtures <- normalizePath(file.path(testthat::test_path(), "..", "..", "..", "..", "tests", "fixtures", "figures", "catalogue"), mustWork = TRUE)
read_fixture <- function(name) utils::read.delim(file.path(fixtures, name), stringsAsFactors = FALSE, na.strings = "NA")

testthat::test_that("V208 coverage-aware correlations equal cor() on the cells observed in both samples", {
  src <- read_fixture("qc_sample_correlation_source.tsv")
  values <- catalogue_wide(src, "feature", "sample", "value")
  observed <- !is.na(values)
  r <- catalogue_correlations(values, observed)
  for (a in colnames(values)) for (b in colnames(values)) {
    shared <- observed[, a] & observed[, b]
    want <- if (sum(shared) >= 3L) stats::cor(values[shared, a], values[shared, b]) else NA_real_
    testthat::expect_equal(r[a, b], want, tolerance = 1e-10)
  }
})

testthat::test_that("V208 negative: a source without the observed mask is refused as imputed", {
  src <- read_fixture("qc_sample_correlation_imputed.tsv")
  testthat::expect_error(catalogue_qc_correlation(src), "E_FIGURE_SOURCE_IMPUTED")
})

testthat::test_that("V209 PCA scores agree with prcomp on the median-filled display matrix, and the ellipse is the 95 per cent contour", {
  m <- as.matrix(read.delim(file.path(fixtures, "qc_pca_matrix.tsv"), check.names = FALSE, row.names = 1))
  filled <- t(apply(m, 1, function(x) { x[is.na(x)] <- stats::median(x, na.rm = TRUE); x }))
  p <- stats::prcomp(t(filled), center = TRUE, scale. = FALSE)
  src <- read_fixture("qc_pca_ellipses_source.tsv")
  for (k in 1:2) {
    col <- paste0("PC", k)
    s <- sign(sum(p$x[src$sample, k] * src[[col]]))
    testthat::expect_equal(as.numeric(s * p$x[src$sample, k]), src[[col]], tolerance = 1e-10)
  }
  d <- src[src$group == "G1", ]
  ell <- catalogue_ellipse(d$PC1, d$PC2)
  cov_m <- stats::cov(cbind(d$PC1, d$PC2))
  centre <- c(mean(d$PC1), mean(d$PC2))
  mahal <- apply(as.matrix(ell[, c("x", "y")]), 1, function(pt) { z <- pt - centre; drop(t(z) %*% solve(cov_m, z)) })
  testthat::expect_equal(mahal, rep(stats::qchisq(0.95, df = 2), length(mahal)), tolerance = 1e-8)
})

testthat::test_that("V214 exclusive region counts equal the brute-force enumeration of the 2^k - 1 patterns", {
  src <- read_fixture("diff_upset_venn_k5.tsv")
  membership <- as.matrix(src[, c("A", "B", "C", "D", "E")])
  regions <- catalogue_regions(membership)
  testthat::expect_equal(nrow(regions), 2^5 - 1)
  for (i in seq_len(nrow(regions))) {
    pattern <- as.integer(strsplit(regions$pattern[i], "")[[1]])
    want <- sum(apply(membership, 1, function(row) all(row == pattern)))
    testthat::expect_identical(regions$count[i], as.integer(want))
  }
  testthat::expect_identical(sum(regions$count), as.integer(sum(rowSums(membership) > 0)))
})

testthat::test_that("V214 negative: a Venn for k = 5 fails E_VENN_K and the UpSet still draws", {
  src <- read_fixture("diff_upset_venn_k5.tsv")
  testthat::expect_error(catalogue_upset_venn(src, c("A", "B", "C", "D", "E"), venn = TRUE), "E_VENN_K")
  testthat::expect_s3_class(catalogue_upset_venn(src, c("A", "B", "C", "D", "E"), venn = FALSE), "ggplot")
})

testthat::test_that("V213 negative: a non-DEP feature in the heatmap is refused as a set mismatch", {
  src <- read_fixture("diff_dep_heatmap_with_nondep.tsv")
  testthat::expect_error(catalogue_dep_heatmap(src), "E_FIGURE_SET_MISMATCH")
})

testthat::test_that("V210 negative: a label for a feature absent from the source is refused", {
  src <- read_fixture("diff_volcano_labelled_source.tsv")
  absent <- read_fixture("diff_volcano_labels_absent.tsv")$feature
  testthat::expect_error(catalogue_volcano(src, absent), "E_FIGURE_SOURCE_MISSING")
})

testthat::test_that("V206 negative: an intensity source without the observed mask is refused as imputed", {
  src <- read_fixture("qc_intensity_distributions_imputed.tsv")
  testthat::expect_error(catalogue_qc_intensity(src), "E_FIGURE_SOURCE_IMPUTED")
})

# ---------------------------------------------------------------- stage 2 (V215-V222)

testthat::test_that("V215 the ROC band equals pROC::ci.se with the recorded seed; a training-fold ROC is refused", {
  src <- read_fixture("bio_roc_ci_band_source.tsv")
  scores <- utils::read.delim(file.path(fixtures, "bio_roc_scores.tsv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  r <- pROC::roc(response = scores$class, predictor = scores$score, levels = c("Control", "Case"), direction = "<", quiet = TRUE)
  set.seed(src$seed[1])
  se <- pROC::ci.se(r, specificities = seq(0, 1, 0.02), boot.n = src$B[1])
  band <- src[src$component == "band", ]
  testthat::expect_equal(band$band_lower, unname(se[, 1]), tolerance = 1e-10)
  testthat::expect_equal(band$band_upper, unname(se[, 3]), tolerance = 1e-10)
  testthat::expect_s3_class(catalogue_roc_ci_band(src), "ggplot")
  testthat::expect_error(catalogue_roc_ci_band(read_fixture("bio_roc_ci_band_training_source.tsv")), "E_BIOMARKER_CLAIM")
})

testthat::test_that("V216 forest estimates and intervals equal stats::t.test; a panel without a claim label is refused", {
  src <- read_fixture("bio_panel_forest_source.tsv")
  vals <- utils::read.delim(file.path(fixtures, "bio_panel_values.tsv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  for (i in seq_len(nrow(src))) {
    x <- vals[vals$feature == src$feature[i], ]
    t <- stats::t.test(x$value[x$group == "Case"], x$value[x$group == "Control"])
    testthat::expect_equal(src$estimate[i], unname(t$estimate[1] - t$estimate[2]), tolerance = 1e-10)
    testthat::expect_equal(src$ci_lower[i], t$conf.int[1], tolerance = 1e-10)
    testthat::expect_equal(src$ci_upper[i], t$conf.int[2], tolerance = 1e-10)
  }
  testthat::expect_s3_class(catalogue_forest(src), "ggplot")
  testthat::expect_error(catalogue_forest(read_fixture("bio_panel_forest_unlabelled_source.tsv")), "E_BIOMARKER_CLAIM")
})

testthat::test_that("V217 the displayed permutation P equals (k + 1) / (B + 1); a wrong displayed P is refused", {
  src <- read_fixture("bio_cv_auc_vs_permutation_source.tsv")
  observed <- mean(src$value[src$series == "cv_repeat"])
  null <- src$value[src$series == "null"]
  k <- sum(null >= observed)
  testthat::expect_equal(src$P_display[1], (k + 1) / (src$B[1] + 1), tolerance = 1e-12)
  testthat::expect_s3_class(catalogue_cv_permutation(src), "ggplot")
  testthat::expect_error(catalogue_cv_permutation(read_fixture("bio_cv_auc_vs_permutation_wrong_p_source.tsv")), "E_FIGURE_PERMUTATION_P")
})

testthat::test_that("V218 the gene ratio is k / n, a set without a foreground hit is refused as enriched", {
  src <- read_fixture("path_enrichment_source.tsv")
  hits <- src[src$k >= 1, ]
  hits <- hits[order(hits$q_value, hits$set_id), ]          # the rows are drawn most significant first
  built <- ggplot2::ggplot_build(catalogue_enrichment(src))$data[[1]]
  testthat::expect_equal(built$x, hits$k / hits$n, tolerance = 1e-10)
  testthat::expect_error(catalogue_enrichment(src, c("SET_A", "SET_E")), "E_FIGURE_SET_MISMATCH")
})

testthat::test_that("V218 the colour follows -log10 q and the rows follow q, the most significant set at the top", {
  src <- read_fixture("path_enrichment_source.tsv")
  hits <- src[src$k >= 1, ]
  hits <- hits[order(hits$q_value, hits$set_id), ]
  p <- catalogue_enrichment(src)
  built <- ggplot2::ggplot_build(p)$data[[1]]
  # darker (lower summed RGB) as -log10 q rises: the colours are ordered by the statistic
  rgb_sum <- colSums(grDevices::col2rgb(built$fill))
  neglog <- -log10(hits$q_value)
  testthat::expect_true(all(diff(rgb_sum[order(neglog)]) <= 0))
  # the most significant set (first row) is drawn at the top: the largest factor position
  testthat::expect_equal(which.max(built$y), 1L)
  testthat::expect_equal(built$y[order(hits$q_value)], rev(sort(built$y)))
  # black outline layer: exactly the sets with q < 0.05
  outline <- ggplot2::ggplot_build(p)$data[[2]]
  testthat::expect_equal(nrow(outline), sum(hits$q_value < 0.05))
})

testthat::test_that("V221 the PCoA plot keeps its group legend, titled group", {
  src <- read_fixture("perm_pcoa_source.tsv")
  p <- catalogue_pcoa(src)
  testthat::expect_false(identical(as.character(ns$completed_theme(p)$legend.position), "none"))
  built <- ggplot2::ggplot_build(p)
  testthat::expect_equal(sort(built$plot$scales$get_scales("colour")$get_labels()), c("A", "B"))
})

testthat::test_that("V219 module cells equal the PC1 oracle; an undeclared module is refused", {
  src <- read_fixture("net_module_heatmap_source.tsv")
  m <- utils::read.delim(file.path(fixtures, "net_module_matrix.tsv"), stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")
  s <- utils::read.delim(file.path(fixtures, "net_module_samples.tsv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  for (i in seq_len(nrow(src))) {
    X <- as.matrix(m[m$module == src$module[i], s$sample])
    Z <- t(apply(X, 1, function(x) (x - mean(x)) / stats::sd(x)))
    e <- stats::prcomp(t(Z), center = FALSE, scale. = FALSE)$x[, 1]
    if (stats::cor(e, colMeans(Z)) < 0) e <- -e
    testthat::expect_equal(src$eigengene[i], mean(e[s$group == src$group[i]]), tolerance = 1e-10)
  }
  testthat::expect_s3_class(catalogue_module_heatmap(src, c("M1", "M2", "M3")), "ggplot")
  testthat::expect_error(catalogue_module_heatmap(read_fixture("net_module_heatmap_undeclared_source.tsv"), c("M1", "M2", "M3")),
                         "E_FIGURE_SET_MISMATCH")
})

testthat::test_that("V220 the layout is deterministic and an edge with an unknown endpoint is refused", {
  nodes <- data.frame(id = c("N01", "N02", "N03", "N04"), module = c("C1", "C1", "C2", "C2"), stringsAsFactors = FALSE)
  testthat::expect_identical(catalogue_network_layout(nodes), catalogue_network_layout(nodes))
  testthat::expect_equal(nrow(catalogue_network_layout(nodes)), 4L)
  src <- read_fixture("net_network_graph_source.tsv")
  testthat::expect_s3_class(catalogue_network_graph(src), "ggplot")
  bad <- src
  bad$to[bad$kind == "edge"][1] <- "N99"
  testthat::expect_error(catalogue_network_graph(bad), "E_FIGURE_SOURCE_MISSING")
})

testthat::test_that("V221 the coordinates equal cmdscale and each drawn centroid is the mean of its own group", {
  src <- read_fixture("perm_pcoa_source.tsv")
  m <- utils::read.delim(file.path(fixtures, "perm_pcoa_matrix.tsv"), row.names = 1, check.names = FALSE, fileEncoding = "UTF-8")
  X <- as.matrix(m[, setdiff(names(m), "group"), drop = FALSE])
  cm <- stats::cmdscale(stats::dist(X), k = 2, eig = TRUE)
  for (k in 1:2) {
    col <- paste0("axis", k)
    s <- sign(sum(cm$points[src$observation_id, k] * src[[col]]))
    testthat::expect_equal(as.numeric(s * cm$points[src$observation_id, k]), src[[col]], tolerance = 1e-8)
  }
  built <- ggplot2::ggplot_build(catalogue_pcoa(src))$data
  centroid <- built[[which(vapply(built, function(d) "shape" %in% names(d) && all(d$shape == 23), logical(1)))[1]]]
  means_x <- tapply(src$axis1, src$group, mean)
  means_y <- tapply(src$axis2, src$group, mean)
  testthat::expect_equal(as.numeric(centroid$x), as.numeric(means_x[c("A", "B")]), tolerance = 1e-10)
  testthat::expect_equal(as.numeric(centroid$y), as.numeric(means_y[c("A", "B")]), tolerance = 1e-10)
  testthat::expect_gt(max(abs(centroid$x - mean(src$axis1))), 1e-6)   # a full-sample centroid would not match
})

testthat::test_that("V222 distances equal betadisper; a P without its seed is refused", {
  src <- read_fixture("perm_permdisp_source.tsv")
  m <- utils::read.delim(file.path(fixtures, "perm_permdisp_matrix.tsv"), row.names = 1, check.names = FALSE, fileEncoding = "UTF-8")
  X <- as.matrix(m[, setdiff(names(m), "group"), drop = FALSE])
  bd <- vegan::betadisper(stats::dist(X), factor(m$group))
  testthat::expect_equal(src$distance_to_centroid, unname(bd$distances), tolerance = 1e-10)
  testthat::expect_s3_class(catalogue_permdisp(src), "ggplot")
  testthat::expect_error(catalogue_permdisp(read_fixture("perm_permdisp_no_seed_source.tsv")), "E_FIGURE_SOURCE_MISSING")
})

# ---------------------------------------------------------------- retrofit of existing figures (D-72)

testthat::test_that("retrofit: a scatter of a registered table has one point per row and an unknown kind is refused", {
  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(c("feature_id\teffect_log2\tneg_log10_p\tq_at_or_below_cutoff",
               "F1\t1.5\t3.2\ttrue", "F2\t-0.4\t0.8\tfalse", "F3\tNA\t2.0\tfalse"), tmp)
  src <- utils::read.delim(tmp, stringsAsFactors = FALSE, na.strings = "NA", fileEncoding = "UTF-8")
  opts <- list(x = "effect_log2", y = "neg_log10_p", flag_col = "q_at_or_below_cutoff", flag_values = list("true"),
               xlab = "effect", ylab = "-log10 P", title = "t")
  built <- ggplot2::ggplot_build(ns$retro_plot_scatter(src, opts))
  testthat::expect_equal(nrow(built$data[[1]]), 2L)                 # the row with NA effect is not drawn
  testthat::expect_equal(sort(built$data[[1]]$y), c(0.8, 3.2))
  testthat::expect_error(ns$render_retrofit_png("no_such_kind", tmp, tempfile(fileext = ".png")), "E_FIGURE_SOURCE_MISSING")
})

testthat::test_that("D1 dot plot: one jittered point per sample, a mean bar and a SD bar per group, and the top bracket label inside the y range", {
  values <- read_fixture("diff_dotplot_values.tsv")
  top <- tapply(values$value, values$protein, max)
  brackets <- data.frame(protein = names(top), group1 = "Control", group2 = "Chronic", label = "FR-198",
                         y.position = as.numeric(top) + 0.5, stringsAsFactors = FALSE)
  built <- ggplot2::ggplot_build(catalogue_dotplot_brackets(values, brackets))
  points <- built$data[[1]]
  testthat::expect_equal(nrow(points), nrow(values))            # one point per sample
  testthat::expect_true(all(!is.na(points$colour)))             # every point has a group colour
  sd_bars <- built$data[[3]]                                    # the +/- SD error bars: one per group and protein
  testthat::expect_equal(nrow(sd_bars), length(unique(values$group)) * length(unique(values$protein)))
  tip <- 0.02 * diff(range(values$value))
  for (i in seq_len(nrow(built$layout$layout))) {
    protein <- as.character(built$layout$layout$protein[i])
    label_y <- brackets$y.position[brackets$protein == protein] + tip
    testthat::expect_gte(built$layout$panel_params[[i]]$y.range[2], label_y)   # the top label is not clipped
  }
})

testthat::test_that("D2 Venn: each circle carries its contrast name and total, the wrapped name is rejoined, and no A/B letter is drawn", {
  src <- read_fixture("diff_upset_venn_k3.tsv")
  labels <- c("disease vs control", "treated vs disease", "treated vs control")
  text <- gsub("\\s+", " ", catalogue_text_labels(catalogue_venn(src, c("A", "B", "C"), set_labels = labels)))
  counts <- as.integer(colSums(src[, c("A", "B", "C")]))
  for (i in seq_along(labels)) {
    testthat::expect_true(any(grepl(sprintf("%s (n = %d)", labels[i], counts[i]), text, fixed = TRUE)))
  }
  testthat::expect_false(any(text %in% c("A", "B", "C")))
})
