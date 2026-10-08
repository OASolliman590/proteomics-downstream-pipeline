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
