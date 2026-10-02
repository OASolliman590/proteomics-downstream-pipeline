# R03 QC acceptance tests (V026).  Independent oracle: base-R prcomp/SVD after
# the declared display fill and centering, and cor() on shared observed cells.

ns <- asNamespace("proteomicsCore")
fn <- function(name) get(name, envir = ns)

testthat::test_that("V026 PCA matches an independent SVD up to sign", {
  x <- matrix(c(1, 2, 3, 4,  2, 4, 5, 9,  0, 1, 0, 2,  3, 3, 4, 1), 4, byrow = TRUE, dimnames = list(paste0("F", 1:4), paste0("O", 1:4)))
  pca <- fn("qc_pca")(x)
  testthat::expect_identical(pca$state, "COMPLETED")
  reference <- stats::prcomp(t(x), center = TRUE, scale. = FALSE)
  k <- ncol(pca$scores)
  testthat::expect_equal(pca$variance$eigenvalue, reference$sdev[seq_len(k)]^2, tolerance = 1e-10)
  testthat::expect_equal(sum(pca$variance$variance_explained), 1, tolerance = 1e-12)
  for (j in seq_len(k)) testthat::expect_equal(abs(unname(pca$scores[, j])), abs(unname(reference$x[, j])), tolerance = 1e-10)
})

testthat::test_that("V026 one-NA variant uses the declared median display fill", {
  x <- matrix(c(1, 2, 3, 4,  2, 4, 5, 9,  0, 1, 0, 2,  3, 3, 4, 1), 4, byrow = TRUE, dimnames = list(paste0("F", 1:4), paste0("O", 1:4)))
  x[2, 3] <- NA
  filled <- x; filled[2, 3] <- stats::median(c(2, 4, 9))
  pca <- fn("qc_pca")(x)
  reference <- stats::prcomp(t(filled), center = TRUE, scale. = FALSE)
  testthat::expect_equal(pca$variance$eigenvalue, reference$sdev[seq_len(nrow(pca$variance))]^2, tolerance = 1e-10)
  testthat::expect_true(is.na(x[2, 3]))             # caller's matrix never filled
})

testthat::test_that("V026 constant and empty inputs give explanatory states, no fake variance", {
  constant <- matrix(5, 3, 4, dimnames = list(paste0("F", 1:3), paste0("O", 1:4)))
  p <- fn("qc_pca")(constant)
  testthat::expect_identical(p$state, "INAPPLICABLE"); testthat::expect_identical(p$reason_code, "E_QC_CONSTANT"); testthat::expect_null(p$variance)
  empty <- matrix(NA_real_, 2, 4, dimnames = list(c("F1", "F2"), paste0("O", 1:4)))
  e <- fn("qc_pca")(empty)
  testthat::expect_identical(e$reason_code, "E_QC_NO_FEATURES")
})

testthat::test_that("V026 correlations use shared observed cells and retain n", {
  x <- matrix(c(1, 2, 3, 4, 5,  2, 4, 6, 8, NA,  5, 3, 2, 1, 0), 5, dimnames = list(paste0("F", 1:5), c("A", "B", "C")))
  usable <- !is.na(x); usable[1, "C"] <- FALSE
  r <- fn("qc_correlations")(x, usable)
  ab <- r[r$observation_a == "A" & r$observation_b == "B", ]
  testthat::expect_equal(ab$n_shared, 4L); testthat::expect_equal(ab$pearson_r, stats::cor(x[1:4, "A"], x[1:4, "B"]))
  ac <- r[r$observation_a == "A" & r$observation_b == "C", ]
  testthat::expect_equal(ac$n_shared, 4L); testthat::expect_equal(ac$pearson_r, stats::cor(x[2:5, "A"], x[2:5, "C"]))
})
