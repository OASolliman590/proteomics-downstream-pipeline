# R04 acceptance tests (V035 algebra, V037 exact contrast covariance, V038 availability).
# The oracle is a direct featurewise GLS written here, independent of the adapter.

ns <- asNamespace("proteomicsCore")
fn <- function(name) get(name, envir = ns)

gls_oracle <- function(y, X, contrast, w = rep(1, length(y)), block = NULL, rho = 0) {
  use <- !is.na(y); y <- y[use]; X <- X[use, , drop = FALSE]; w <- w[use]
  n <- length(y); V <- diag(n)
  if (!is.null(block)) { b <- block[use]; same <- outer(b, b, "=="); V[same] <- rho; diag(V) <- 1 }
  Sigma_inv <- diag(sqrt(w), n) %*% solve(V) %*% diag(sqrt(w), n)
  info <- t(X) %*% Sigma_inv %*% X
  beta <- solve(info, t(X) %*% Sigma_inv %*% y)
  r <- y - X %*% beta
  list(effect = sum(contrast * beta), unscaled = sqrt(as.numeric(t(contrast) %*% solve(info) %*% contrast)),
       sigma = sqrt(as.numeric(t(r) %*% Sigma_inv %*% r) / (n - ncol(X))))
}

fixture <- function() {
  group <- c(rep(0, 5), rep(1, 5))
  x <- c(1.0, 2.5, 3.1, 4.0, 2.2, 3.5, 5.1, 4.4, 6.0, 5.2)          # nonorthogonal covariate correlated with group
  X <- cbind(group.A = 1 - group, group.B = group, continuous.x = x - mean(x))
  Y <- rbind(F1 = 10 + 0.8 * group + 0.3 * x + c(0.1, -0.2, 0.05, 0.12, -0.07, 0.2, -0.1, 0.03, -0.15, 0.08),
             F2 = 12 - 0.5 * group + 0.1 * x + c(-0.05, 0.11, -0.2, 0.07, 0.02, -0.09, 0.14, -0.03, 0.06, -0.12),
             F3 = 8 + 0.2 * x + c(0.3, -0.1, 0.2, -0.25, 0.05, 0.1, -0.2, 0.15, -0.05, -0.2))
  Y["F2", 7] <- NA                                                   # feature-specific missing row
  list(X = X, Y = Y, w = c(1, 2, 0.5, 1.5, 1, 3, 0.7, 1, 2, 1), contrast = c(-1, 1, 0))
}

testthat::test_that("V037 contrast-as-coefficient refit matches direct GLS (weights + missing row)", {
  testthat::skip_if_not_installed("limma")
  f <- fixture()
  ex <- fn("exact_contrast_fit")(f$Y, f$X, f$contrast, weights = f$w)
  testthat::expect_identical(ex$path, "contrast_as_coefficient_refit")
  for (i in seq_len(nrow(f$Y))) {
    o <- gls_oracle(f$Y[i, ], f$X, f$contrast, f$w)
    testthat::expect_equal(unname(ex$fit$coefficients[i, 1]), o$effect, tolerance = 1e-10)
    testthat::expect_equal(unname(ex$fit$stdev.unscaled[i, 1]), o$unscaled, tolerance = 1e-10)
    testthat::expect_equal(unname(ex$fit$sigma[i]), o$sigma, tolerance = 1e-10)
  }
  # the contrasts.fit shortcut is only approximate for the sparse feature
  approx <- fn("approximate_contrast_fit")(f$Y, f$X, matrix(f$contrast, ncol = 1), weights = f$w)
  testthat::expect_gt(abs(unname(approx$fit$stdev.unscaled["F2", 1]) - gls_oracle(f$Y["F2", ], f$X, f$contrast, f$w)$unscaled), 1e-6)
})

testthat::test_that("V037 blocked variant matches direct GLS with the frozen correlation", {
  testthat::skip_if_not_installed("limma")
  f <- fixture(); block <- rep(1:5, times = 2); rho <- 0.4
  ex <- fn("exact_contrast_fit")(f$Y, f$X, f$contrast, weights = f$w, block = block, correlation = rho)
  for (i in seq_len(nrow(f$Y))) {
    o <- gls_oracle(f$Y[i, ], f$X, f$contrast, f$w, block, rho)
    testthat::expect_equal(unname(ex$fit$coefficients[i, 1]), o$effect, tolerance = 1e-8)
    testthat::expect_equal(unname(ex$fit$stdev.unscaled[i, 1]), o$unscaled, tolerance = 1e-8)
  }
  testthat::expect_error(fn("exact_contrast_fit")(f$Y, f$X, f$contrast, block = block), "E_BLOCKING_CONFLICT")
})

testthat::test_that("V035/V050 contrast covariance algebra: Cov(d,t) = -Var(U)", {
  X <- cbind(C = rep(c(1, 0, 0), each = 4), U = rep(c(0, 1, 0), each = 4), T = rep(c(0, 0, 1), each = 4))
  C <- cbind(d = c(-1, 1, 0), t = c(0, -1, 1), r = c(-1, 0, 1))
  cov <- fn("contrast_unscaled_covariance")(X, C)
  # group-mean covariance diag(.25,.25,.25) for n=4 per group
  testthat::expect_equal(cov["d", "t"], -0.25); testthat::expect_equal(cov["r", "r"], 0.5)
  testthat::expect_equal(cov["d", "d"] + cov["t", "t"] + 2 * cov["d", "t"], cov["r", "r"])
})

testthat::test_that("V036 contrast estimability uses the observed-row space", {
  X <- cbind(A = c(1, 1, 1, 0, 0, 0), B = c(0, 0, 0, 1, 1, 1))
  testthat::expect_true(fn("contrast_estimable")(X, c(-1, 1)))
  testthat::expect_false(fn("contrast_estimable")(X[1:3, , drop = FALSE], c(-1, 1)))
  testthat::expect_true(fn("contrast_estimable")(X[1:3, , drop = FALSE], c(1, 0)))
})

testthat::test_that("V038 absent adapter software is NOT_RUN, not inapplicable", {
  missing <- fn("operational_availability")("deqms", installed = function(pkg) FALSE)
  testthat::expect_identical(missing$state, "NOT_RUN")
  testthat::expect_true(missing$reason_code %in% c("E_CAPABILITY_NOT_IMPLEMENTED", "E_ENGINE_NOT_AVAILABLE"))
  limma <- fn("operational_availability")("limma", installed = function(pkg) TRUE)
  testthat::expect_identical(limma$state, "AVAILABLE")
})
