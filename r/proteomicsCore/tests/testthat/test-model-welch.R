# Welch engine at the R level (packet R15a): V168, V171, V172 and V173.  Each oracle is a direct
# stats::t.test() or limma call written in this file; each negative case is stated next to its oracle.

ns <- asNamespace("proteomicsCore")
fn <- function(name) get(name, envir = ns)

welch_fixture <- function(seed, n_num, n_den, features) {
  set.seed(seed)
  sd_num <- 0.4; sd_den <- 1.6
  x <- matrix(0, features, n_num); y <- matrix(0, features, n_den)
  for (f in seq_len(features)) {
    x[f, ] <- 10 + 0.02 * f + rnorm(n_num, 0, sd_num)
    y[f, ] <- 10 + 0.02 * f + (if (f <= 3) 1.5 else 0) + rnorm(n_den, 0, sd_den)
  }
  list(x = x, y = y)
}

testthat::test_that("V168 Welch per-contrast statistics equal stats::t.test(var.equal = FALSE) within 1e-10", {
  fx <- welch_fixture(168, 6, 9, 200)
  for (f in seq_len(200)) {
    row <- fn("welch_two_sample_row")(fx$x[f, ], fx$y[f, ], 0.95)
    oracle <- stats::t.test(fx$x[f, ], fx$y[f, ], var.equal = FALSE, conf.level = 0.95)
    testthat::expect_lt(abs(row[["t"]] - unname(oracle$statistic)), 1e-10)
    testthat::expect_lt(abs(row[["df"]] - unname(oracle$parameter)), 1e-10)
    testthat::expect_lt(abs(row[["p"]] - oracle$p.value), 1e-10)
    testthat::expect_lt(abs(row[["effect"]] - (mean(fx$x[f, ]) - mean(fx$y[f, ]))), 1e-10)
    testthat::expect_lt(abs(row[["lower"]] - oracle$conf.int[1]), 1e-10)
    testthat::expect_lt(abs(row[["upper"]] - oracle$conf.int[2]), 1e-10)
  }
})

testthat::test_that("V168 negative: the pooled-variance form and a reversed effect sign are different outputs", {
  fx <- welch_fixture(168, 6, 9, 200)
  row <- fn("welch_two_sample_row")(fx$x[10, ], fx$y[10, ], 0.95)
  pooled <- stats::t.test(fx$x[10, ], fx$y[10, ], var.equal = TRUE)
  testthat::expect_false(isTRUE(all.equal(unname(pooled$parameter), row[["df"]])))   # pooled df = n - 2 = 13
  testthat::expect_false(isTRUE(all.equal(row[["effect"]], -(mean(fx$x[10, ]) - mean(fx$y[10, ])))))
  reversed <- fn("welch_two_sample_row")(fx$y[10, ], fx$x[10, ], 0.95)
  testthat::expect_equal(reversed[["effect"]], -row[["effect"]], tolerance = 1e-12)   # the sign follows the numerator
  testthat::expect_equal(reversed[["p"]], row[["p"]], tolerance = 1e-12)
})

testthat::test_that("V171 paired variant equals stats::t.test(paired = TRUE) with df = n - 1 = 7", {
  set.seed(171)
  x <- rnorm(8, 10, 1); y <- x + rnorm(8, 0.6, 0.4)
  row <- fn("welch_paired_row")(x, y, 0.95)
  oracle <- stats::t.test(x, y, paired = TRUE, conf.level = 0.95)
  testthat::expect_equal(unname(row[["df"]]), 7)
  testthat::expect_lt(abs(row[["t"]] - unname(oracle$statistic)), 1e-10)
  testthat::expect_lt(abs(row[["effect"]] - unname(oracle$estimate)), 1e-10)
  testthat::expect_lt(abs(row[["p"]] - oracle$p.value), 1e-10)
  testthat::expect_lt(abs(row[["lower"]] - oracle$conf.int[1]), 1e-10)
  testthat::expect_lt(abs(row[["upper"]] - oracle$conf.int[2]), 1e-10)
})

testthat::test_that("V171 negative: a two-sample Welch form on a declared paired design is refused (E_WELCH_PAIRED_UNIT)", {
  obs <- data.frame(observation_id = c(paste0("S", 1:8, "A"), paste0("S", 1:8, "B")), group = rep(c("A", "B"), each = 8),
                    subject_id = rep(paste0("subj", 1:8), 2), stringsAsFactors = FALSE)
  blocked <- list(model_id = "m", blocked = TRUE, subject_column = "subject_id", paired = "false")
  testthat::expect_error(fn("welch_variant")(blocked, obs, obs$group, "B", "A"), "E_WELCH_PAIRED_UNIT")
  declared_true_unblocked <- list(model_id = "m", blocked = FALSE, paired = "true")
  testthat::expect_error(fn("welch_variant")(declared_true_unblocked, obs, obs$group, "B", "A"), "E_WELCH_PAIRED_UNIT")
  auto <- fn("welch_variant")(list(model_id = "m", blocked = TRUE, subject_column = "subject_id", paired = "auto"), obs, obs$group, "B", "A")
  testthat::expect_identical(auto$variant, "paired")
  testthat::expect_length(auto$pairs$num, 8L)
})

testthat::test_that("V172 covariate adaptation runs the unadjusted two-sample Welch, not the sex-adjusted fit", {
  # Sex is unbalanced across the two groups, so the sex-adjusted group effect differs from the unadjusted one.
  sex_x <- c("F", "F", "F", "F", "M", "M"); sex_y <- c("M", "M", "M", "F", "F", "F")
  x <- 10 + 0.5 * (sex_x == "M") + c(0.1, -0.2, 0.3, -0.1, 0.2, -0.3)
  y <- 11 + 0.5 * (sex_y == "M") + c(-0.1, 0.2, -0.3, 0.1, -0.2, 0.3)
  row <- fn("welch_two_sample_row")(x, y, 0.95)
  oracle <- stats::t.test(x, y, var.equal = FALSE)
  testthat::expect_lt(abs(row[["t"]] - unname(oracle$statistic)), 1e-10)
  testthat::expect_lt(abs(row[["df"]] - unname(oracle$parameter)), 1e-10)
  testthat::expect_lt(abs(row[["effect"]] - unname(oracle$estimate[1] - oracle$estimate[2])), 1e-10)
  # negative: a silently sex-adjusted effect would differ from the unadjusted Welch effect (orientation: group x minus group y)
  grp <- rep(c(0, 1), each = 6)
  adjusted <- stats::coef(stats::lm(c(x, y) ~ grp + c(sex_x, sex_y)))["grp"]
  testthat::expect_false(isTRUE(all.equal(unname(-adjusted), row[["effect"]])))
})

testthat::test_that("V173 Welch and limma give identical effects and distinct P values on the same heterogeneous fixture", {
  testthat::skip_if_not_installed("limma")
  fx <- welch_fixture(173, 12, 12, 200)
  Y <- cbind(fx$x, fx$y); colnames(Y) <- paste0("o", seq_len(24))
  group <- factor(rep(c("A", "B"), each = 12), levels = c("A", "B"))
  design <- model.matrix(~ 0 + group); colnames(design) <- c("A", "B")
  fit <- limma::lmFit(Y, design)
  effect_oracle <- fit$coefficients[, "B"] - fit$coefficients[, "A"]   # lmFit with ~ group and no moderation
  welch_effect <- vapply(seq_len(200), function(f) fn("welch_two_sample_row")(Y[f, 13:24], Y[f, 1:12], 0.95)[["effect"]], numeric(1))
  testthat::expect_lt(max(abs(welch_effect - (rowMeans(Y[, 13:24]) - rowMeans(Y[, 1:12])))), 1e-10)
  testthat::expect_lt(max(abs(welch_effect - effect_oracle)), 1e-10)
  # limma moderated P (the frozen settings: trend and robust) from the same lmFit
  cm <- limma::makeContrasts(B - A, levels = design)
  moderated <- limma::eBayes(limma::contrasts.fit(fit, cm), trend = TRUE, robust = TRUE)
  limma_p <- moderated$p.value[, 1]
  welch_p <- vapply(seq_len(200), function(f) fn("welch_two_sample_row")(Y[f, 13:24], Y[f, 1:12], 0.95)[["p"]], numeric(1))
  welch_oracle_p <- vapply(seq_len(200), function(f) stats::t.test(Y[f, 13:24], Y[f, 1:12], var.equal = FALSE)$p.value, numeric(1))
  testthat::expect_lt(max(abs(welch_p - welch_oracle_p)), 1e-10)
  testthat::expect_gt(max(abs(welch_p - limma_p)), 1e-6)   # negative: identical P values would mean limma pooled the variances
})
