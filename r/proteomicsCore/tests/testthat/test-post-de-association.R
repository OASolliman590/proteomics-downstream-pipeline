# R14c association acceptance at the R level (V145, V147).  Oracles are cor(), cor.test(), lm() residuals and a
# seeded sample.int() loop written here.

ns <- asNamespace("proteomicsCore")
fn <- function(name) get(name, envir = ns)

testthat::test_that("V145 correlations equal cor()/cor.test() and the permutation P is (k+1)/(B+1) over units", {
  set.seed(5); x <- rnorm(8); age <- rnorm(8); y <- 0.7 * x + 0.3 * age + rnorm(8, 0, 0.5)
  perms <- fn("pd_unit_permutations")(8L, 99L, 11L)
  set.seed(11, kind = "L'Ecuyer-CMRG"); oracle_perms <- lapply(1:99, function(b) sample.int(8))
  testthat::expect_identical(perms, oracle_perms)
  for (method in c("pearson", "spearman")) {
    res <- fn("pd_correlation")(y, x, method, NULL, perms)
    testthat::expect_equal(res$r, unname(stats::cor.test(y, x, method = method, exact = FALSE)$estimate), tolerance = 1e-12)
    k <- sum(abs(vapply(oracle_perms, function(p) stats::cor(y, x[p], method = method), numeric(1))) >= abs(res$r) - 1e-12)
    testthat::expect_identical(res$k, k); testthat::expect_equal(res$p, (k + 1) / 100)
  }
  part <- fn("pd_correlation")(y, x, "partial", cbind(1, age), perms)
  testthat::expect_equal(part$r, stats::cor(stats::residuals(stats::lm(y ~ age)), stats::residuals(stats::lm(x ~ age))), tolerance = 1e-12)
})

testthat::test_that("V145 negative: subject-blocked permutations keep subjects whole or permute within subjects only", {
  subjects <- rep(c("s1", "s2", "s3", "s4"), each = 2)
  within <- fn("pd_unit_permutations")(8L, 50L, 3L, subjects)
  testthat::expect_true(all(vapply(within, function(p) all(subjects[p] == subjects), logical(1))))
  whole <- fn("pd_subject_permutations")(subjects, 50L, 3L)
  testthat::expect_true(all(vapply(whole, function(p) all(tapply(subjects[p], subjects, function(v) length(unique(v))) == 1L), logical(1))))
  testthat::expect_true(any(vapply(whole, function(p) any(subjects[p] != subjects), logical(1))))
})

testthat::test_that("V147 within-group slopes equal lm() per group (Simpson check input)", {
  set.seed(8); g <- rep(c("A", "B"), each = 6); x <- c(runif(6, 1, 7), runif(6, 4, 10)); y <- 6 * (g == "B") - 0.4 * x + rnorm(12, 0, 0.1)
  sl <- fn(".pd_within_slopes")(y, x, g)
  testthat::expect_equal(unname(sl["A"]), unname(stats::coef(stats::lm(y[g == "A"] ~ x[g == "A"]))[2]), tolerance = 1e-12)
  testthat::expect_true(all(sl < 0)); testthat::expect_gt(unname(stats::coef(stats::lm(y ~ x))[2]), 0)
})
