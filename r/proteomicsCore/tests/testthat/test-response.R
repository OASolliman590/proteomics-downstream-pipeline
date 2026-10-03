# R09 acceptance at the R level (V083-V088, V090).  Oracles are literal arithmetic.
ns <- asNamespace("proteomicsCore"); fn <- function(name) get(name, envir = ns)

testthat::test_that("V083 bounded response classes with literal boundaries", {
  ri <- c(-0.1, 0, 0.299999, 0.3, 0.799999, 0.8, 1, 1.2, 1.200001, 3)
  expected <- c("same_direction", "limited_return", "limited_return", "partial_return", "partial_return", "near_restoration", "near_restoration", "near_restoration", "overshoot", "overshoot")
  testthat::expect_identical(vapply(ri, fn("response_class"), ""), expected)
  computed <- fn("reversal_index")(rep(1, length(ri)), -ri, 0.25)
  testthat::expect_equal(computed, ri)
  testthat::expect_true(is.na(fn("reversal_index")(0, -1, 0.25))); testthat::expect_true(is.na(fn("reversal_index")(0.1, -1, 0.25)))
  testthat::expect_identical(fn("response_class")(NA_real_), "not_evaluable")
  testthat::expect_equal(1 * (1 - 3), -2)   # residual d*(1-RI) for RI = 3 crosses control and exceeds disease
})

testthat::test_that("V084 Fieller sets match the quadratic-inequality oracle, including unbounded/disjoint", {
  oracle <- function(a, b, va, vb, cab, q) {
    A <- b^2 - q^2 * vb; B <- -2 * (a * b - q^2 * cab); C <- a^2 - q^2 * va
    sort(Re(polyroot(c(C, B, A))))
  }
  q <- stats::qt(0.975, 20)
  f <- fn("fieller_interval")(1, 2, 0.04, 0.04, -0.01, q)
  testthat::expect_identical(f$kind, "bounded"); testthat::expect_equal(c(f$lower, f$upper), oracle(1, 2, 0.04, 0.04, -0.01, q), tolerance = 1e-10)
  g <- fn("fieller_interval")(1, 0.2, 0.04, 0.09, 0, q)
  testthat::expect_true(g$kind %in% c("disjoint", "unbounded")); testthat::expect_true(any(is.infinite(unlist(g$intervals))))
})

testthat::test_that("V085 TOST passes only for the narrow interval", {
  narrow <- fn("tost")(0, 0.05, 20, 0.2, 0.05); wide <- fn("tost")(0, 1, 20, 0.2, 0.05)
  testthat::expect_equal(narrow$p_lower, stats::pt(0.2 / 0.05, 20, lower.tail = FALSE)); testthat::expect_equal(narrow$p_upper, stats::pt(-0.2 / 0.05, 20))
  testthat::expect_lt(narrow$p_value, 0.05); testthat::expect_gt(wide$p_value, 0.05)
  testthat::expect_equal(narrow$ci_upper, stats::qt(0.95, 20) * 0.05); testthat::expect_equal(narrow$ci_level, 0.9)
  testthat::expect_error(fn("tost")(0, 1, 20, 0, 0.05), "E_EQUIVALENCE_INELIGIBLE")
})

testthat::test_that("V086 conjunction joint P is the intersection-union maximum", {
  testthat::expect_equal(fn("intersection_union")(rbind(c(0.01, 0.02, 0.03, 0.04))), 0.04)
  comp <- fn("rescue_components")(1, 2, 0.2, -1.8, 0.2, 0.2, 0.15, 30, 0, 0, 0.6)
  testthat::expect_equal(comp$p_value, max(comp$p_disease, comp$p_treatment, comp$p_residual_lower, comp$p_residual_upper))
})

testthat::test_that("V088/V090 exact randomization denominators and Monte Carlo protocol", {
  s <- 1:8; a <- s > 4
  e <- fn("exact_label_randomization")(s, a)
  testthat::expect_identical(c(e$k, e$N), c(2L, 70L)); testthat::expect_equal(e$p_value, 2 / 70); testthat::expect_identical(e$observed_allocation_count, 1L)
  s2 <- c(1, 2, 4, 8, 16, 32, 64, 128); a2 <- s2 %in% c(4, 32, 64, 128)
  all <- utils::combn(8, 4); obs <- mean(s2[a2]) - mean(s2[!a2])
  manual <- sum(apply(all, 2, function(i) { x <- logical(8); x[i] <- TRUE; abs(mean(s2[x]) - mean(s2[!x])) >= abs(obs) - (1e-12 + 1e-10 * abs(obs)) }))
  e2 <- fn("exact_label_randomization")(s2, a2)
  testthat::expect_identical(e2$k, manual); testthat::expect_equal(e2$p_value, 6 / 70)
  mc <- fn("monte_carlo_label_randomization")(s, a, 999, 11)
  testthat::expect_equal(mc$p_value, (mc$k + 1) / 1000); testthat::expect_identical(mc$forced_observed, 1L)
  testthat::expect_identical(fn("monte_carlo_label_randomization")(s, a, 999, 11)$k, mc$k)
  sf <- fn("exact_sign_flip")(c(1, 2, 3, 4))
  testthat::expect_identical(sf$N, 16L); testthat::expect_equal(sf$p_value, 2 / 16)
})

testthat::test_that("V082/V087 eligibility: overlap or unknown provenance degrades to descriptive", {
  m <- list(training_subject_ids = list("t1", "t2"), selection_subject_ids = list("t3"))
  testthat::expect_true(fn("score_eligibility")(m, c("s1", "s2"), c("u1"), "independent_exact")$inferential)
  testthat::expect_identical(fn("score_eligibility")(m, c("s1", "t2"), c("u1"), "independent_exact")$reason, "E_SCORE_SELECTION_OVERLAP")
  testthat::expect_identical(fn("score_eligibility")(list(), c("s1"), "u1", "independent_exact")$reason, "E_SCORE_INDEPENDENCE_UNVERIFIED")
  testthat::expect_identical(fn("score_eligibility")(m, c("s1"), "u1", "off")$reason, "score_test_off")
  Y <- rbind(F1 = c(2, 4), F2 = c(1, NA)); colnames(Y) <- c("o1", "o2")
  man <- list(features = list(list(feature_id = "F1", weight = 1, center = 1, scale = 2), list(feature_id = "F2", weight = -1, center = 0, scale = 1)), missing_feature_policy = "require_all")
  sc <- fn("apply_fixed_score")(Y, man)
  testthat::expect_equal(sc$value, c((2 - 1) / 2 - 1, NA)); testthat::expect_equal(sc$observed_feature_fraction, c(1, 0.5))
})
