# R14d biomarker acceptance at the R level (V150, V151, V152, V155, V156, V158).  Oracles are pROC, direct glmnet and
# e1071 fits and hand computations written here.

ns <- asNamespace("proteomicsCore")
fn <- function(name) get(name, envir = ns)
skip_bm <- function() { testthat::skip_if_not_installed("pROC"); testthat::skip_if_not_installed("glmnet"); testthat::skip_if_not_installed("e1071") }

testthat::test_that("V150 DeLong and stratified bootstrap CIs equal pROC; AUC is the Mann-Whitney statistic", {
  skip_bm()
  set.seed(2); pos <- rnorm(15, 1); neg <- rnorm(15); neg[3] <- pos[2]          # one tie
  r <- pROC::roc(controls = neg, cases = pos, direction = "<", quiet = TRUE)
  d <- fn("pd_delong")(pos, neg, 0.95); ci <- pROC::ci.auc(r, method = "delong")
  testthat::expect_equal(d$auc, as.numeric(r$auc), tolerance = 1e-12)
  testthat::expect_equal(c(d$lower, d$upper), c(ci[1], ci[3]), tolerance = 1e-10)
  testthat::expect_equal(fn("pd_auc")(c(pos, neg), rep(1:0, each = 15)), as.numeric(r$auc), tolerance = 1e-12)
  b <- fn("pd_boot_auc")(pos, neg, 300L, 9L)
  set.seed(9, kind = "L'Ecuyer-CMRG"); cb <- pROC::ci.auc(r, method = "bootstrap", boot.n = 300, boot.stratified = TRUE)
  testthat::expect_equal(unname(c(b$lower, b$upper)), unname(c(cb[1], cb[3])), tolerance = 1e-12)
})

testthat::test_that("V151 folds are stratified and grouped: every fold holds both classes when each class has >= k units", {
  set.seed(4); cls <- c(rep("1", 7), rep("0", 12))
  f <- fn("pd_bm_folds")(cls, 5L)
  testthat::expect_true(all(vapply(1:5, function(k) all(c("0", "1") %in% cls[f == k]), logical(1))))
  testthat::expect_true(max(table(f)) - min(table(f)) <= 1L)
})

testthat::test_that("V155 scores are oriented from the fitted model's class coding for glmnet and SVM", {
  skip_bm()
  set.seed(6); X <- matrix(rnorm(80), 40, 2, dimnames = list(NULL, c("a", "b"))); y <- rep(1:0, each = 20); X[y == 1, 1] <- X[y == 1, 1] + 3
  lab <- ifelse(y == 1, "P", "N")
  for (levels in list(c("P", "N"), c("N", "P"))) {
    m <- fn(".bm_fit_classifier")(X, lab, levels, "P", "penalized_logistic", 0.1)
    s <- fn(".bm_predict")(m, X)$score[, 1]
    direct <- glmnet::glmnet(X, factor(lab, levels = levels), family = "binomial", alpha = 0, lambda = 0.1, standardize = FALSE)
    p <- as.numeric(stats::predict(direct, X, s = 0.1, type = "response")); oracle <- if (direct$classnames[2] == "P") p else 1 - p
    testthat::expect_equal(s, oracle, tolerance = 1e-12)
    testthat::expect_gt(fn("pd_auc")(s, y), 0.95)
    sv <- fn(".bm_fit_classifier")(X, lab, levels, "P", "svm_linear", 1)
    ds <- fn(".bm_predict")(sv, X)$score[, 1]
    e <- e1071::svm(X, factor(lab, levels = levels), kernel = "linear", cost = 1, scale = FALSE)
    dv <- attr(stats::predict(e, X, decision.values = TRUE), "decision.values"); orient <- if (strsplit(colnames(dv)[1], "/")[[1]][1] == "P") 1 else -1
    testthat::expect_equal(ds, as.numeric(dv[, 1]) * orient, tolerance = 1e-10)
    testthat::expect_gt(fn("pd_auc")(ds, y), 0.95)                           # never 1 - AUC
  }
})

testthat::test_that("V152 the transform audit refuses a transform fitted on held-out samples", {
  ok <- list(list(repeat_id = 1L, fold = 2L, transform = "scaling", fitted_on = c("a", "b"), applied_to = "c"))
  testthat::expect_silent(fn("pd_bm_audit_check")(ok, list(`1 2` = "c")))
  bad <- list(list(repeat_id = 1L, fold = 2L, transform = "selection", fitted_on = c("a", "c"), applied_to = "c"))
  testthat::expect_error(fn("pd_bm_audit_check")(bad, list(`1 2` = "c")), "E_BIOMARKER_LEAKAGE")
})

testthat::test_that("V156 label permutations respect units: whole subjects or within-subject swaps", {
  units <- rep(paste0("s", 1:5), each = 2)
  set.seed(1, kind = "L'Ecuyer-CMRG")
  between <- rep(c(1L, 1L, 0L, 0L, 1L), each = 2)
  for (i in 1:30) { p <- fn("pd_bm_permute_labels")(between, units); testthat::expect_true(all(tapply(p, units, function(v) length(unique(v))) == 1L)); testthat::expect_identical(sum(p), sum(between)) }
  paired <- rep(c(1L, 0L), 5)
  for (i in 1:30) { p <- fn("pd_bm_permute_labels")(paired, units); testthat::expect_true(all(tapply(p, units, sum) == 1L)) }
})

testthat::test_that("V158 the Youden threshold maximises sensitivity + specificity on the given (training) scores", {
  s <- c(0.1, 0.3, 0.35, 0.6, 0.8, 0.9); y <- c(0L, 0L, 1L, 0L, 1L, 1L)
  t <- fn("pd_youden")(s, y)
  J <- function(t) mean(s[y == 1] >= t) + mean(s[y == 0] < t) - 1
  testthat::expect_equal(J(t), max(vapply(c(0, s, 1), J, numeric(1))))
  testthat::expect_error(fn(".bm_threshold")("youden_test", s, y), "E_BIOMARKER_THRESHOLD_LEAKAGE")
})
