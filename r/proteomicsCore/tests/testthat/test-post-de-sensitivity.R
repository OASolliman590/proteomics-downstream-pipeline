# R14b sensitivity acceptance at the R level (V138, V140, V141).  Oracles are direct limma fits, base R set
# operations, cor()/lm() and a seeded sample() loop written here.

ns <- asNamespace("proteomicsCore")
fn <- function(name) get(name, envir = ns)

sim <- function(seed = 3, n = 6, p = 20) {
  set.seed(seed)
  g <- rep(c("A", "B"), each = n); sex <- rep(c("F", "M"), n)
  Y <- sapply(seq_along(g), function(j) 10 + rnorm(p, 0, 0.3) + ifelse(g[j] == "B" & seq_len(p) <= 5, 1, 0) + ifelse(sex[j] == "M" & seq_len(p) <= 8, 0.8, 0))
  dimnames(Y) <- list(sprintf("F%02d", seq_len(p)), paste0(g, seq_along(g)))
  obs <- data.frame(observation_id = colnames(Y), biological_unit_id = colnames(Y), group = g, sex = sex, stringsAsFactors = FALSE)
  list(Y = Y, obs = obs, g = g, sex = sex)
}

testthat::test_that("V138 pd_fit_model equals a direct limma fit with ~ group + covariate and BH within its own family", {
  testthat::skip_if_not_installed("limma")
  d <- sim()
  X <- cbind("(Intercept)" = 1, "group.B" = as.numeric(d$g == "B"), "factor.sex.M" = as.numeric(d$sex == "M")); rownames(X) <- colnames(d$Y)
  contrasts <- list(list(contrast_id = "B-A", required_groups = list("B", "A"), weights = list(0, 1, 0), role = "primary"))
  observed <- matrix(TRUE, nrow(d$Y), ncol(d$Y), dimnames = dimnames(d$Y))
  settings <- list(group_column = "group", subject_column = NULL, blocking_mode = "none", trend = TRUE, robust = TRUE, ci_level = 0.95)
  fam <- list(list(family_id = "primary__s", hypothesis_type = "protein_zero_null", role = "sensitivity", adjustment = "BH", q_cutoff = 0.05,
                   dependence_assumption = "BH", members = list(list(model_id = "s", contrast_id = "B-A"))))
  res <- fn("pd_fit_model")("s", d$Y, observed, d$obs, X, contrasts, settings, fam, list(run_id = "r", plan_hash = "h"),
                            list(policy = "available_case", minimum_observed_per_group = 2, minimum_fraction = 0.5))
  eb <- limma::eBayes(limma::lmFit(d$Y, model.matrix(~ factor(d$g) + factor(d$sex))), trend = TRUE, robust = TRUE)
  r <- res$rows[match(rownames(d$Y), res$rows$feature_id), ]
  testthat::expect_equal(r$effect, unname(eb$coefficients[, 2]), tolerance = 1e-10)
  testthat::expect_equal(r$statistic, unname(eb$t[, 2]), tolerance = 1e-8)
  testthat::expect_equal(unname(r$q_value), unname(stats::p.adjust(eb$p.value[, 2], "BH")), tolerance = 1e-10)
  testthat::expect_identical(unique(r$family_id), "primary__s")
})

testthat::test_that("V140 comparison equals set operations, cor(), lm(y ~ 0 + x) and Jaccard", {
  mk <- function(q, e) data.frame(contrast_id = "c", feature_id = sprintf("F%d", 1:6), eligibility = "tested", q_value = q, p_value = q, effect = e, stringsAsFactors = FALSE)
  a <- mk(c(0.01, 0.02, 0.2, 0.01, 0.5, 0.04), c(1, 2, 0.1, -1, 0.2, 0.5)); b <- mk(c(0.01, 0.3, 0.01, 0.02, 0.5, 0.2), c(0.8, 1.5, 0.3, -0.9, 0.1, 0.2))
  cmp <- fn("pd_compare")(a, b, "c", list(criterion = "family_q", threshold = 0.05))
  A <- a$feature_id[a$q_value < 0.05]; Bm <- b$feature_id[b$q_value < 0.05]
  testthat::expect_setequal(cmp$retained, intersect(A, Bm)); testthat::expect_setequal(cmp$lost, setdiff(A, Bm)); testthat::expect_setequal(cmp$gained, setdiff(Bm, A))
  testthat::expect_equal(cmp$jaccard, length(intersect(A, Bm)) / length(union(A, Bm)))
  testthat::expect_equal(cmp$effect_correlation, stats::cor(a$effect, b$effect), tolerance = 1e-12)
  testthat::expect_equal(cmp$attenuation_slope, unname(stats::coef(stats::lm(b$effect ~ 0 + a$effect))[1]), tolerance = 1e-12)
  testthat::expect_match(get("PD_INTERACTION_WARNING", envir = ns), "not evidence of an interaction")
})

testthat::test_that("V141 matched-n draws equal a direct seeded loop and contain whole biological units only", {
  units <- c(rep(paste0("A", 1:9), each = 2), paste0("B", 1:4))      # A units carry two observations each
  groups <- c(rep("A", 18), rep("B", 4))
  draws <- fn("pd_matched_draws")(units, groups, list(A = 4), 5L, 99L)
  set.seed(99, kind = "L'Ecuyer-CMRG")
  oracle <- lapply(1:5, function(d) c(sample(paste0("A", 1:9), 4), paste0("B", 1:4)))
  testthat::expect_identical(draws, oracle)
  for (d in draws) testthat::expect_true(all(d %in% unique(units)))   # units, never individual observations
})

testthat::test_that("mirrored sensitivity families keep the primary adjustment and only the model's contrasts", {
  fams <- list(list(family_id = "p", hypothesis_type = "protein_zero_null", adjustment = "BH", q_cutoff = 0.05, dependence_assumption = "BH",
                    members = list(list(model_id = "m", contrast_id = "c1"), list(model_id = "m", contrast_id = "c2"))))
  out <- fn("pd_mirror_families")(fams, "s", "c2", "s")
  testthat::expect_identical(out[[1]]$family_id, "p__s"); testthat::expect_identical(out[[1]]$role, "sensitivity")
  testthat::expect_identical(length(out[[1]]$members), 1L); testthat::expect_identical(out[[1]]$members[[1]]$contrast_id, "c2")
})
