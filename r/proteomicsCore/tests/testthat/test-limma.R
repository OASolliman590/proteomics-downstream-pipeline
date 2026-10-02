# R05 acceptance tests at the R API level (V043 negative, V046, V048 negative).

ns <- asNamespace("proteomicsCore")
fn <- function(name) get(name, envir = ns)

rows_for <- function() {
  data.frame(model_id = "m", contrast_id = c("c1", "c1", "c2", "c2", "aux"), feature_id = c("F1", "F2", "F1", "F2", "F1"),
             hypothesis_type = "protein_zero_null", eligibility = "tested", p_value = c(.01, .5, .02, .8, .001),
             q_value = NA_real_, family_id = NA_character_, stringsAsFactors = FALSE)
}
family <- function(id, contrasts, adjustment = "BH", role = "primary") list(family_id = id, hypothesis_type = "protein_zero_null", role = role, adjustment = adjustment,
  dependence_assumption = "test", q_cutoff = 0.05, members = lapply(contrasts, function(cc) list(model_id = "m", contrast_id = cc)))

testthat::test_that("V046 primary BH pools every planned contrast once; secondary stays separate", {
  rows <- rows_for()
  primary <- fn("adjust_family")(rows, family("primary", c("c1", "c2")))
  rows <- primary$rows
  secondary <- fn("adjust_family")(rows, family("secondary", "aux", role = "secondary"))
  rows <- secondary$rows
  # hand BH over p = (.01, .5, .02, .8): q = (.04, 2/3, .04, .8)
  testthat::expect_equal(rows$q_value[1:4], c(.04, 2 / 3, .04, .8), tolerance = 1e-12)
  testthat::expect_equal(rows$q_value[5], .001)
  testthat::expect_identical(primary$summary$n_planned, 4L); testthat::expect_identical(primary$summary$n_finite, 4L)
  testthat::expect_identical(primary$summary$completeness, "complete"); testthat::expect_identical(primary$summary$rejection_count, 2L)
  # BY uses the harmonic factor 1 + 1/2 + 1/3 + 1/4
  by <- fn("adjust_family")(rows_for(), family("primary", c("c1", "c2"), "BY"))$rows
  testthat::expect_equal(by$q_value[1:4], pmin(1, c(.04, 2 / 3, .04, .8) * (1 + 1 / 2 + 1 / 3 + 1 / 4)), tolerance = 1e-12)
  # negative: per-contrast BH gives different values and would fail this oracle
  per_contrast <- c(stats::p.adjust(c(.01, .5), "BH"), stats::p.adjust(c(.02, .8), "BH"))
  testthat::expect_false(isTRUE(all.equal(per_contrast, c(.04, 2 / 3, .04, .8))))
})

testthat::test_that("V048 numerical failure marks the family incomplete, never P=1", {
  rows <- rows_for(); rows$eligibility[2] <- "numerical_failure"; rows$p_value[2] <- NA
  out <- fn("adjust_family")(rows, family("primary", c("c1", "c2")))
  testthat::expect_identical(out$summary$completeness, "incomplete"); testthat::expect_identical(out$summary$n_numerical_failure, 1L)
  testthat::expect_true(is.na(out$rows$q_value[2]))
  testthat::expect_equal(out$rows$q_value[c(1, 3, 4)], stats::p.adjust(c(.01, .02, .8), "BH"), tolerance = 1e-12)
})

testthat::test_that("V048 negative: table verification rejects missing rows, NA-to-zero and orphan q", {
  base <- fn("endpoint_frame")(2)
  base$model_id <- "m"; base$contrast_id <- "c"; base$feature_id <- c("F1", "F2"); base$hypothesis_type <- "protein_zero_null"
  base$eligibility <- c("tested", "excluded"); base$reason_code <- c(NA, "coverage_below_minimum:A"); base$p_value <- c(0.2, NA); base$family_id <- "fam"; base$q_value <- c(0.2, NA)
  keys <- paste("m", "c", c("F1", "F2"), "protein_zero_null", sep = "|")
  testthat::expect_silent(fn("verify_differential_table")(base, keys))
  testthat::expect_error(fn("verify_differential_table")(base[1, ], keys), "E_TABLE_INCOMPLETE")
  zeroed <- base; zeroed$p_value[2] <- 0
  testthat::expect_error(fn("verify_differential_table")(zeroed, keys), "E_TABLE_SCHEMA")
  orphan <- base; orphan$family_id[1] <- NA
  testthat::expect_error(fn("verify_differential_table")(orphan, keys), "central q")
})

testthat::test_that("V043 negative: robust=TRUE without statmod fails preflight instead of switching robust off", {
  testthat::expect_error(fn("limma_preflight")(TRUE, installed = function(pkg) pkg != "statmod"), "E_DEPENDENCY_UNAVAILABLE")
  testthat::expect_silent(fn("limma_preflight")(FALSE, installed = function(pkg) pkg != "statmod"))
})

testthat::test_that("V045 interval uses the moderated-t quantile, not a normal quantile", {
  ci <- fn("ci_moderated_t")(1, 0.5, 6, 0.95)
  testthat::expect_equal(unname(ci[1, "upper"]), 1 + stats::qt(0.975, 6) * 0.5)
  testthat::expect_false(isTRUE(all.equal(unname(ci[1, "upper"]), 1 + stats::qnorm(0.975) * 0.5)))
})
