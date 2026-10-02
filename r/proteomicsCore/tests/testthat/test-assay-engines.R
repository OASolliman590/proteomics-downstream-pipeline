# R06 acceptance at the R level (V051, V055).
ns <- asNamespace("proteomicsCore"); fn <- function(name) get(name, envir = ns)
counts <- function(values, source = "search_engine_export") data.frame(feature_id = paste0("F", 1:8), count_type = "peptide", count_value = as.character(values),
                                                                         count_source = source, count_aggregation = "protein_level", pseudocount_policy = "none", stringsAsFactors = FALSE)

testthat::test_that("V051 genuine counts accepted; zeros preserved; proxies and hidden +1 refused", {
  ok <- fn("validate_count_evidence")(counts(1:8), paste0("F", 1:8), "protein_level", list(mode = "reject"))
  testthat::expect_equal(unname(ok$fit), 1:8)
  testthat::expect_error(fn("validate_count_evidence")(counts(c(0, 2:8)), paste0("F", 1:8), "protein_level", list(mode = "reject")), "E_DEQMS_COUNT_NONPOSITIVE")
  off <- fn("validate_count_evidence")(counts(c(0, 2:8)), paste0("F", 1:8), "protein_level", list(mode = "explicit_offset", offset = 1, justification = "declared single-hit policy"))
  testthat::expect_equal(unname(off$original[1]), 0); testthat::expect_equal(unname(off$fit[1]), 1)
  testthat::expect_error(fn("validate_count_evidence")(counts(1:8, "observed_sample_count"), paste0("F", 1:8), "protein_level", list(mode = "reject")), "E_DEQMS_COUNT_EVIDENCE")
  testthat::expect_error(fn("validate_count_evidence")(counts(1:8), paste0("F", 1:9), "protein_level", list(mode = "reject")), "E_DEQMS_COUNT_EVIDENCE")
  testthat::expect_error(fn("validate_count_evidence")(counts(1:8), paste0("F", 1:8), "minimum_across_observations", list(mode = "reject")), "E_DEQMS_COUNT_EVIDENCE")
})

testthat::test_that("V055 proDA / DEqMS eligibility is independent of installation", {
  m <- list(engine = "proda", hypotheses = list("zero_null"))
  testthat::expect_identical(fn("assay_scientific_eligibility")(m, "none", FALSE, "tmt", "none_documented"), "E_PRODA_LFQ_UNIMPUTED_REQUIRED")
  testthat::expect_identical(fn("assay_scientific_eligibility")(m, "none", FALSE, "lfq_dda", "unknown"), "E_PRODA_LFQ_UNIMPUTED_REQUIRED")
  testthat::expect_identical(fn("assay_scientific_eligibility")(m, "duplicate_correlation", FALSE, "lfq_dda", "none_documented"), "E_PRODA_DESIGN_UNSUPPORTED")
  testthat::expect_null(fn("assay_scientific_eligibility")(m, "none", FALSE, "lfq_dia", "none_documented"))
  d <- list(engine = "deqms", hypotheses = list("zero_null", "treat"))
  testthat::expect_identical(fn("assay_scientific_eligibility")(d, "none", FALSE, "lfq_dda", "none_documented"), "E_DEQMS_HYPOTHESIS_UNSUPPORTED")
  testthat::expect_identical(fn("assay_scientific_eligibility")(d, "none", TRUE, "lfq_dda", "none_documented"), "E_DEQMS_DESIGN_UNSUPPORTED")
})
