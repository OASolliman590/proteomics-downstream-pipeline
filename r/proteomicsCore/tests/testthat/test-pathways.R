# R08 acceptance at the R level (V075 negative, V076, V077, V079).
ns <- asNamespace("proteomicsCore"); fn <- function(name) get(name, envir = ns)

testthat::test_that("V076 exact hypergeometric ORA keeps the k = 0 set through BH", {
  universe <- paste0("g", 1:10); sets <- list(hit = c("g1", "g2"), miss = c("g3", "g4"))
  o <- fn("ora_test")(universe, c("g1", "g2"), sets, 1, 10)
  testthat::expect_equal(o$p_value, c(1 / 45, 1)); testthat::expect_identical(o$k, c(2L, 0L))
  q <- fn("adjust_pvalues")(o$p_value, "BH"); testthat::expect_equal(q, c(2 / 45, 1))
  testthat::expect_false(isTRUE(all.equal(fn("adjust_pvalues")(o$p_value[1], "BH"), 2 / 45)))     # dropping the zero-hit set changes q1
  e <- fn("ora_test")(universe, character(), sets, 1, 10)
  testthat::expect_equal(e$p_value, c(1, 1)); testthat::expect_identical(unique(e$status), "empty_foreground")
})

testthat::test_that("V075 invalid ranks are refused", {
  testthat::expect_error(fn("validate_ranks")(c(a = 1, b = NA)), "E_FGSEA_RANK_INVALID")
  testthat::expect_error(fn("validate_ranks")(c(a = 1, a = 2)), "E_FGSEA_RANK_INVALID")
  testthat::expect_identical(names(fn("validate_ranks")(c(b = 1, a = 1, c = 2))), c("c", "a", "b"))   # ties broken by stable gene ID
})

testthat::test_that("V077 central family adjustment never pools different nulls", {
  rows <- data.frame(hypothesis_type = c("competitive_enrichment", "competitive_enrichment", "self_contained_mixed", "self_contained_mixed"),
                     contrast_id = "c", resource_id = "sets", p_value = c(0.01, 0.04, 0.02, 0.5), stringsAsFactors = FALSE)
  fams <- list(list(id = "cam", hypothesis_type = "competitive_enrichment", contrast_ids = list("c"), collection_ids = list("sets"), adjustment = "BH"),
               list(id = "mix", hypothesis_type = "self_contained_mixed", contrast_ids = list("c"), collection_ids = list("sets"), adjustment = "BY"))
  out <- fn("pathway_family_adjust")(rows, fams)
  testthat::expect_equal(out$q_value[1:2], stats::p.adjust(c(0.01, 0.04), "BH")); testthat::expect_equal(out$q_value[3:4], stats::p.adjust(c(0.02, 0.5), "BY"))
  testthat::expect_identical(out$family_id, c("cam", "cam", "mix", "mix"))
})

testthat::test_that("V079 leading-edge overlap is displayed without changing tests", {
  fg <- data.frame(contrast_id = "c", set_id = c("A", "B", "C"), p_value = c(0.01, 0.02, 0.5),
                   leading_edge = c('["g1","g2","g3"]', '["g2","g3","g4"]', '["g9"]'), stringsAsFactors = FALSE)
  ov <- fn("leading_edge_overlap")(fg)
  ab <- ov[ov$set_a == "c|A" & ov$set_b == "c|B", ]
  testthat::expect_identical(ab$intersection, 2L); testthat::expect_equal(ab$jaccard, 2 / 4)
  testthat::expect_identical(fg$p_value, c(0.01, 0.02, 0.5)); testthat::expect_identical(nrow(ov), 3L)
})

testthat::test_that("V080 all-nonfinite or empty pathway results fail instead of reporting zero discoveries", {
  testthat::expect_error(fn("pathway_require_finite")(data.frame(p_value = c(NA_real_, NaN))), "E_PATHWAY_NO_FINITE_TESTS")
  testthat::expect_error(fn("pathway_require_finite")(NULL), "E_PATHWAY_NO_ELIGIBLE_SETS")
  testthat::expect_silent(fn("pathway_require_finite")(data.frame(p_value = c(NA, 0.3))))
})
