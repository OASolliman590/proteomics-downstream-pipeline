# Offline GMT libraries at the R level (packet R15b): V175 (licence rule), V176 (ORA with the measured background),
# V177 (pre-ranked enrichment on Welch statistics, hand-computed and enumerated) and V178 (family discipline).
# Oracles: stats::phyper, a base-R loop over every subset, hand-computed running sums, and p.adjust-by-family.

ns <- asNamespace("proteomicsCore")
fn <- function(name) get(name, envir = ns)

testthat::test_that("V176 ORA on the measured background equals stats::phyper(k - 1, K, N - K, n, lower.tail = FALSE)", {
  universe <- sprintf("G%04d", 1:500)
  set_a <- universe[1:30]
  set_b <- universe[100:129]                       # no foreground overlap: k = 0 is retained through BH
  foreground <- c(set_a[1:8], universe[200:231])   # n = 40, k = 8 in set A
  sets <- list(A = set_a, B = set_b)
  rows <- fn("ora_test")(universe, foreground, sets, 1L, 500L)
  testthat::expect_equal(sort(rows$set_id), c("A", "B"))
  a <- rows[rows$set_id == "A", ]; b <- rows[rows$set_id == "B", ]
  oracle_a <- stats::phyper(8 - 1, 30, 500 - 30, 40, lower.tail = FALSE)
  testthat::expect_equal(a$N, 500L); testthat::expect_equal(a$K, 30L); testthat::expect_equal(a$n, 40L); testthat::expect_equal(a$k, 8L)
  testthat::expect_lt(abs(a$p_value - oracle_a) / oracle_a, 1e-12)
  testthat::expect_equal(b$k, 0L); testthat::expect_equal(b$p_value, 1)
  bh <- fn("adjust_pvalues")(rows$p_value, "BH")
  testthat::expect_equal(bh, stats::p.adjust(rows$p_value, method = "BH"), tolerance = 1e-12)
})

testthat::test_that("V176 negative: a whole-genome background (N = 20,000) is refused (E_ORA_UNIVERSE)", {
  universe <- sprintf("G%04d", 1:500)
  testthat::expect_error(fn("ora_background_check")(list(kind = "genome", size = 20000), universe), "E_ORA_UNIVERSE")
  testthat::expect_error(fn("ora_background_check")(list(kind = "measured", size = 20000), universe), "E_ORA_UNIVERSE")
  testthat::expect_invisible(fn("ora_background_check")(list(kind = "measured", size = 500), universe))
})

testthat::test_that("V177 running-sum enrichment score matches a hand computation on a six-gene ranking", {
  # Ranks t = 5, 4, 3, 2, 1, 0.5; set = positions 1 and 3 (|t| = 5 and 3).  Hit steps 5/8 and 3/8; miss step 1/4.
  # Running sum: 0.625, 0.375, 0.75, 0.5, 0.25, 0.0 -> maximum deviation 0.75 (hand computed).
  ranked <- stats::setNames(c(5, 4, 3, 2, 1, 0.5), c("g1", "g2", "g3", "g4", "g5", "g6"))
  es <- fn("prerank_running_es")(ranked, c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE))
  testthat::expect_equal(es, 0.75, tolerance = 1e-12)
})

testthat::test_that("V177 exact permutation P equals a base-R loop over every subset of the set size", {
  set.seed(177)
  t_values <- stats::setNames(round(stats::rnorm(10, 0, 2), 6), sprintf("G%02d", 1:10))
  ranked <- fn("prerank_welch_ranks")(t_values, tie_break = "sm15_stable_id")
  members <- c(2L, 5L, 9L)                         # set of size 3
  res <- fn("prerank_exact")(ranked, members)
  # oracle: enumerate every 3-subset of the 10 ranked positions with explicit loops
  subsets <- utils::combn(10, 3); scores <- numeric(ncol(subsets))
  for (j in seq_len(ncol(subsets))) {
    mask <- logical(10); mask[subsets[, j]] <- TRUE
    hit <- ifelse(mask, abs(ranked), 0); hit <- hit / sum(hit[mask])
    miss <- ifelse(mask, 0, 1 / 7); run <- cumsum(hit - miss)
    scores[j] <- run[which.max(abs(run))]
  }
  mask_obs <- logical(10); mask_obs[members] <- TRUE
  hit <- ifelse(mask_obs, abs(ranked), 0); hit <- hit / sum(hit[mask_obs]); miss <- ifelse(mask_obs, 0, 1 / 7)
  run <- cumsum(hit - miss); es_obs <- run[which.max(abs(run))]
  p_oracle <- if (es_obs >= 0) mean(scores >= es_obs - 1e-12) else mean(scores <= es_obs + 1e-12)
  testthat::expect_equal(res$es, es_obs, tolerance = 1e-12)
  testthat::expect_equal(res$p_value, p_oracle, tolerance = 1e-12)
  testthat::expect_equal(res$subsets, 120L)
})

testthat::test_that("V177 negative: duplicated gene IDs (E_RANK_DUPLICATE) and unbroken ties (E_RANK_TIE) are refused", {
  dup <- stats::setNames(c(1, 2, 3), c("g1", "g1", "g2"))
  testthat::expect_error(fn("prerank_welch_ranks")(dup, "sm15_stable_id"), "E_RANK_DUPLICATE")
  tied <- stats::setNames(c(1, 1, 2), c("g1", "g2", "g3"))
  testthat::expect_error(fn("prerank_welch_ranks")(tied, "none"), "E_RANK_TIE")
  broken <- fn("prerank_welch_ranks")(tied, "sm15_stable_id")
  testthat::expect_identical(names(broken), c("g3", "g1", "g2"))   # ties resolved by the smallest stable ID
})

testthat::test_that("V178 each library x method x contrast is its own family and BH is applied inside it", {
  set.seed(178)
  rows <- expand.grid(library_id = c("libA", "libB"), method = c("ora_up", "prerank_welch"), contrast_id = c("c1", "c2"),
                      set_id = c("s1", "s2", "s3"), stringsAsFactors = FALSE)
  rows$p_value <- round(stats::runif(nrow(rows)), 4)
  adjusted <- fn("enrichment_family_adjust")(rows, fn("gene_set_family_key"))
  testthat::expect_equal(length(unique(adjusted$family_id)), 8L)                          # 2 libraries x 2 methods x 2 contrasts
  testthat::expect_true(all(table(adjusted$family_id) == 3L))
  for (key in unique(adjusted$family_id)) {
    sel <- adjusted$family_id == key
    testthat::expect_equal(adjusted$q_value[sel], stats::p.adjust(adjusted$p_value[sel], method = "BH"), tolerance = 1e-12)
  }
})

testthat::test_that("V178 negative: a BH adjustment spanning libraries is refused (E_ENRICHMENT_FAMILY_POOLED)", {
  rows <- data.frame(library_id = c("libA", "libB"), method = "ora_up", contrast_id = "c1", set_id = c("s1", "s1"), p_value = c(0.01, 0.2),
                     stringsAsFactors = FALSE)
  pooled_key <- function(r) paste(r$method, r$contrast_id, sep = "|")
  testthat::expect_error(fn("enrichment_family_adjust")(rows, pooled_key), "E_ENRICHMENT_FAMILY_POOLED")
})

testthat::test_that("V175 the licence rule refuses restricted and non-synthetic redistribution and allows the synthetic open fixture", {
  testthat::expect_null(fn("library_redistribution_refusal")("SYNTHETIC-CC0", TRUE, TRUE))
  testthat::expect_null(fn("library_redistribution_refusal")("KEGG-restricted", FALSE, FALSE))
  testthat::expect_match(fn("library_redistribution_refusal")("KEGG-restricted", TRUE, TRUE), "restricted")
  testthat::expect_match(fn("library_redistribution_refusal")("MSigDB-restricted", TRUE, FALSE), "restricted")
  testthat::expect_match(fn("library_redistribution_refusal")("CC-BY-4.0", TRUE, FALSE), "non-synthetic")
})

testthat::test_that("V175 an artifact that contains a restricted library is refused before anything is written", {
  tmp <- tempfile(fileext = ".gmt"); on.exit(unlink(tmp), add = TRUE)
  writeLines("KEGG_X\tSYNTHETIC_FIXTURE: restricted stand-in\tG0001\tG0002", tmp)
  sha <- proteomicsCore:::sha256_file(tmp)
  libraries <- list(list(id = "kegg", sha256 = sha, licence_id = "KEGG-restricted", redistributable = FALSE, synthetic = TRUE))
  testthat::expect_error(fn("check_artifact_redistribution")(sha, libraries), "E_LIBRARY_REDISTRIBUTION")
  testthat::expect_invisible(fn("check_artifact_redistribution")("0000", list(libraries[[1]])))
})

testthat::test_that("gmt_parse reads sets and refuses duplicated names", {
  tmp <- tempfile(fileext = ".gmt"); on.exit(unlink(tmp), add = TRUE)
  writeLines(c("setA\tSYNTHETIC_FIXTURE: a\tG1\tG2\tG2", "setB\tSYNTHETIC_FIXTURE: b\tG3"), tmp)
  sets <- fn("gmt_parse")(tmp)
  testthat::expect_identical(names(sets), c("setA", "setB"))
  testthat::expect_identical(sets$setA, c("G1", "G2"))
  writeLines(c("setA\tx\tG1", "setA\ty\tG2"), tmp)
  testthat::expect_error(fn("gmt_parse")(tmp), "E_LIBRARY_FORMAT")
})
