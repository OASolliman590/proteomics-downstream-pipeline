# R14a result-structure acceptance at the R level (V131, V133-V136).  Oracles are base R set operations,
# brute-force enumeration, cor(), lm() and phyper() written here, independent of the adapter.

ns <- asNamespace("proteomicsCore")
fn <- function(name) get(name, envir = ns)

endpoint_table <- function(seed = 1, n = 40) {
  set.seed(seed)
  make <- function(contrast) data.frame(model_id = "m", contrast_id = contrast, hypothesis_type = "protein_zero_null", family_id = "fam",
                                        feature_id = sprintf("F%02d", seq_len(n)), eligibility = c(rep("tested", n - 3), rep("excluded", 3)),
                                        effect = rnorm(n), p_value = runif(n), q_value = runif(n, 0, 0.2), stringsAsFactors = FALSE)
  rbind(make("c1"), make("c2"))
}
leaf <- function(contrast, criterion = "family_q", threshold = 0.05, direction = "any")
  list(family_id = "fam", model_id = "m", contrast_id = contrast, hypothesis_type = "protein_zero_null", criterion = criterion, threshold = threshold, direction = direction)

testthat::test_that("V131 set algebra equals base R set operations; raw-P leaves use p", {
  z <- endpoint_table()
  t1 <- z[z$contrast_id == "c1" & z$eligibility == "tested", ]; t2 <- z[z$contrast_id == "c2" & z$eligibility == "tested", ]
  a <- t1$feature_id[t1$q_value < 0.05]; b <- t2$feature_id[t2$q_value < 0.05]
  ev <- function(rule) sort(fn("pd_eval_rule")(rule, z, NULL)$members)
  testthat::expect_identical(ev(list(op = "union", args = list(leaf("c1"), leaf("c2")))), sort(union(a, b)))
  testthat::expect_identical(ev(list(op = "intersect", args = list(leaf("c1"), leaf("c2")))), sort(intersect(a, b)))
  testthat::expect_identical(ev(list(op = "difference", args = list(leaf("c1"), leaf("c2")))), sort(setdiff(a, b)))
  testthat::expect_identical(ev(list(op = "complement_within_tested", args = list(leaf("c1")))), sort(setdiff(t1$feature_id, a)))
  testthat::expect_identical(ev(leaf("c1", direction = "up")), sort(t1$feature_id[t1$q_value < 0.05 & t1$effect > 0]))
  testthat::expect_identical(ev(leaf("c1", "exploratory_raw_p", 0.1)), sort(t1$feature_id[t1$p_value < 0.1]))
  testthat::expect_error(fn("pd_eval_rule")(list(op = "xor", args = list(leaf("c1"), leaf("c2"))), z, NULL), "E_SETRULE_GRAMMAR")
  testthat::expect_error(fn("pd_leaf")(leaf("c9"), z, NULL), "E_SETRULE_UNKNOWN_CONTRAST")
})

testthat::test_that("V133 exclusive regions equal brute-force enumeration for k = 3 and k = 5 and sum to the union", {
  set.seed(4)
  universe <- sprintf("P%03d", 1:60)
  for (k in c(3L, 5L)) {
    sets <- stats::setNames(lapply(seq_len(k), function(i) sample(universe, 20)), paste0("S", seq_len(k)))
    genes <- stats::setNames(lapply(seq_along(universe), function(i) if (i %% 7 == 0) character() else paste0("G", (i + 1) %/% 2)), universe)
    regions <- fn("pd_regions")(sets, genes)
    testthat::expect_identical(nrow(regions), as.integer(2^k - 1))
    all_members <- unique(unlist(sets))
    for (r in seq_len(nrow(regions))) {
      inside <- strsplit(regions$defining_sets[r], ";")[[1]]
      oracle <- all_members[vapply(all_members, function(f) all(vapply(names(sets), function(s) (f %in% sets[[s]]) == (s %in% inside), logical(1))), logical(1))]
      testthat::expect_identical(regions$n_protein_groups[r], length(oracle))
      testthat::expect_setequal(Filter(nzchar, strsplit(regions$members[r], ";")[[1]]), oracle)
      testthat::expect_identical(regions$n_unique_genes[r], length(unique(unlist(genes[oracle]))))
    }
    testthat::expect_identical(sum(regions$n_protein_groups), length(all_members))
    testthat::expect_true(any(regions$n_unique_genes != regions$n_protein_groups))      # shared genes: counts differ
  }
})

testthat::test_that("V134 direction-aware overlap equals hand labels; discordant features are not concordant overlap", {
  ea <- c(F1 = 1, F2 = 2, F3 = -1, F4 = 1, F5 = 3); eb <- c(F1 = 0.5, F2 = 1, F3 = -2, F4 = -1, F6 = 1)
  ov <- fn("pd_direction_overlap")(c("F1", "F2", "F3", "F4", "F5"), c("F1", "F2", "F3", "F4", "F6"), ea, eb)
  testthat::expect_setequal(ov$concordant_up, c("F1", "F2")); testthat::expect_identical(ov$concordant_down, "F3")
  testthat::expect_identical(ov$discordant, "F4"); testthat::expect_identical(ov$only_a, "F5"); testthat::expect_identical(ov$only_b, "F6")
})

testthat::test_that("V135 concordance equals cor() and lm(y ~ 0 + x) within 1e-10 and reports no P", {
  set.seed(9); x <- stats::setNames(rnorm(50), sprintf("F%02d", 1:50)); y <- 0.6 * x + rnorm(50, 0, 0.3); names(y) <- names(x)
  y["F03"] <- NA
  cc <- fn("pd_concordance")(x, y, c("F01", "F02", "F04"))
  keep <- names(x) != "F03"
  testthat::expect_equal(cc$pearson, stats::cor(x[keep], y[keep]), tolerance = 1e-10)
  testthat::expect_equal(cc$spearman, stats::cor(x[keep], y[keep], method = "spearman"), tolerance = 1e-10)
  testthat::expect_equal(cc$origin_slope, unname(stats::coef(stats::lm(y[keep] ~ 0 + x[keep]))[1]), tolerance = 1e-10)
  testthat::expect_identical(cc$n_sign_agree, sum(sign(x[c("F01", "F02", "F04")]) == sign(y[c("F01", "F02", "F04")])))
  testthat::expect_false(any(grepl("^p_|p_value|q_value", names(cc))))
})

testthat::test_that("V136 overlap P is the phyper upper tail within the universe tested in both; shared units are detected", {
  universe <- sprintf("F%02d", 1:40); a <- c(sprintf("F%02d", 1:10), "X1"); b <- c(sprintf("F%02d", 6:15), "X2")
  o <- fn("pd_overlap_p")(a, b, universe)
  testthat::expect_identical(o$overlap, 5L); testthat::expect_identical(o$universe_n, 40L)
  testthat::expect_equal(o$p_value, stats::phyper(4, 10, 30, 10, lower.tail = FALSE), tolerance = 1e-14)
  obs <- data.frame(observation_id = c("a1", "a2", "b1", "b2", "c1", "c2"), biological_unit_id = c("a1", "a2", "b1", "b2", "c1", "c2"),
                    group = c("A", "A", "B", "B", "C", "C"), stringsAsFactors = FALSE)
  ba <- fn("pd_contrast_units")(list(required_groups = list("B", "A")), obs, "group", NULL)
  ca <- fn("pd_contrast_units")(list(required_groups = list("C", "A")), obs, "group", NULL)
  cb <- fn("pd_contrast_units")(list(required_groups = list("C")), obs, "group", NULL)
  testthat::expect_setequal(intersect(ba$units, ca$units), c("a1", "a2"))
  testthat::expect_length(intersect(ba$units, cb$units), 0L)
})

testthat::test_that("claim labels outside the frozen vocabulary are refused", {
  testthat::expect_silent(fn(".pd_check_claims")(data.frame(claim_label = c("descriptive", "module_level"))))
  testthat::expect_error(fn(".pd_check_claims")(data.frame(claim_label = "validated")), "E_BIOMARKER_CLAIM")
  testthat::expect_error(fn(".pd_check_claims")(data.frame(x = 1)), "E_BIOMARKER_CLAIM")
})
