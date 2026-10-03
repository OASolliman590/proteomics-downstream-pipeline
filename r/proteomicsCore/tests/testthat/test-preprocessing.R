# R03 preprocessing acceptance tests (V022, V023, V024, V028, V029).
# Oracles are hand arithmetic or direct package calls written in the test.

ns <- asNamespace("proteomicsCore")
fn <- function(name) get(name, envir = ns)

testthat::test_that("V022 median centering matches the analytic factors and preserves NA", {
  x <- matrix(c(10, 11, 12, NA, 12, 13, 14, 15), nrow = 4, dimnames = list(paste0("F", 1:4), c("A", "B")))
  # A median of observed = 11; B median = 13.5; median of medians = 12.25
  out <- fn("normalize_median")(x)
  testthat::expect_equal(out$factors$factor, c(11 - 12.25, 13.5 - 12.25), tolerance = 1e-12)
  testthat::expect_true(is.na(out$values["F4", "A"]))
  testthat::expect_equal(unname(out$values[, "B"]), c(12, 13, 14, 15) - 1.25, tolerance = 1e-12)
  offset <- matrix(c(5, 6, 7, 8, 7, 8, 9, 10), 4, dimnames = list(paste0("F", 1:4), c("A", "B")))   # +2 log2 offset
  centred <- fn("normalize_median")(offset)$values
  testthat::expect_equal(unname(centred[, "A"]), unname(centred[, "B"]), tolerance = 1e-12)
})

testthat::test_that("V022 reference centering uses declared proteins and refuses missing coverage", {
  x <- matrix(c(5, 6, 20, 7, 8, 30), 3, dimnames = list(c("R1", "R2", "X"), c("A", "B")))
  out <- fn("normalize_reference")(x, NULL, list("R1", "R2"), 1)
  # reference medians 5.5 and 7.5; centre 6.5 -> factors -1, +1 (X does not influence)
  testthat::expect_equal(out$factors$factor, c(-1, 1), tolerance = 1e-12)
  x2 <- x; x2["R1", "B"] <- NA; x2["R2", "B"] <- NA
  testthat::expect_error(fn("normalize_reference")(x2, NULL, list("R1", "R2"), 1), "E_REFERENCE_COVERAGE")
  testthat::expect_error(fn("normalize_reference")(x, NULL, list(), 1), "E_REFERENCE_COVERAGE")
})

testthat::test_that("V022 quantile normalization is a complete-case sensitivity via limma", {
  testthat::skip_if_not_installed("limma")
  x <- matrix(c(1, 2, 3, NA, 2, 4, 6, 8), 4, dimnames = list(paste0("F", 1:4), c("A", "B")))
  q <- fn("normalize_quantile_sensitivity")(x)
  testthat::expect_identical(q$lost_features, "F4")
  testthat::expect_equal(q$values, limma::normalizeBetweenArrays(x[1:3, ], method = "quantile"))
})

testthat::test_that("V023 TMT loading and bridge recover hand-calculated values", {
  # two plexes, three channels each: two sample channels and one bridge; channel loadings +0, +1, -1 and plex offset +3
  base <- c(F1 = 10, F2 = 12, F3 = 14)
  load <- c(0, 1, -1, 0, 1, -1)
  plex_shift <- c(0, 0, 0, 3, 3, 3)
  x <- sapply(seq_along(load), function(j) base + load[j] + plex_shift[j])
  colnames(x) <- c("p1_s1", "p1_s2", "p1_bridge", "p2_s1", "p2_s2", "p2_bridge")
  obs <- data.frame(observation_id = colnames(x), plex_id = rep(c("P1", "P2"), each = 3), channel_role = rep(c("sample", "sample", "bridge"), 2), group = c("a", "b", NA, "a", "b", NA), stringsAsFactors = FALSE)
  loaded <- fn("tmt_loading")(x, obs)
  # channel medians: base median 12 + load + shift; plex median of medians 12 + shift -> factor = load
  testthat::expect_equal(loaded$factors$factor, load, tolerance = 1e-12)
  bridged <- fn("tmt_bridge")(loaded$values, obs, obs$channel_role == "bridge")
  # after loading each plex equals base + shift; b_fp = base + shift; b_f = base + 1.5 -> all channels become base + 1.5
  testthat::expect_equal(unname(bridged$values[, "p2_s1"]), unname(base) + 1.5, tolerance = 1e-12)
  testthat::expect_equal(unname(bridged$values[, "p1_s2"]), unname(base) + 1.5, tolerance = 1e-12)
  testthat::expect_true(all(bridged$eligibility$bridge_eligible))
  partial <- loaded$values; partial["F2", "p2_bridge"] <- NA
  testthat::expect_identical(fn("tmt_bridge")(partial, obs, obs$channel_role == "bridge")$eligibility$reason[2], "tmt_bridge_measurement_missing_in_plex")
  testthat::expect_error(fn("tmt_bridge")(loaded$values, obs, c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)), "E_TMT_BRIDGE_REQUIRED")
  testthat::expect_silent(fn("tmt_no_bridge_check")(obs, "group"))
  confounded <- obs; confounded$group <- c("a", "a", NA, "b", "b", NA)
  testthat::expect_error(fn("tmt_no_bridge_check")(confounded, "group"), "E_TMT_PLEX_CONFOUNDED")
})

testthat::test_that("V024 coverage counts genuine observation, not numeric availability", {
  obs <- data.frame(observation_id = paste0("O", 1:8), group = rep(c("ctl", "trt"), each = 4), stringsAsFactors = FALSE)
  observed <- matrix(TRUE, 5, 8, dimnames = list(paste0("F", 1:5), obs$observation_id))
  observed["F2", c("O5", "O6")] <- FALSE                  # 2/4 in trt
  observed["F3", c("O5", "O6", "O7")] <- FALSE            # 1/4
  observed["F4", c("O5", "O6", "O7", "O8")] <- FALSE      # 0/4
  observed["F5", ] <- FALSE                               # all-study missing (values may be prior-imputed)
  rule <- list(model_id = "m", contrast_id = "c", required_groups = list("trt", "ctl"), policy = "available_case", minimum_observed_per_group = 2, minimum_fraction = 0.5)
  cov <- fn("coverage_tables")(observed, obs, "group", list(rule))$summary
  testthat::expect_identical(cov$eligibility, c("eligible", "eligible", "excluded", "nonestimable", "excluded"))
  testthat::expect_identical(cov$reason[3:5], c("coverage_below_minimum:trt", "all_missing_required_group:trt", "all_study_missing"))
  testthat::expect_identical(cov$n_obs_by_required_group[3], "{\"trt\":1,\"ctl\":4}")
  deferred <- rule; deferred$policy <- "native_dropout"
  testthat::expect_identical(unique(fn("coverage_tables")(observed, obs, "group", list(deferred))$summary$eligibility), "not_run")
  testthat::expect_error(fn("coverage_tables")(NULL, obs, "group", list(rule)), "E_ORIGINAL_MASK_REQUIRED")
})

testthat::test_that("V025 missingness separates observed, prior-imputed and numeric availability", {
  values <- matrix(c(1, 2, NA, 4, 5, 6), 2, dimnames = list(c("F1", "F2"), c("A", "B", "C")))
  observed <- matrix(c(TRUE, FALSE, FALSE, TRUE, TRUE, TRUE), 2, dimnames = dimnames(values))
  prior <- !is.na(values) & !observed
  obs <- data.frame(observation_id = c("A", "B", "C"), group = c("g", "g", "h"), stringsAsFactors = FALSE)
  m <- fn("missingness_tables")(values, observed, prior, obs, "group")
  testthat::expect_equal(m$by_observation$n_numeric, c(2, 1, 2))
  testthat::expect_equal(m$by_observation$n_original_observed, c(1, 1, 2))
  testthat::expect_equal(m$by_observation$n_prior_imputed, c(1, 0, 0))
  testthat::expect_identical(m$statement$mechanism_assertion, "none")
  unknown <- fn("missingness_tables")(values, NULL, NULL, obs, "group")
  testthat::expect_true(all(is.na(unknown$by_observation$n_original_observed)))
  testthat::expect_identical(unknown$statement$mask_state, "unknown")
})

testthat::test_that("V028 sensitivities follow SM05 and never touch the primary values", {
  x <- matrix(c(10, NA, 12, 13, 9, 11, NA, 14, 8, 10, 12, 15), 4, dimnames = list(paste0("F", 1:4), c("A", "B", "C")))
  original <- x
  md <- fn("run_sensitivity")(x, list(id = "md", method = "min_deterministic", model_id = "m"), 1)
  testthat::expect_equal(md$values["F2", "A"], min(x["F2", ], na.rm = TRUE))
  testthat::expect_equal(md$values["F3", "B"], min(x["F3", ], na.rm = TRUE))
  testthat::expect_error(fn("run_sensitivity")(x, list(id = "md", method = "min_deterministic", model_id = "m", shift_sd = 1.8, scale_sd = 0.3), 1), "E_SENSITIVITY_PARAMETERS")
  g1 <- fn("run_sensitivity")(x, list(id = "g", method = "left_shifted_gaussian", model_id = "m", shift_sd = 1.8, scale_sd = 0.3), 42)
  g2 <- fn("run_sensitivity")(x, list(id = "g", method = "left_shifted_gaussian", model_id = "m", shift_sd = 1.8, scale_sd = 0.3), 42)
  testthat::expect_identical(g1$values, g2$values)
  testthat::expect_true(all(g1$imputed == is.na(x)))
  flat <- matrix(c(5, 5, NA, 1, 2, 3), 3, dimnames = list(paste0("F", 1:3), c("A", "B")))
  testthat::expect_error(fn("run_sensitivity")(flat, list(id = "g", method = "left_shifted_gaussian", model_id = "m", shift_sd = 1, scale_sd = 1), 1), "E_IMPUTATION_DISTRIBUTION")
  testthat::skip_if_not_installed("impute")
  big <- matrix(c(10, 11, 12, 13, 14, 15, 10.5, NA, 12.5, 13.5, 14.5, 15.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5), 6, dimnames = list(paste0("F", 1:6), c("A", "B", "C")))
  spec <- list(id = "k", method = "knn", model_id = "m", k = 2, rowmax = 0.5, colmax = 0.8, maxp = 1500)
  kn <- fn("run_sensitivity")(big, spec, 7)
  direct <- suppressWarnings(utils::capture.output(ref <- impute::impute.knn(big, k = 2, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed = 7)))
  testthat::expect_equal(kn$values, ref$data)
  cc <- fn("run_sensitivity")(x, list(id = "cc", method = "complete_case", model_id = "m"), 1)
  testthat::expect_identical(rownames(cc$values), c("F1", "F4")); testthat::expect_identical(cc$excluded, c("F2", "F3"))
  testthat::expect_identical(x, original)
})

testthat::test_that("V029 detection Fisher P and detection-only BH match enumeration", {
  obs <- data.frame(observation_id = paste0("O", 1:8), group = rep(c("a", "b"), each = 4), stringsAsFactors = FALSE)
  observed <- rbind(F1 = c(TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE),       # 3/1 vs 1/3
                    F2 = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE))      # 4/0 vs 0/4
  colnames(observed) <- obs$observation_id
  # fixed-margin enumeration: P(x) = choose(4,x) choose(4,4-x) / choose(8,4) = (1,16,36,16,1)/70
  p_f1 <- (1 + 16 + 16 + 1) / 70; p_f2 <- (1 + 1) / 70
  bh <- c(min(1, p_f2 * 2 / 1), p_f1 * 2 / 2)
  res <- fn("detection_tests")(observed, obs, "group", list(list(contrast_id = "a-b", required_groups = list("a", "b"))), "none", "detection-family")
  testthat::expect_equal(res$p_value, c(p_f1, p_f2), tolerance = 1e-12)
  testthat::expect_equal(res$q_value[order(res$p_value)], sort(pmin(bh, 1)), tolerance = 1e-12)
  testthat::expect_true(all(res$family_id == "detection-family" & res$hypothesis_type == "detection"))
  testthat::expect_error(fn("detection_tests")(observed, obs, "group", list(list(contrast_id = "a-b", required_groups = list("a", "b"))), "fixed_subject"), "E_DETECTION_DESIGN_UNSUPPORTED")
})

# audit 2026-10-02: the vectorized coverage_tables must equal the former per-feature loop (kept here verbatim as the oracle).
old_coverage_tables <- function(observed, observations, group_column, rules, extra_exclusions = NULL) {
  if (is.null(observed)) stop("E_ORIGINAL_MASK_REQUIRED: primary observed-coverage rules need a genuine original-observed mask", call. = FALSE)
  groups <- observations[[group_column]]
  all_missing <- rowSums(observed) == 0
  summary_rows <- list(); group_rows <- list()
  for (rule in rules) {
    required <- unlist(rule$required_groups)
    min_n <- if (is.null(rule$minimum_observed_per_group)) 2L else as.integer(rule$minimum_observed_per_group)
    min_fraction <- if (is.null(rule$minimum_fraction)) 0.5 else as.numeric(rule$minimum_fraction)
    policy <- if (is.null(rule$policy)) "available_case" else rule$policy
    n_units <- vapply(required, function(g) sum(groups == g), integer(1))
    counts <- vapply(required, function(g) rowSums(observed[, groups == g, drop = FALSE]), numeric(nrow(observed)))
    counts <- matrix(counts, nrow = nrow(observed), dimnames = list(rownames(observed), required))
    for (i in seq_len(nrow(observed))) {
      feature <- rownames(observed)[i]
      n_obs <- counts[i, ]
      fraction <- ifelse(n_units > 0, n_obs / n_units, NA_real_)
      if (policy == "native_dropout") { eligibility <- "not_run"; reason <- "native_dropout_policy_deferred_to_R06" }
      else if (!is.null(extra_exclusions) && feature %in% names(extra_exclusions)) { eligibility <- "excluded"; reason <- extra_exclusions[[feature]] }
      else if (all_missing[i]) { eligibility <- "excluded"; reason <- "all_study_missing" }
      else if (any(n_units == 0)) { eligibility <- "nonestimable"; reason <- paste0("required_group_has_no_units:", paste(required[n_units == 0], collapse = ",")) }
      else if (any(n_obs == 0)) { eligibility <- "nonestimable"; reason <- paste0("all_missing_required_group:", paste(required[n_obs == 0], collapse = ",")) }
      else if (any(n_obs < min_n | fraction < min_fraction)) { eligibility <- "excluded"; reason <- paste0("coverage_below_minimum:", paste(required[n_obs < min_n | fraction < min_fraction], collapse = ",")) }
      else { eligibility <- "eligible"; reason <- "eligible" }
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(model_id = rule$model_id, contrast_id = rule$contrast_id, feature_id = feature,
        eligibility = eligibility, reason = reason, n_obs_by_required_group = as.character(jsonlite::toJSON(as.list(stats::setNames(as.integer(n_obs), required)), auto_unbox = TRUE)),
        policy = policy, minimum_observed_per_group = min_n, minimum_fraction = min_fraction, stringsAsFactors = FALSE)
    }
    group_rows[[length(group_rows) + 1L]] <- data.frame(model_id = rule$model_id, contrast_id = rule$contrast_id,
      feature_id = rep(rownames(observed), times = length(required)), group = rep(required, each = nrow(observed)),
      n_observed = as.integer(as.vector(counts)), n_units = rep(as.integer(n_units), each = nrow(observed)),
      observed_fraction = as.vector(counts) / rep(n_units, each = nrow(observed)), stringsAsFactors = FALSE)
  }
  list(summary = if (length(summary_rows)) do.call(rbind, summary_rows) else data.frame(),
       by_group = if (length(group_rows)) do.call(rbind, group_rows) else data.frame())
}
testthat::test_that("audit: vectorized coverage_tables equals the former per-feature loop", {
  set.seed(5); g <- rep(c("A", "B", "C"), c(4, 4, 3))
  obs <- matrix(runif(80 * 11) > 0.3, 80, dimnames = list(sprintf("P%02d", 1:80), paste0("O", 1:11))); obs[3, ] <- FALSE; obs[4, g == "A"] <- FALSE
  observations <- data.frame(observation_id = colnames(obs), group = g, stringsAsFactors = FALSE)
  rules <- list(list(model_id = "m", contrast_id = "B-A", required_groups = list("B", "A"), minimum_observed_per_group = 2, minimum_fraction = 0.5),
                list(model_id = "m", contrast_id = "D-A", required_groups = list("D", "A")),
                list(model_id = "m", contrast_id = "C-B", required_groups = list("C", "B"), minimum_observed_per_group = 3, minimum_fraction = 0.9),
                list(model_id = "p", contrast_id = "B-A", required_groups = list("B", "A"), policy = "native_dropout"))
  extra <- list(P07 = "prespecified_exclusion")
  new <- get("coverage_tables", envir = asNamespace("proteomicsCore"))(obs, observations, "group", rules, extra)
  old <- old_coverage_tables(obs, observations, "group", rules, extra)
  rownames(new$summary) <- NULL; rownames(old$summary) <- NULL
  testthat::expect_identical(new$summary, old$summary); testthat::expect_identical(new$by_group, old$by_group)
})
