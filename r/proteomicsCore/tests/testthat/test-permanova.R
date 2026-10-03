# R13 PERMANOVA acceptance tests at the R level (V121-V128).  Oracles are direct
# base-R / vegan / permute calls written here, independent of the adapter.

ns <- asNamespace("proteomicsCore")
fn <- function(name) get(name, envir = ns)
skip_vegan <- function() { testthat::skip_if_not_installed("vegan"); testthat::skip_if_not_installed("permute") }

make_location <- function(seed = 11, n = 6, p = 8, groups = c("A", "B", "C"), shift = 1.2) {
  set.seed(seed)
  g <- rep(groups, each = n)
  Y <- sapply(seq_along(g), function(j) 10 + seq_len(p) * 0.3 + rnorm(p, 0, 0.4) + ifelse(g[j] == "B" & seq_len(p) <= p / 2, shift, 0))
  dimnames(Y) <- list(sprintf("F%02d", seq_len(p)), paste0(g, rep(seq_len(n), length(groups))))
  list(Y = Y, g = g, obs = data.frame(observation_id = colnames(Y), group = g, stringsAsFactors = FALSE))
}

settings_for <- function(levels, ...) {
  s <- list(group_column = "group", group_levels = levels, subject_column = NULL, metric = "euclidean", scaling = "zscore", permutations = 199L,
            seed = 7L, alpha = 0.05, pairwise = TRUE, adjustment = "holm", covariates = list(), interaction = FALSE)
  extra <- list(...)
  for (name in names(extra)) s[[name]] <- extra[[name]]   # plain replacement (modifyList would merge nested lists)
  s
}

testthat::test_that("V121 z-scored Euclidean distance equals dist(scale(t(X))) and PCA space; exclusions have reasons", {
  skip_vegan()
  d <- make_location(); Y <- d$Y; Y["F02", 3] <- NA; Y["F03", ] <- 5
  use <- fn("pm_usable_features")(Y, rownames(Y), rep(TRUE, ncol(Y)))
  testthat::expect_identical(use$reason[use$feature_id == "F02"], "incomplete_in_test_observations")
  testthat::expect_identical(use$reason[use$feature_id == "F03"], "constant_within_test")
  feats <- use$feature_id[use$used]
  X <- fn("pm_scaled")(Y, feats, rep(TRUE, ncol(Y)), "zscore")
  D <- fn("pm_distance")(X, "euclidean")
  oracle <- stats::dist(scale(t(Y[feats, ])))
  testthat::expect_equal(as.vector(D), as.vector(oracle), tolerance = 1e-10)
  pcs <- stats::prcomp(t(Y[feats, ]), center = TRUE, scale. = TRUE)$x
  testthat::expect_equal(as.vector(D), as.vector(stats::dist(pcs)), tolerance = 1e-10)
  testthat::expect_equal(as.vector(fn("pm_distance")(X, "manhattan")), as.vector(vegan::vegdist(X, "manhattan")), tolerance = 1e-12)
  testthat::expect_error(fn("pm_resolve_sets")(list(list(id = "p", kind = "declared_panel", feature_ids = list("F01", "NOPE"), selection_provenance = "same_data")), rownames(Y), NULL), "E_PERMANOVA_FEATURE_SET")
  refused <- fn("pm_feature_set")(Y, list(id = "tiny", features = c("F01", "F03")), d$obs, settings_for(c("A", "B", "C")))
  testthat::expect_true(refused$refused)
})

testthat::test_that("V122 Df, SS, R2, pseudo-F match hand sums of squares; P matches a direct seeded adonis2", {
  skip_vegan()
  d <- make_location(groups = c("A", "B"))
  X <- scale(t(d$Y)); D <- stats::dist(X); n <- nrow(X); g <- factor(d$g)
  D2 <- as.matrix(D)^2
  ss_total <- sum(D2[upper.tri(D2)]) / n
  ss_within <- sum(vapply(levels(g), function(l) { m <- D2[g == l, g == l]; sum(m[upper.tri(m)]) / sum(g == l) }, numeric(1)))
  ss_between <- ss_total - ss_within; k <- nlevels(g)
  f_oracle <- (ss_between / (k - 1)) / (ss_within / (n - k))
  res <- fn("pm_feature_set")(d$Y, list(id = "all", features = rownames(d$Y)), d$obs, settings_for(c("A", "B")))
  row <- res$tests[res$tests$analysis == "global", ]
  testthat::expect_equal(row$df, 1); testthat::expect_equal(row$ss, ss_between, tolerance = 1e-10)
  testthat::expect_equal(row$r2, ss_between / ss_total, tolerance = 1e-10); testthat::expect_equal(row$pseudo_f, f_oracle, tolerance = 1e-10)
  set.seed(7); direct <- vegan::adonis2(D ~ g, data = data.frame(g = g), permutations = permute::how(nperm = 199), by = "terms")
  testthat::expect_equal(row$p_value, direct["g", "Pr(>F)"]); testthat::expect_identical(row$nperm, 199L)
  testthat::expect_identical(fn("pm_display")(1 / 200, 199), "< 0.005"); testthat::expect_identical(fn("pm_display")(0.0421, 199), "0.0421")
  testthat::expect_true(all(res$tests$p_value >= res$tests$p_floor))
})

testthat::test_that("V123 global and pairwise tests with Holm (and BH) adjustment within a feature set", {
  skip_vegan()
  d <- make_location()
  res <- fn("pm_feature_set")(d$Y, list(id = "all", features = rownames(d$Y)), d$obs, settings_for(c("A", "B", "C")))
  pw <- res$tests[res$tests$analysis == "pairwise", ]
  testthat::expect_equal(nrow(pw), 3L); testthat::expect_equal(sum(res$tests$analysis == "global"), 1L)
  for (i in seq_len(nrow(pw))) {
    pair <- rev(strsplit(pw$comparison[i], " vs ")[[1]]); keep <- d$g %in% pair
    X <- scale(t(d$Y[, keep])); gg <- factor(d$g[keep], levels = pair)
    set.seed(7); direct <- vegan::adonis2(stats::dist(X) ~ gg, data = data.frame(gg = gg), permutations = permute::how(nperm = 199), by = "terms")
    testthat::expect_equal(pw$p_value[i], direct["gg", "Pr(>F)"]); testthat::expect_equal(pw$r2[i], direct["gg", "R2"], tolerance = 1e-10)
  }
  testthat::expect_equal(pw$p_adjusted, stats::p.adjust(pw$p_value, "holm"))
  testthat::expect_true(is.na(res$tests$p_adjusted[res$tests$analysis == "global"]))
  bh <- fn("pm_feature_set")(d$Y, list(id = "all", features = rownames(d$Y)), d$obs, settings_for(c("A", "B", "C"), adjustment = "BH"))
  testthat::expect_equal(bh$tests$p_adjusted[bh$tests$analysis == "pairwise"], stats::p.adjust(pw$p_value, "BH"))
})

testthat::test_that("V124 covariate terms use within-group permutations; the group term never comes from the blocked model", {
  skip_vegan()
  d <- make_location(groups = c("A", "B"), n = 8); d$obs$sex <- rep(c("F", "M"), 8)
  res <- fn("pm_feature_set")(d$Y, list(id = "all", features = rownames(d$Y)), d$obs, settings_for(c("A", "B"), covariates = list(list(column = "sex", type = "categorical"))))
  X <- scale(t(d$Y)); D <- stats::dist(X); data <- data.frame(group = factor(d$g), cov_sex = factor(d$obs$sex))
  set.seed(7); free <- vegan::adonis2(D ~ group + cov_sex, data = data, permutations = permute::how(nperm = 199), by = "margin")
  set.seed(7); blocked <- vegan::adonis2(D ~ group + cov_sex, data = data, permutations = permute::how(nperm = 199, blocks = data$group), by = "margin")
  adj <- res$tests[res$tests$analysis == "group_adjusted", ]; cov <- res$tests[res$tests$analysis == "covariate", ]
  testthat::expect_equal(adj$p_value, free["group", "Pr(>F)"]); testthat::expect_identical(adj$permutation_scheme, "free")
  testthat::expect_equal(cov$p_value, blocked["cov_sex", "Pr(>F)"]); testthat::expect_identical(cov$permutation_scheme, "within:group")
  testthat::expect_equal(cov$ss, blocked["cov_sex", "SumOfSqs"], tolerance = 1e-10)
  testthat::expect_false(any(res$tests$term == "group" & grepl("^within:group", res$tests$permutation_scheme)))
  # audit 2026-10-02: check the pipeline's own covariate permutation design, not a hand-built permute::how
  set.seed(1); perms <- permute::shuffleSet(16, 50, control = fn("pm_covariate_how")(d$obs$group, 199L))
  testthat::expect_true(all(apply(perms, 1L, function(idx) all(data$group[idx] == data$group))))
})

testthat::test_that("V125 interaction only when every cell is filled; empty cell is refused", {
  skip_vegan()
  d <- make_location(groups = c("A", "B"), n = 8); d$obs$sex <- rep(c("F", "M"), 8)
  res <- fn("pm_feature_set")(d$Y, list(id = "all", features = rownames(d$Y)), d$obs, settings_for(c("A", "B"), covariates = list(list(column = "sex", type = "categorical")), interaction = TRUE))
  X <- scale(t(d$Y)); data <- data.frame(group = factor(d$g), cov_sex = factor(d$obs$sex))
  set.seed(7); direct <- vegan::adonis2(stats::dist(X) ~ group * cov_sex, data = data, permutations = permute::how(nperm = 199, blocks = data$group), by = "margin")
  it <- res$tests[res$tests$analysis == "interaction", ]
  testthat::expect_equal(it$p_value, direct["group:cov_sex", "Pr(>F)"])
  d$obs$sex[d$g == "B"] <- "M"
  # audit 2026-10-02: an empty cell refuses only the interaction term (typed), never the other results
  refused <- fn("pm_feature_set")(d$Y, list(id = "all", features = rownames(d$Y)), d$obs, settings_for(c("A", "B"), covariates = list(list(column = "sex", type = "categorical")), interaction = TRUE))
  testthat::expect_identical(refused$refusals$reason_code, "E_PERMANOVA_INTERACTION_NONESTIMABLE")
  testthat::expect_false("interaction" %in% refused$tests$analysis)
})

testthat::test_that("V126 PERMDISP matches permutest and labels dispersion-only versus location shift", {
  skip_vegan()
  set.seed(3); n <- 8; p <- 10
  dev <- matrix(rnorm(n * p, 0, 0.3), p); dev <- dev - rowMeans(dev)
  Y <- cbind(10 + dev, 10 - 4 * dev); dimnames(Y) <- list(sprintf("F%02d", 1:p), c(paste0("A", 1:n), paste0("B", 1:n)))
  obs <- data.frame(observation_id = colnames(Y), group = rep(c("A", "B"), each = n), stringsAsFactors = FALSE)
  res <- fn("pm_feature_set")(Y, list(id = "all", features = rownames(Y)), obs, settings_for(c("A", "B")))
  row <- res$tests[res$tests$analysis == "global", ]
  X <- scale(t(Y)); bd <- vegan::betadisper(stats::dist(X), factor(obs$group)); set.seed(7); pt <- vegan::permutest(bd, permutations = permute::how(nperm = 199))
  testthat::expect_equal(row$permdisp_p, pt$tab[1, "Pr(>F)"]); testthat::expect_equal(row$permdisp_f, pt$tab[1, "F"], tolerance = 1e-10)
  testthat::expect_lt(row$permdisp_p, 0.05); testthat::expect_gt(row$p_value, 0.05)
  testthat::expect_identical(row$interpretation, "dispersion_difference_location_not_established")
  loc <- make_location(groups = c("A", "B"), shift = 2)
  lrow <- fn("pm_feature_set")(loc$Y, list(id = "all", features = rownames(loc$Y)), loc$obs, settings_for(c("A", "B")))$tests[1, ]
  testthat::expect_lt(lrow$p_value, 0.05)   # exact location label: see "audit m2: exact interpretation labels" below
  testthat::expect_identical(fn("pm_interpret")(0.01, 0.01, 0.05), "location_and_or_dispersion")
})

testthat::test_that("V127 per-feature R2 match lm and their mean equals the multivariate R2", {
  skip_vegan()
  d <- make_location()
  res <- fn("pm_feature_set")(d$Y, list(id = "all", features = rownames(d$Y)), d$obs, settings_for(c("A", "B", "C")))
  Z <- scale(t(d$Y)); g <- factor(d$g)
  oracle <- apply(Z, 2L, function(z) { a <- stats::anova(stats::lm(z ~ g)); a[["Sum Sq"]][1] / sum(a[["Sum Sq"]]) })
  testthat::expect_equal(res$feature_r2$univariate_group_r2, unname(oracle), tolerance = 1e-10)
  testthat::expect_identical(res$identity$state, "pass"); testthat::expect_equal(mean(oracle), res$global_r2, tolerance = 1e-8)
  manhattan <- fn("pm_feature_set")(d$Y, list(id = "all", features = rownames(d$Y)), d$obs, settings_for(c("A", "B", "C"), metric = "manhattan"))
  testthat::expect_identical(manhattan$identity$state, "not_applicable")
  original <- fn("pm_univariate_r2")
  utils::assignInNamespace("pm_univariate_r2", function(X, groups) original(X, groups) + 0.01, ns = "proteomicsCore")
  on.exit(utils::assignInNamespace("pm_univariate_r2", original, ns = "proteomicsCore"), add = TRUE)
  testthat::expect_error(fn("pm_feature_set")(d$Y, list(id = "all", features = rownames(d$Y)), d$obs, settings_for(c("A", "B", "C"))), "E_PERMANOVA_IDENTITY")
})

testthat::test_that("V128 random equal-size null, best-possible set and optimism flag", {
  r2 <- stats::setNames(seq(0.05, 0.95, length.out = 20), sprintf("F%02d", 1:20))
  a <- fn("pm_selection_context")("panel", 0.9, 3L, r2, 300L, 99L, TRUE)
  b <- fn("pm_selection_context")("panel", 0.9, 3L, r2, 300L, 99L, TRUE)
  testthat::expect_identical(a$null, b$null); testthat::expect_equal(nrow(a$null), 300L)
  testthat::expect_equal(a$summary$best_possible_r2, mean(sort(r2, decreasing = TRUE)[1:3]))
  set.seed(99); manual <- replicate(300, mean(sample(r2, 3)))
  testthat::expect_equal(a$null$r2, manual)
  testthat::expect_true(a$summary$in_sample_optimistic)
  testthat::expect_equal(a$summary$fraction_null_at_or_above_set, mean(manual >= 0.9 - 1e-12))
})

# ----------------------------------------------------------------------------- audit 2026-10-02
# Known-answer fixtures for subject-blocked designs. 'between': group constant within each subject
# (whole subjects must be permuted); 'within': every subject observed in every group (permute within subject).
make_subjects <- function(kind, shift = 1.5, n_subjects = 8, reps = 2, p = 10, seed = 5) {
  set.seed(seed)
  subjects <- sprintf("S%d", seq_len(n_subjects))
  if (kind %in% c("between", "between_null")) {
    obs <- data.frame(subject = rep(subjects, each = reps), group = rep(ifelse(seq_len(n_subjects) <= n_subjects / 2, "A", "B"), each = reps), stringsAsFactors = FALSE)
  } else if (kind %in% c("within", "within_null")) {
    obs <- data.frame(subject = rep(subjects, each = 2), group = rep(c("A", "B"), n_subjects), stringsAsFactors = FALSE)
  } else {  # mixed: half the subjects in both groups, half in one group only
    obs <- data.frame(subject = rep(subjects, each = 2), group = c(rep(c("A", "B"), n_subjects / 2), rep(c("A", "A", "B", "B"), n_subjects / 4)), stringsAsFactors = FALSE)
  }
  obs$observation_id <- paste0(obs$subject, "_", stats::ave(seq_len(nrow(obs)), obs$subject, FUN = seq_along))
  subject_effect <- matrix(rnorm(p * n_subjects, 0, 0.3), p, dimnames = list(NULL, subjects))
  eff <- if (kind %in% c("between", "within")) shift else 0
  Y <- sapply(seq_len(nrow(obs)), function(j) 10 + seq_len(p) * 0.2 + subject_effect[, obs$subject[j]] + rnorm(p, 0, 0.3) +
                ifelse(obs$group[j] == "B" & seq_len(p) <= p / 2, eff, 0))
  dimnames(Y) <- list(sprintf("F%02d", seq_len(p)), obs$observation_id)
  list(Y = Y, obs = obs)
}
blocked <- function(...) settings_for(c("A", "B"), subject_column = "subject", permutations = 999L, ...)

testthat::test_that("audit M2: group constant within subjects permutes whole subjects (known answer, not degenerate)", {
  skip_vegan()
  d <- make_subjects("between")
  row <- fn("pm_feature_set")(d$Y, list(id = "all", features = rownames(d$Y)), d$obs, blocked())$tests[1, ]
  testthat::expect_identical(row$permutation_scheme, "between_subjects:subject")
  testthat::expect_identical(row$permdisp_scheme, "between_subjects:subject")
  testthat::expect_lt(row$p_value, 0.05)                       # was P = 1 with within-subject shuffling
  testthat::expect_false(identical(row$interpretation, "dispersion_difference_location_not_established"))
  n0 <- fn("pm_feature_set")(make_subjects("between_null")$Y, list(id = "all", features = rownames(d$Y)), make_subjects("between_null")$obs, blocked())$tests[1, ]
  testthat::expect_gt(n0$p_value, 0.2)
  # the pipeline's own permutation design: every permutation moves whole subjects and is not the identity
  how <- fn("pm_group_scheme")(d$obs$group, d$obs$subject, 999L, "subject")
  testthat::expect_identical(how$scheme, "between_subjects:subject")
  set.seed(1); P <- permute::shuffleSet(nrow(d$obs), 200, control = how$how)
  testthat::expect_true(all(apply(P, 1L, function(i) all(tapply(d$obs$group[i], d$obs$subject, function(v) length(unique(v))) == 1))))
  testthat::expect_lt(mean(apply(P, 1L, function(i) all(d$obs$group[i] == d$obs$group))), 0.1)
  # unbalanced subjects cannot be permuted as whole plots: typed refusal
  ub <- d; ub$obs <- ub$obs[-1, ]; ub$Y <- ub$Y[, -1]
  testthat::expect_error(fn("pm_feature_set")(ub$Y, list(id = "all", features = rownames(ub$Y)), ub$obs, blocked()), "E_PERMANOVA_BLOCKING_UNBALANCED")
})

testthat::test_that("audit M2: group varying within every subject permutes within subject (known answer)", {
  skip_vegan()
  d <- make_subjects("within")
  row <- fn("pm_feature_set")(d$Y, list(id = "all", features = rownames(d$Y)), d$obs, blocked())$tests[1, ]
  testthat::expect_identical(row$permutation_scheme, "within:subject"); testthat::expect_identical(row$permdisp_scheme, "within:subject")
  testthat::expect_lt(row$p_value, 0.05)
  n0 <- make_subjects("within_null")
  testthat::expect_gt(fn("pm_feature_set")(n0$Y, list(id = "all", features = rownames(n0$Y)), n0$obs, blocked())$tests[1, ]$p_value, 0.2)
  how <- fn("pm_group_scheme")(d$obs$group, d$obs$subject, 999L, "subject")
  set.seed(1); P <- permute::shuffleSet(nrow(d$obs), 100, control = how$how)
  testthat::expect_true(all(apply(P, 1L, function(i) all(d$obs$subject[i] == d$obs$subject))))
})

testthat::test_that("audit M2: mixed constant/varying subjects are refused with a typed error", {
  skip_vegan()
  d <- make_subjects("mixed")
  testthat::expect_error(fn("pm_feature_set")(d$Y, list(id = "all", features = rownames(d$Y)), d$obs, blocked()), "E_PERMANOVA_BLOCKING_MIXED")
})

testthat::test_that("audit m2: exact interpretation labels for constructed location-only and dispersion-only data", {
  skip_vegan()
  set.seed(9); n <- 8; p <- 10
  dev <- matrix(rnorm(n * p, 0, 0.3), p); dev <- dev - rowMeans(dev)
  Y <- cbind(10 + dev, 12 + dev); dimnames(Y) <- list(sprintf("F%02d", 1:p), c(paste0("A", 1:n), paste0("B", 1:n)))   # identical spread, shifted centroid
  obs <- data.frame(observation_id = colnames(Y), group = rep(c("A", "B"), each = n), stringsAsFactors = FALSE)
  row <- fn("pm_feature_set")(Y, list(id = "all", features = rownames(Y)), obs, settings_for(c("A", "B")))$tests[1, ]
  testthat::expect_identical(row$interpretation, "location_shift")
  testthat::expect_gt(row$permdisp_p, 0.5)
})

testthat::test_that("audit m2: covariate terms are tested with the pipeline's within-group permutations", {
  skip_vegan()
  d <- make_location(groups = c("A", "B"), n = 8); d$obs$sex <- rep(c("F", "M"), 8)
  how <- fn("pm_covariate_how")(d$obs$group, 199L)
  set.seed(1); P <- permute::shuffleSet(16, 100, control = how)
  testthat::expect_true(all(apply(P, 1L, function(i) all(d$obs$group[i] == d$obs$group))))
  testthat::expect_true(any(apply(P, 1L, function(i) any(d$obs$sex[i] != d$obs$sex))))
})

testthat::test_that("audit m3: empty declared groups and refused pairs are recorded; Holm family size is stated", {
  skip_vegan()
  d <- make_location(groups = c("A", "B", "C"))
  res <- fn("pm_feature_set")(d$Y, list(id = "all", features = rownames(d$Y)), d$obs, settings_for(c("A", "B", "C", "D")))
  pw <- res$tests[res$tests$analysis == "pairwise", ]
  testthat::expect_equal(nrow(pw), 3L)
  testthat::expect_true(all(pw$adjustment_family_size == 3L)); testthat::expect_true(all(pw$adjustment_family_planned == 6L))
  testthat::expect_true("group_without_observations" %in% res$refusals$reason)
  testthat::expect_setequal(res$refusals$comparison[res$refusals$analysis == "pairwise"], c("D vs A", "D vs B", "D vs C"))
  testthat::expect_true(any(res$refusals$analysis == "group" & res$refusals$comparison == "D"))
})

testthat::test_that("audit m4: an empty interaction cell refuses only the interaction term", {
  skip_vegan()
  d <- make_location(groups = c("A", "B"), n = 8); d$obs$sex <- rep(c("F", "M"), 8); d$obs$sex[d$g == "B"] <- "M"
  res <- fn("pm_feature_set")(d$Y, list(id = "all", features = rownames(d$Y)), d$obs,
                              settings_for(c("A", "B"), covariates = list(list(column = "sex", type = "categorical")), interaction = TRUE))
  testthat::expect_true(all(c("global", "group_adjusted", "covariate") %in% res$tests$analysis))
  testthat::expect_false("interaction" %in% res$tests$analysis)
  testthat::expect_identical(res$refusals$reason_code[res$refusals$analysis == "interaction"], "E_PERMANOVA_INTERACTION_NONESTIMABLE")
})

testthat::test_that("audit m6: non-syntactic covariate names are used safely", {
  skip_vegan()
  d <- make_location(groups = c("A", "B"), n = 8); d$obs[["age (years)"]] <- seq(20, 50, length.out = 16)
  res <- fn("pm_feature_set")(d$Y, list(id = "all", features = rownames(d$Y)), d$obs,
                              settings_for(c("A", "B"), covariates = list(list(column = "age (years)", type = "continuous"))))
  cov <- res$tests[res$tests$analysis == "covariate", ]
  testthat::expect_identical(cov$term, "age (years)")
  X <- scale(t(d$Y)); data <- data.frame(group = factor(d$g), age = d$obs[["age (years)"]])
  set.seed(7); direct <- vegan::adonis2(stats::dist(X) ~ group + age, data = data, permutations = permute::how(nperm = 199, blocks = data$group), by = "margin")
  testthat::expect_equal(cov$ss, direct["age", "SumOfSqs"], tolerance = 1e-10)
})

testthat::test_that("audit m5: the R side verifies the DEP table hash like the Python side", {
  dir <- tempfile(); dir.create(dir)
  plan <- list(plan_hash = "h1", artifacts = list())
  jsonlite::write_json(plan, file.path(dir, "plan.json"), auto_unbox = TRUE)
  writeLines(c("feature_id\tq_value", "F01\t0.01"), file.path(dir, "zero_null.tsv"))
  req <- list(plan_hash = "h1", inputs = list(list(artifact_id = "plan", path = file.path(dir, "plan.json"), sha256 = fn("sha256_file")(file.path(dir, "plan.json"))),
                                               list(artifact_id = "dea_zero_null", path = file.path(dir, "zero_null.tsv"), sha256 = strrep("0", 64))))
  testthat::expect_error(fn(".pm_verify_inputs")(req), "E_INTEGRITY: differential table")
  req$inputs[[2]]$sha256 <- fn("sha256_file")(file.path(dir, "zero_null.tsv"))
  testthat::expect_silent(fn(".pm_verify_inputs")(req))
})
