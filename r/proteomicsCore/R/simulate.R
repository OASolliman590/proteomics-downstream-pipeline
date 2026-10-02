# Release calibration simulator (packet R11, V105).  Uses the production
# adapter functions (exact contrast refit, eBayes with the frozen settings,
# moderated-t intervals and the central BH API) on independently simulated
# datasets.  The dataset, not the protein, is the replication unit.

simulate_dataset <- function(n_features, n_per_group, scenario, seed, effect = 1, prop_alt = 0.1) {
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  g <- rep(c(0, 1), each = n_per_group)
  sd <- sqrt(1 / stats::rchisq(n_features, df = 6) * 6) * 0.4                 # feature-specific variances (moderation is meaningful)
  noise <- switch(scenario,
                  heavy_tail = matrix(stats::rt(n_features * length(g), df = 3), n_features) / sqrt(3),
                  matrix(stats::rnorm(n_features * length(g)), n_features))
  true <- rep(0, n_features)
  if (scenario %in% c("mixture", "heavy_tail", "mnar")) { alt <- sample.int(n_features, round(prop_alt * n_features)); true[alt] <- sample(c(-1, 1), length(alt), TRUE) * effect * sd[alt] * 2 }
  Y <- 20 + noise * sd + outer(true, g)
  if (scenario == "mnar") { cut <- stats::quantile(Y, 0.08); Y[Y < cut & stats::runif(length(Y)) < 0.8] <- NA }
  rownames(Y) <- sprintf("F%05d", seq_len(n_features))
  list(Y = Y, X = cbind(A = 1 - g, B = g), true = true)
}

calibrate_once <- function(sim, alpha = 0.05) {
  keep <- rowSums(!is.na(sim$Y[, sim$X[, 1] == 1])) >= 2 & rowSums(!is.na(sim$Y[, sim$X[, 2] == 1])) >= 2
  ex <- exact_contrast_fit(sim$Y[keep, , drop = FALSE], sim$X, c(-1, 1))
  eb <- limma::eBayes(ex$fit, trend = TRUE, robust = TRUE)
  p <- eb$p.value[, 1]; q <- adjust_pvalues(p, "BH")
  se <- sqrt(eb$s2.post) * eb$stdev.unscaled[, 1]
  ci <- ci_moderated_t(eb$coefficients[, 1], se, eb$df.total, 0.95)
  truth <- sim$true[keep]
  rejected <- q <= alpha
  list(any_rejection = any(rejected), fdp = if (any(rejected)) mean(truth[rejected] == 0) else 0,
       power = if (any(truth != 0)) mean(rejected[truth != 0]) else NA_real_, coverage = mean(ci[, "lower"] <= truth & truth <= ci[, "upper"]),
       bias = mean(eb$coefficients[truth != 0, 1] - truth[truth != 0]), n_tested = sum(keep))
}

run_core_calibration <- function(scenario, n_datasets, n_features, n_per_group, seed_base) {
  rows <- lapply(seq_len(n_datasets), function(i) { r <- calibrate_once(simulate_dataset(n_features, n_per_group, scenario, seed_base + i)); r$dataset <- i; r })
  data.frame(scenario = scenario, dataset = vapply(rows, `[[`, 0, "dataset"), any_rejection = vapply(rows, `[[`, TRUE, "any_rejection"),
             fdp = vapply(rows, `[[`, 0, "fdp"), power = vapply(rows, `[[`, 0, "power"), coverage = vapply(rows, `[[`, 0, "coverage"),
             bias = vapply(rows, function(r) if (is.null(r$bias) || is.nan(r$bias)) NA_real_ else r$bias, 0), n_tested = vapply(rows, `[[`, 0, "n_tested"))
}
