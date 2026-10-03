# Release calibration simulator (packet R11, V105).  Each simulated dataset is
# analysed by the production path: R03 coverage_tables (available-case rule,
# >= 2 observed and >= 50 % per required group), R04 featurewise_estimability,
# R05 fit_limma_model (exact contrast refit, eBayes trend = TRUE, robust = TRUE,
# the example configurations' settings, moderated-t intervals) and the central
# family adjustment adjust_family (BH over finite eligible tests).  The dataset,
# not the protein, is the replication unit.  (Audit 2026-10-02: the former
# version re-implemented the fit and a simplified eligibility rule.)

simulate_dataset <- function(n_features, n_per_group, scenario, seed, effect = if (identical(scenario, "mixture_high_power")) 2 else 1,
                             prop_alt = if (identical(scenario, "mixture_high_power")) 0.2 else 0.1) {
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  g <- rep(c(0, 1), each = n_per_group)
  sd <- sqrt(1 / stats::rchisq(n_features, df = 6) * 6) * 0.4                 # feature-specific variances (moderation is meaningful)
  noise <- switch(scenario,
                  heavy_tail = matrix(stats::rt(n_features * length(g), df = 3), n_features) / sqrt(3),
                  matrix(stats::rnorm(n_features * length(g)), n_features))
  true <- rep(0, n_features)
  if (scenario %in% c("mixture", "mixture_high_power", "heavy_tail", "mnar")) { alt <- sample.int(n_features, round(prop_alt * n_features)); true[alt] <- sample(c(-1, 1), length(alt), TRUE) * effect * sd[alt] * 2 }
  Y <- 20 + noise * sd + outer(true, g)
  if (scenario == "mnar") { cut <- stats::quantile(Y, 0.08); Y[Y < cut & stats::runif(length(Y)) < 0.8] <- NA }
  rownames(Y) <- sprintf("F%05d", seq_len(n_features))
  list(Y = Y, X = cbind(A = 1 - g, B = g), true = true)
}

calibrate_once <- function(sim, alpha = 0.05) {
  g <- ifelse(sim$X[, 2] == 1, "B", "A")
  ids <- sprintf("O%03d", seq_along(g))
  Y <- sim$Y; colnames(Y) <- ids
  X <- sim$X; rownames(X) <- ids
  observed <- !is.na(Y)
  observations <- data.frame(observation_id = ids, group = g, stringsAsFactors = FALSE)
  rule <- list(model_id = "cal", contrast_id = "B-A", required_groups = list("B", "A"), minimum_observed_per_group = 2L, minimum_fraction = 0.5, policy = "available_case")
  coverage <- coverage_tables(observed, observations, "group", list(rule))$summary
  contrast <- list(contrast_id = "B-A", design_id = "joint", role = "primary", weights = list(-1, 1), required_groups = list("B", "A"))
  estimability <- featurewise_estimability(Y, observed, X, list(contrast), coverage, "cal", g)
  model <- list(model_id = "cal", design_id = "joint", role = "primary", trend = TRUE, robust = TRUE, hypotheses = list("zero_null"), effect_threshold = NULL)
  fit <- fit_limma_model(model, Y, X, list(contrast), estimability, observations, NULL, NULL, NULL, 0.95, list(run_id = "calibration", plan_hash = "calibration"))
  family <- list(family_id = "cal", hypothesis_type = "protein_zero_null", role = "primary", adjustment = "BH", dependence_assumption = "PRDS",
                 q_cutoff = alpha, members = list(list(model_id = "cal", contrast_id = "B-A")))
  rows <- adjust_family(fit$rows, family)$rows
  tested <- rows[rows$eligibility == "tested", , drop = FALSE]
  truth <- sim$true[match(tested$feature_id, rownames(sim$Y))]
  rejected <- tested$q_value <= alpha
  list(any_rejection = any(rejected), fdp = if (any(rejected)) mean(truth[rejected] == 0) else 0,
       power = if (any(truth != 0)) mean(rejected[truth != 0]) else NA_real_, coverage = mean(tested$ci_lower <= truth & truth <= tested$ci_upper),
       bias = mean(tested$effect[truth != 0] - truth[truth != 0]), n_tested = nrow(tested), n_rejected = sum(rejected), fit_path = "production:fit_limma_model")
}

run_core_calibration <- function(scenario, n_datasets, n_features, n_per_group, seed_base) {
  rows <- lapply(seq_len(n_datasets), function(i) { r <- calibrate_once(simulate_dataset(n_features, n_per_group, scenario, seed_base + i)); r$dataset <- i; r })
  data.frame(scenario = scenario, dataset = vapply(rows, `[[`, 0, "dataset"), any_rejection = vapply(rows, `[[`, TRUE, "any_rejection"),
             fdp = vapply(rows, `[[`, 0, "fdp"), power = vapply(rows, `[[`, 0, "power"), coverage = vapply(rows, `[[`, 0, "coverage"),
             bias = vapply(rows, function(r) if (is.null(r$bias) || is.nan(r$bias)) NA_real_ else r$bias, 0), n_tested = vapply(rows, `[[`, 0, "n_tested"),
             n_rejected = vapply(rows, `[[`, 0, "n_rejected"))
}
