# proDA adapter (packet R06, SM10, V055-V057).
#
# Native proDA fit on unimputed LFQ data with genuine NA, independent fixed
# design; test_diff gives the native zero-null result.  A target group with
# no observations follows proDA's native dropout model: n_obs and the
# prior/dropout dependence are retained, never invented zeros.

fit_proda_model <- function(model, values, X, contrasts, estimability, context, obs, group_column, seed) {
  version <- as.character(utils::packageVersion("proDA"))
  keep <- rowSums(!is.na(values)) > 0
  set.seed(seed, kind = "L'Ecuyer-CMRG")   # proDA fitting uses random numbers; the recorded seed makes it reproducible
  fit <- proDA::proDA(values[keep, , drop = FALSE], design = X)
  conv <- tryCatch(fit$convergence, error = function(e) NULL)
  if (!is.null(conv) && !isTRUE(conv$successful)) stop("E_ENGINE_FAILED: proDA did not converge", call. = FALSE)
  rows <- list()
  groups <- obs[[group_column]]
  for (contrast in contrasts) {
    base <- .assay_base_rows(model, contrast, rownames(values), estimability, context, "proda", version)
    out <- base$rows
    td <- proDA::test_diff(fit, contrast = as.numeric(unlist(contrast$weights)))
    td <- td[match(rownames(values), td$name), , drop = FALSE]
    tested <- keep & !is.na(td$pval)
    required <- unlist(contrast$required_groups)
    target_missing <- vapply(seq_len(nrow(values)), function(i) any(vapply(required, function(g) all(is.na(values[i, groups == g])), logical(1))), logical(1))
    out$eligibility <- ifelse(tested, "tested", "excluded"); out$reason_code <- ifelse(tested, NA_character_, ifelse(keep, "proda_no_result", "all_study_missing"))
    out$effect <- ifelse(tested, td$diff, NA); out$effect_se <- ifelse(tested, td$se, NA); out$statistic <- ifelse(tested, td$t_statistic, NA)
    out$statistic_type <- ifelse(tested, "proda_t", NA); out$p_value <- ifelse(tested, td$pval, NA); out$df_inference <- ifelse(tested, td$df, NA)
    out$native_q_value <- ifelse(tested, td$adj_pval, NA); out$mean_abundance <- ifelse(tested, td$avg_abundance, NA)
    out$null_region <- ifelse(tested, "log2 effect = 0", NA)
    ci <- ci_moderated_t(td$diff, td$se, td$df, 0.95)
    out$ci_lower <- ifelse(tested, ci[, "lower"], NA); out$ci_upper <- ifelse(tested, ci[, "upper"], NA); out$ci_level <- ifelse(tested, 0.95, NA)
    out$ci_method <- ifelse(tested, "diff +/- qt(.975, native df) * native se (proDA test_diff)", NA)
    out$estimable <- tested
    out$native_n_obs <- td$n_obs; out$native_n_approx <- td$n_approx
    out$dropout_prior_dependent <- tested & target_missing
    rows[[length(rows) + 1L]] <- out
  }
  hp <- tryCatch(proDA::hyper_parameters(fit), error = function(e) NULL)
  curve <- if (!is.null(hp) && !is.null(hp$dropout_curve_position)) data.frame(model_id = model$model_id, observation_id = colnames(values),
             dropout_curve_position = hp$dropout_curve_position, dropout_curve_scale = hp$dropout_curve_scale, stringsAsFactors = FALSE) else NULL
  list(rows = do.call(rbind, rows), curve = curve,
       diagnostics = list(engine = "proDA", version = version, convergence = conv, seed = seed, rng_kind = "L'Ecuyer-CMRG", call = "proDA::proDA(data, design = numeric design) then test_diff(contrast = weights)",
                          note = "estimates for all-missing target groups depend on the dropout prior (dropout_prior_dependent = TRUE)"))
}
