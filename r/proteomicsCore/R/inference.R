# Typed endpoint rows: zero-null, TREAT and omnibus (packet R05, SM11/SM13).

ci_moderated_t <- function(effect, se, df, level) {
  quantile <- stats::qt(1 - (1 - level) / 2, df)
  cbind(lower = effect - quantile * se, upper = effect + quantile * se)
}

# Exact GLS moderated F for H0: beta_S = 0 (per feature, observed rows).
omnibus_moderated_f <- function(Y, X, coefficient_index, eb_full, weights = NULL, block = NULL, correlation = NULL) {
  k <- length(coefficient_index)
  out <- data.frame(feature_id = rownames(Y), statistic = NA_real_, df1 = k, df2 = NA_real_, p_value = NA_real_, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(Y))) {
    use <- !is.na(Y[i, ])
    Xi <- X[use, , drop = FALSE]
    if (qr(Xi)$rank < ncol(X)) next
    n <- sum(use)
    V <- diag(n)
    if (!is.null(block)) { b <- block[use]; same <- outer(b, b, "=="); V[same] <- correlation; diag(V) <- 1 }
    wi <- if (is.null(weights)) rep(1, n) else if (is.matrix(weights)) weights[i, use] else weights[use]
    Wh <- diag(sqrt(wi), n)
    unscaled <- solve(t(Xi) %*% Wh %*% solve(V) %*% Wh %*% Xi)
    beta <- eb_full$coefficients[i, coefficient_index]
    S <- unscaled[coefficient_index, coefficient_index, drop = FALSE]
    statistic <- as.numeric(t(beta) %*% solve(S) %*% beta) / (k * eb_full$s2.post[i])
    # limma convention for moderated F: df2 = df.prior + df.residual (classifyTestsF), not the capped df.total
    df_prior <- if (length(eb_full$df.prior) == 1L) eb_full$df.prior else eb_full$df.prior[i]
    df2 <- df_prior + eb_full$df.residual[i]
    out$statistic[i] <- statistic; out$df2[i] <- df2
    out$p_value[i] <- stats::pf(statistic, k, df2, lower.tail = FALSE)
  }
  out
}

endpoint_frame <- function(n) {
  data.frame(schema_version = rep("1.2.0", n), run_id = NA_character_, plan_hash = NA_character_, result_type = NA_character_, hypothesis_type = NA_character_,
             family_id = NA_character_, engine = "limma", engine_version = NA_character_, model_id = NA_character_, design_id = NA_character_,
             contrast_id = NA_character_, feature_id = NA_character_, estimable = NA, n_obs_by_required_group = NA_character_, eligibility = NA_character_,
             reason_code = NA_character_, effect = NA_real_, effect_scale = "log2_difference", effect_se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
             ci_level = NA_real_, ci_method = NA_character_, statistic = NA_real_, statistic_type = NA_character_, df_residual = NA_real_, df_inference = NA_real_,
             p_value = NA_real_, q_value = NA_real_, native_q_value = NA_real_, effect_threshold = NA_real_, null_region = NA_character_,
             mean_abundance = NA_real_, role = NA_character_, display_effect_filter = NA, exactness_path = NA_character_, stringsAsFactors = FALSE)
}
