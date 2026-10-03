# DEqMS adapter (packet R06, SM10, V052-V054/V057).
#
# Official sequence on the exact contrast-as-coefficient refit:
#   lmFit -> eBayes -> fit$count <- counts -> spectraCounteBayes -> outputResult(coef_col = 1).
# Rows report sca.P.Value as the DEqMS P value; the original limma fields are
# kept side by side under native_limma_* names, never substituted.

.assay_base_rows <- function(model, contrast, features, est, context, engine, version) {
  out <- endpoint_frame(length(features))
  out$feature_id <- features; out$model_id <- model$model_id; out$design_id <- model$design_id; out$contrast_id <- contrast$contrast_id
  out$run_id <- context$run_id; out$plan_hash <- context$plan_hash; out$engine <- engine; out$engine_version <- version
  out$result_type <- "ProteinZeroNullResult"; out$hypothesis_type <- "protein_zero_null"; out$role <- "secondary"
  e <- est[est$model_id == model$model_id & est$contrast_id == contrast$contrast_id, , drop = FALSE]; e <- e[match(features, e$feature_id), , drop = FALSE]
  out$estimable <- e$estimable == "true"; out$n_obs_by_required_group <- e$n_obs_by_required_group
  list(rows = out, est = e)
}

fit_deqms_model <- function(model, values, X, contrasts, estimability, counts, context) {
  version <- as.character(utils::packageVersion("DEqMS"))
  rows <- list(); curves <- list()
  for (contrast in contrasts) {
    base <- .assay_base_rows(model, contrast, rownames(values), estimability, context, "deqms", version)
    out <- base$rows; e <- base$est
    eligible <- e$eligibility == "eligible"
    out$eligibility <- ifelse(eligible, "tested", e$eligibility); out$reason_code <- ifelse(eligible, NA_character_, e$reason)
    U <- rownames(values)[eligible]
    ex <- exact_contrast_fit(values[U, , drop = FALSE], X, as.numeric(unlist(contrast$weights)))
    fit <- limma::eBayes(ex$fit)
    fit$count <- counts$fit[rownames(fit$coefficients)]
    fit4 <- DEqMS::spectraCounteBayes(fit, coef_col = 1)
    res <- DEqMS::outputResult(fit4, coef_col = 1)
    res <- res[match(U, rownames(res)), , drop = FALSE]
    idx <- match(U, out$feature_id)
    sca_se <- unname(fit4$stdev.unscaled[U, 1] * sqrt(fit4$sca.postvar[U]))
    sca_df <- unname(fit4$sca.dfprior + fit4$df.residual[match(U, rownames(fit4$coefficients))])
    out$effect[idx] <- res$logFC; out$effect_se[idx] <- sca_se; out$statistic[idx] <- res$sca.t; out$statistic_type[idx] <- "deqms_sca_t"
    out$p_value[idx] <- res$sca.P.Value; out$df_inference[idx] <- sca_df; out$df_residual[idx] <- fit4$df.residual[match(U, rownames(fit4$coefficients))]
    ci <- ci_moderated_t(res$logFC, sca_se, sca_df, 0.95)
    out$ci_lower[idx] <- ci[, "lower"]; out$ci_upper[idx] <- ci[, "upper"]; out$ci_level[idx] <- 0.95
    out$ci_method[idx] <- "logFC +/- qt(.975, sca.dfprior + df.residual) * stdev.unscaled * sqrt(sca.postvar)"
    out$null_region[idx] <- "log2 effect = 0"; out$mean_abundance[idx] <- res$AveExpr; out$exactness_path[idx] <- ex$path
    out$native_q_value[idx] <- res$sca.adj.pval
    out$native_limma_p_value <- NA_real_; out$native_limma_t <- NA_real_; out$native_sca_postvar <- NA_real_; out$count_original <- NA_real_; out$count_fit <- NA_real_
    out$native_limma_p_value[idx] <- res$P.Value; out$native_limma_t[idx] <- res$t; out$native_sca_postvar[idx] <- fit4$sca.postvar[U]
    out$count_original[idx] <- counts$original[U]; out$count_fit[idx] <- counts$fit[U]
    rows[[length(rows) + 1L]] <- out
    curves[[length(curves) + 1L]] <- data.frame(model_id = model$model_id, contrast_id = contrast$contrast_id, feature_id = U, log2_count = log2(counts$fit[U]),
                                                log_residual_variance = log(fit4$sigma[match(U, rownames(fit4$coefficients))]^2), sca_priorvar = fit4$sca.priorvar[U], stringsAsFactors = FALSE)
  }
  list(rows = do.call(rbind, rows), curve = do.call(rbind, curves),
       diagnostics = list(engine = "DEqMS", version = version, sequence = "lmFit(exact contrast refit) -> eBayes -> spectraCounteBayes(fit.method='loess') -> outputResult(coef_col=1)",
                          count_zero_policy = counts$policy, count_offset = counts$offset))
}
