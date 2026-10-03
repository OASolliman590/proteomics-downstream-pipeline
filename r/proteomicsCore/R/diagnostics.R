# Moderation diagnostics (packet R05, SM09, V043).
#
# Records requested and actual trend/robust settings and the prior and
# posterior variances.  robust=TRUE is robust empirical-Bayes hyperparameter
# estimation (limma eBayes), not robust sample regression.

limma_preflight <- function(robust, installed = function(pkg) requireNamespace(pkg, quietly = TRUE)) {
  if (!installed("limma")) stop("E_DEPENDENCY_UNAVAILABLE: limma is not installed", call. = FALSE)
  if (isTRUE(robust) && !installed("statmod")) stop("E_DEPENDENCY_UNAVAILABLE: robust=TRUE requires statmod; robust is never silently switched off", call. = FALSE)
  invisible(TRUE)
}

moderation_diagnostics <- function(eb, model_id, trend, robust) {
  prior <- if (length(eb$s2.prior) == 1L) rep(eb$s2.prior, nrow(eb$coefficients)) else eb$s2.prior
  df_prior <- if (length(eb$df.prior) == 1L) rep(eb$df.prior, nrow(eb$coefficients)) else eb$df.prior
  data.frame(model_id = model_id, feature_id = rownames(eb$coefficients), mean_abundance = eb$Amean, sigma = eb$sigma, s2_residual = eb$sigma^2,
             df_residual = eb$df.residual, s2_prior = prior, df_prior = df_prior, s2_posterior = eb$s2.post, df_total = eb$df.total,
             trend_requested = trend, robust_requested = robust, stringsAsFactors = FALSE)
}

moderation_settings <- function(eb, trend, robust) {
  list(trend_requested = trend, robust_requested = robust, trend_applied = trend, robust_applied = robust,
       robust_meaning = "robust empirical-Bayes hyperparameter estimation (limma::eBayes robust=TRUE); not robust sample regression",
       df_prior_summary = .pc_summary_list(eb$df.prior), s2_prior_summary = .pc_summary_list(eb$s2.prior),
       limma_version = as.character(utils::packageVersion("limma")),
       statmod_version = if (requireNamespace("statmod", quietly = TRUE)) as.character(utils::packageVersion("statmod")) else NULL)
}

.pc_summary_list <- function(x) { s <- summary(as.numeric(x)); as.list(stats::setNames(as.numeric(s), names(s))) }
