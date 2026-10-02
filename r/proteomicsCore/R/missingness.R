# Missingness diagnostics (packet R03, V025).
#
# Separates numeric availability, genuine original observation and prior
# imputation.  Summaries are descriptive; no MCAR/MAR/MNAR mechanism is
# asserted from any cutoff.

missingness_tables <- function(values, observed, prior_imputed, observations, group_column) {
  numeric <- !is.na(values)
  known <- !is.null(observed)
  by_obs <- data.frame(observation_id = colnames(values), group = observations[[group_column]], n_features = nrow(values),
                       n_numeric = colSums(numeric), n_numeric_missing = colSums(!numeric),
                       n_original_observed = if (known) colSums(observed) else NA_integer_,
                       n_prior_imputed = if (known) colSums(prior_imputed) else NA_integer_,
                       original_mask_state = if (known) "known" else "unknown", stringsAsFactors = FALSE)
  groups <- unique(observations[[group_column]])
  by_group <- do.call(rbind, lapply(groups, function(g) {
    cols <- observations[[group_column]] == g
    data.frame(group = g, n_units = sum(cols), n_cells = nrow(values) * sum(cols), n_numeric = sum(numeric[, cols]),
               n_original_observed = if (known) sum(observed[, cols]) else NA_integer_,
               n_prior_imputed = if (known) sum(prior_imputed[, cols]) else NA_integer_, stringsAsFactors = FALSE)
  }))
  by_feature <- data.frame(feature_id = rownames(values), n_numeric = rowSums(numeric),
                           n_original_observed = if (known) rowSums(observed) else NA_integer_,
                           n_prior_imputed = if (known) rowSums(prior_imputed) else NA_integer_,
                           mean_observed_abundance = vapply(seq_len(nrow(values)), function(i) { use <- if (known) observed[i, ] else numeric[i, ]; if (any(use)) mean(values[i, use]) else NA_real_ }, numeric(1)),
                           stringsAsFactors = FALSE)
  base <- if (known) observed else numeric
  by_feature$missing_fraction <- 1 - rowSums(base) / ncol(values)
  finite <- !is.na(by_feature$mean_observed_abundance)
  abundance <- data.frame()
  if (sum(finite) >= 2L) {
    k <- min(5L, sum(finite))
    breaks <- unique(stats::quantile(by_feature$mean_observed_abundance[finite], probs = seq(0, 1, length.out = k + 1L), names = FALSE))
    bins <- if (length(breaks) > 1L) cut(by_feature$mean_observed_abundance[finite], breaks, include.lowest = TRUE) else factor(rep("all", sum(finite)))
    abundance <- do.call(rbind, lapply(levels(bins), function(b) { sel <- bins == b; data.frame(abundance_bin = b, n_features = sum(sel), mean_missing_fraction = mean(by_feature$missing_fraction[finite][sel]), stringsAsFactors = FALSE) }))
  }
  list(by_observation = by_obs, by_group = by_group, by_feature = by_feature, abundance = abundance,
       statement = list(mask_state = if (known) "known" else "unknown", mechanism_assertion = "none",
                        note = "Missingness summaries are descriptive; no MCAR/MAR/MNAR mechanism is inferred from these tallies."))
}
