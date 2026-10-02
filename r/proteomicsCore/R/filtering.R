# Contrast-aware coverage masks (packet R03, SM04, V024).
#
# Eligibility counts genuinely observed biological units from the original
# observed mask, never numeric availability.  Ordinary available-case rule:
# at least minimum_observed_per_group observed units and at least
# minimum_fraction observed in every required group.  A required group with
# zero observed units is nonestimable for limma; all-study-missing rows are
# excluded.  Nothing is filled.

coverage_tables <- function(observed, observations, group_column, rules, extra_exclusions = NULL) {
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
