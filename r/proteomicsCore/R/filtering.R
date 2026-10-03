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
    # vectorized over features (audit 2026-10-02; same rules and outputs as the former per-feature loop)
    fraction <- sweep(counts, 2L, ifelse(n_units > 0, n_units, NA_real_), "/")
    features <- rownames(observed); nf <- length(features)
    excluded_extra <- if (is.null(extra_exclusions)) rep(FALSE, nf) else features %in% names(extra_exclusions)
    names_list <- function(mask) apply(mask, 1L, function(m) paste(required[m], collapse = ","))
    zero_units <- any(n_units == 0)
    no_obs <- counts == 0
    low <- (counts < min_n) | (fraction < min_fraction); low[is.na(low)] <- FALSE
    eligibility <- rep("eligible", nf); reason <- rep("eligible", nf)
    set <- function(mask, e, r) { mask <- mask & eligibility == "eligible" & reason == "eligible" & !decided; eligibility[mask] <<- e; reason[mask] <<- r[mask]; decided[mask] <<- TRUE }
    decided <- rep(FALSE, nf)
    if (policy == "native_dropout") set(rep(TRUE, nf), "not_run", rep("native_dropout_policy_deferred_to_R06", nf))
    if (any(excluded_extra)) { r <- rep(NA_character_, nf); r[excluded_extra] <- unlist(extra_exclusions[features[excluded_extra]]); set(excluded_extra, "excluded", r) }
    set(all_missing, "excluded", rep("all_study_missing", nf))
    if (zero_units) set(rep(TRUE, nf), "nonestimable", rep(paste0("required_group_has_no_units:", paste(required[n_units == 0], collapse = ",")), nf))
    any_no <- rowSums(no_obs) > 0
    if (any(any_no & !decided)) set(any_no, "nonestimable", paste0("all_missing_required_group:", names_list(no_obs)))
    any_low <- rowSums(low) > 0
    if (any(any_low & !decided)) set(any_low, "excluded", paste0("coverage_below_minimum:", names_list(low)))
    counts_json <- vapply(seq_len(nf), function(i) as.character(jsonlite::toJSON(as.list(stats::setNames(as.integer(counts[i, ]), required)), auto_unbox = TRUE)), "")
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(model_id = rep(rule$model_id, nf), contrast_id = rep(rule$contrast_id, nf), feature_id = features,
      eligibility = eligibility, reason = reason, n_obs_by_required_group = counts_json,
      policy = rep(policy, nf), minimum_observed_per_group = rep(min_n, nf), minimum_fraction = rep(min_fraction, nf), stringsAsFactors = FALSE, row.names = NULL)
    group_rows[[length(group_rows) + 1L]] <- data.frame(model_id = rule$model_id, contrast_id = rule$contrast_id,
      feature_id = rep(rownames(observed), times = length(required)), group = rep(required, each = nrow(observed)),
      n_observed = as.integer(as.vector(counts)), n_units = rep(as.integer(n_units), each = nrow(observed)),
      observed_fraction = as.vector(counts) / rep(n_units, each = nrow(observed)), stringsAsFactors = FALSE)
  }
  list(summary = if (length(summary_rows)) do.call(rbind, summary_rows) else data.frame(),
       by_group = if (length(group_rows)) do.call(rbind, group_rows) else data.frame())
}
