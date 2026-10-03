# Typed differential-table verification before export (packet R05, V048).
#
# Every planned endpoint row must be present with the shared result-row
# envelope; only tested rows carry statistics; no NA is coerced to zero; a
# central q value only exists for rows inside a declared family.

DIFFERENTIAL_FIELDS <- c("schema_version", "run_id", "plan_hash", "result_type", "hypothesis_type", "family_id", "engine", "engine_version",
                         "model_id", "design_id", "contrast_id", "feature_id", "estimable", "n_obs_by_required_group", "eligibility", "reason_code",
                         "effect", "effect_scale", "effect_se", "ci_lower", "ci_upper", "ci_level", "ci_method", "statistic", "statistic_type",
                         "df_residual", "df_inference", "p_value", "q_value", "native_q_value", "effect_threshold", "null_region", "mean_abundance", "role")
ELIGIBILITY_STATES <- c("tested", "excluded", "nonestimable", "inapplicable", "not_run", "numerical_failure")

verify_differential_table <- function(rows, planned_keys) {
  missing_fields <- setdiff(DIFFERENTIAL_FIELDS, names(rows))
  if (length(missing_fields)) stop(sprintf("E_TABLE_SCHEMA: missing fields %s", paste(missing_fields, collapse = ", ")), call. = FALSE)
  if (!all(rows$eligibility %in% ELIGIBILITY_STATES)) stop("E_TABLE_SCHEMA: unknown eligibility state", call. = FALSE)
  keys <- paste(rows$model_id, rows$contrast_id, rows$feature_id, rows$hypothesis_type, sep = "|")
  if (anyDuplicated(keys)) stop("E_TABLE_SCHEMA: duplicated endpoint rows", call. = FALSE)
  absent <- setdiff(planned_keys, keys)
  if (length(absent)) stop(sprintf("E_TABLE_INCOMPLETE: %d planned endpoint row(s) are missing", length(absent)), call. = FALSE)
  untested <- rows$eligibility != "tested"
  if (any(!is.na(rows$p_value[untested])) || any(!is.na(rows$q_value[untested]))) stop("E_TABLE_SCHEMA: untested rows carry P/q values", call. = FALSE)
  if (any(is.na(rows$reason_code[untested]))) stop("E_TABLE_SCHEMA: untested rows need a reason_code", call. = FALSE)
  if (any(!is.na(rows$q_value) & is.na(rows$family_id))) stop("E_TABLE_SCHEMA: central q without a declared family", call. = FALSE)
  tested <- !untested
  if (any(!is.finite(rows$p_value[tested]))) stop("E_TABLE_SCHEMA: tested row without a finite P value", call. = FALSE)
  invisible(TRUE)
}
