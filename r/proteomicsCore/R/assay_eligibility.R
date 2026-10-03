# Count evidence validation and the assay-engine stage (packet R06, SM10, V051/V054/V055/V059).

COUNT_TYPES <- c("peptide", "psm")
COUNT_PROXY_SOURCES <- c("observed_sample_count", "abundance_derived", "intensity_derived", "imputed")

# Accept only genuine declared peptide/PSM evidence, aligned by feature_id; preserve
# zeros/missing; require strictly positive fit counts after an explicit policy.
validate_count_evidence <- function(counts, features, aggregation, zero_policy) {
  for (column in c("feature_id", "count_type", "count_value", "count_source", "count_aggregation", "pseudocount_policy"))
    if (!column %in% names(counts)) stop(sprintf("E_DEQMS_COUNT_EVIDENCE: count evidence lacks column %s", column), call. = FALSE)
  if (any(!counts$count_type %in% COUNT_TYPES)) stop("E_DEQMS_COUNT_EVIDENCE: count_type must be peptide or psm", call. = FALSE)
  if (any(counts$count_source %in% COUNT_PROXY_SOURCES)) stop("E_DEQMS_COUNT_EVIDENCE: observed-sample counts or abundance-derived proxies are not count evidence", call. = FALSE)
  if (length(unique(counts$count_type)) != 1L) stop("E_DEQMS_COUNT_EVIDENCE: mixed peptide and PSM counts", call. = FALSE)
  declared <- unique(counts$count_aggregation)
  if (!identical(declared, aggregation)) stop(sprintf("E_DEQMS_COUNT_EVIDENCE: model count_aggregation %s disagrees with recorded %s", aggregation, paste(declared, collapse = ",")), call. = FALSE)
  value <- suppressWarnings(as.numeric(ifelse(counts$count_value == "NA", NA, counts$count_value)))
  if (any(!is.na(value) & (value < 0 | value != round(value)))) stop("E_DEQMS_COUNT_EVIDENCE: counts must be non-negative integers", call. = FALSE)
  key_cols <- intersect(c("feature_id", "observation_id", "plex_id"), names(counts))
  keys <- do.call(paste, c(counts[key_cols], sep = "\r"))
  if (anyDuplicated(keys)) stop("E_DEQMS_COUNT_EVIDENCE: duplicate count rows for the same key", call. = FALSE)
  per_feature <- if (identical(aggregation, "protein_level")) {
    if (anyDuplicated(counts$feature_id)) stop("E_DEQMS_COUNT_EVIDENCE: protein_level aggregation needs one count per feature", call. = FALSE)
    stats::setNames(value, counts$feature_id)
  } else tapply(value, counts$feature_id, function(v) if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE))
  missing <- setdiff(features, names(per_feature))
  if (length(missing)) stop(sprintf("E_DEQMS_COUNT_EVIDENCE: no count evidence for features %s", paste(utils::head(missing, 5), collapse = ", ")), call. = FALSE)
  original <- per_feature[features]
  fit_counts <- original
  if (identical(zero_policy$mode, "explicit_offset")) {
    if (is.null(zero_policy$offset) || is.null(zero_policy$justification) || !nzchar(zero_policy$justification)) stop("E_DEQMS_COUNT_NONPOSITIVE: explicit_offset needs a positive offset and a justification", call. = FALSE)
    fit_counts <- original + as.numeric(zero_policy$offset)
  }
  if (any(is.na(fit_counts) | fit_counts <= 0)) stop("E_DEQMS_COUNT_NONPOSITIVE: DEqMS fit counts must be finite and strictly positive; no hidden +1 is applied", call. = FALSE)
  list(original = original, fit = fit_counts, policy = zero_policy$mode, offset = if (identical(zero_policy$mode, "explicit_offset")) zero_policy$offset else 0)
}

# Scientific eligibility is decided before and independently of installation (SM10).
assay_scientific_eligibility <- function(model, design_blocking, weighted, assay, prior_imputation) {
  if (identical(model$engine, "deqms")) {
    if (weighted || !identical(design_blocking, "none")) return("E_DEQMS_DESIGN_UNSUPPORTED")
    if (!identical(unique(unlist(model$hypotheses)), "zero_null")) return("E_DEQMS_HYPOTHESIS_UNSUPPORTED")
  }
  if (identical(model$engine, "proda")) {
    if (!assay %in% c("lfq_dda", "lfq_dia") || !identical(prior_imputation, "none_documented")) return("E_PRODA_LFQ_UNIMPUTED_REQUIRED")
    if (weighted || !identical(design_blocking, "none")) return("E_PRODA_DESIGN_UNSUPPORTED")
    if (!identical(unique(unlist(model$hypotheses)), "zero_null")) return("E_PRODA_HYPOTHESIS_UNSUPPORTED")
  }
  NULL
}

assay_engine_stage <- function(request) .pc_run_stage(request, function(out) {
  p <- request$parameters
  outputs <- list(); warnings <- list()
  emit <- function(relative, id, type) outputs[[length(outputs) + 1L]] <<- .pc_output(out, relative, id, type)
  write_tsv <- function(df, relative, id, type) { .pc_write_tsv(df, file.path(out, relative)); emit(relative, id, type) }
  write_json <- function(value, relative, id, type) { .pc_write_json(value, file.path(out, relative)); emit(relative, id, type) }
  plan <- .pm_verify_inputs_response(request)
  context <- list(run_id = request$run_id, plan_hash = request$plan_hash)
  values <- .pc_matrix_from_tsv(.pc_find_input(request, "primary_matrix"), "numeric")
  observed <- .pc_matrix_from_tsv(.pc_find_input(request, "primary_observed_mask"), "logical")[rownames(values), colnames(values), drop = FALSE]
  values[!observed] <- NA
  obs <- .pc_read_tsv(.pc_find_input(request, "primary_observations")); obs <- obs[match(colnames(values), obs$observation_id), , drop = FALSE]
  estimability <- .pc_read_tsv(.pc_find_input(request, "estimability"))
  primary <- .pc_read_tsv(.pc_find_input(request, "dea_zero_null"))
  statuses <- list(); all_rows <- list(); comparisons <- list(); diagnostics <- list()
  for (model in p$models) {
    status <- tryCatch({
      reason <- assay_scientific_eligibility(model, model$blocking_mode, isTRUE(model$weighted), p$assay, p$prior_imputation)
      if (!is.null(reason)) stop(sprintf("%s: model %s is scientifically ineligible", reason, model$model_id), call. = FALSE)
      pkg <- if (identical(model$engine, "deqms")) "DEqMS" else "proDA"
      if (!requireNamespace(pkg, quietly = TRUE)) stop(sprintf("E_ENGINE_NOT_AVAILABLE: R package %s is not installed; the eligible %s request stays in the plan as NOT_RUN", pkg, model$model_id), call. = FALSE)
      X <- .pc_matrix_from_tsv(.pc_find_input(request, paste0("design_", model$design_id)), "numeric")
      contrasts <- Filter(function(c) identical(c$design_id, model$design_id), p$contrasts)
      res <- if (identical(model$engine, "deqms")) {
        counts <- validate_count_evidence(.pc_read_tsv(.pc_find_input(request, paste0("counts_", model$model_id))), rownames(values), model$count_aggregation, model$count_zero_policy)
        fit_deqms_model(model, values, X, contrasts, estimability, counts, context)
      } else fit_proda_model(model, values, X, contrasts, estimability, context, obs, p$group_column, as.integer(request$rng$seed))
      all_rows[[length(all_rows) + 1L]] <- res$rows
      diagnostics[[model$model_id]] <- res$diagnostics
      if (!is.null(res$curve)) write_tsv(res$curve, file.path("diagnostics", paste0(model$model_id, "_", if (identical(model$engine, "deqms")) "count_variance" else "dropout", ".tsv")), paste0("diagnostics_", model$model_id), "EngineDiagnostics")
      comparisons[[length(comparisons) + 1L]] <- compare_with_primary(res$rows, primary, p$primary_model_id)
      list(model_id = model$model_id, engine = model$engine, state = "COMPLETED", required = identical(model$execution_requirement, "required"), reason_code = NULL)
    }, error = function(e) {
      code <- regmatches(conditionMessage(e), regexpr("^E_[A-Z0-9_]+", conditionMessage(e))); if (!length(code)) code <- "E_ENGINE_FAILED"
      if (identical(model$execution_requirement, "required")) stop(if (startsWith(conditionMessage(e), code)) conditionMessage(e) else paste0(code, ": ", conditionMessage(e)), call. = FALSE)
      state <- if (.pc_exit_for(code) == 3L) "NOT_RUN" else if (code %in% c("E_DEQMS_DESIGN_UNSUPPORTED", "E_DEQMS_HYPOTHESIS_UNSUPPORTED", "E_PRODA_LFQ_UNIMPUTED_REQUIRED", "E_PRODA_DESIGN_UNSUPPORTED", "E_PRODA_HYPOTHESIS_UNSUPPORTED")) "INAPPLICABLE" else "FAILED"
      warnings[[length(warnings) + 1L]] <<- .pc_warning(request, code, sprintf("optional model %s: %s", model$model_id, conditionMessage(e)), "assay_engines/model_status.json")
      list(model_id = model$model_id, engine = model$engine, state = state, required = FALSE, reason_code = code, message = conditionMessage(e))
    })
    statuses[[length(statuses) + 1L]] <- status
  }
  rows <- if (length(all_rows)) .pc_rbind_fill(all_rows) else NULL
  family_summaries <- list()
  if (!is.null(rows)) {
    for (family in p$families) { adjusted <- adjust_family(rows, family); rows <- adjusted$rows; family_summaries[[length(family_summaries) + 1L]] <- adjusted$summary }
    write_tsv(rows, "assay_results.tsv", "assay_results", "ProteinZeroNullResult")
    if (length(family_summaries)) write_tsv(do.call(rbind, family_summaries), "families.tsv", "assay_families", "HypothesisFamily")
    write_tsv(do.call(rbind, comparisons), "engine_comparison.tsv", "engine_comparison", "EngineComparison")
  }
  write_json(list(models = statuses, diagnostics = diagnostics, primary_model_id = p$primary_model_id,
                  note = "Alternative engines are separate models with their own families; they never replace the frozen primary limma model and never fall back to limma."),
             "model_status.json", "assay_model_status", "EngineEligibility")
  list(outputs = outputs, warnings = warnings, message = sprintf("assay engines: %d model(s)", length(p$models)))
})

# Row-bind tables whose engine-native columns differ; absent native fields become NA.
.pc_rbind_fill <- function(tables) {
  columns <- unique(unlist(lapply(tables, names)))
  do.call(rbind, lapply(tables, function(t) { for (c in setdiff(columns, names(t))) t[[c]] <- NA; t[columns] }))
}
