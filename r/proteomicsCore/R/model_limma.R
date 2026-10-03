# Observed-data limma adapter and stage handler (packet R05: FR-041 to FR-050).
#
# Fits only genuinely observed cells (no new imputation), with frozen
# trend/robust settings, frozen weights and frozen block correlation.  Every
# contrast is estimated by the exact contrast-as-coefficient refit; the
# moderation universe is identical for every contrast of a model.  Numerical
# failures are explicit typed rows, never P=1.

.pc_verify_plan_inputs <- function(request) {
  plan_path <- .pc_find_input(request, "plan")
  plan <- jsonlite::fromJSON(plan_path, simplifyVector = FALSE)
  if (is.null(request$plan_hash) || !identical(plan$plan_hash, request$plan_hash)) stop("E_PLAN_CHANGED: request plan_hash does not match the frozen plan", call. = FALSE)
  planned <- list(); for (a in plan$artifacts) planned[[a$artifact_id]] <- a
  for (item in request$inputs) {
    if (identical(item$artifact_id, "plan")) next
    reference <- planned[[item$artifact_id]]
    if (is.null(reference)) stop(sprintf("E_PLAN_CHANGED: input %s is not a frozen plan artifact", item$artifact_id), call. = FALSE)
    if (identical(reference$result_type, "DisplayOnlyMatrix")) stop(sprintf("E_DISPLAY_MATRIX_REJECTED: %s is a display-only QC artifact and cannot feed inference", item$artifact_id), call. = FALSE)
    if (!identical(reference$sha256, item$sha256)) stop(sprintf("E_PLAN_CHANGED: input %s hash differs from the frozen plan", item$artifact_id), call. = FALSE)
  }
  plan
}

.pc_read_weights <- function(path, kind, features, observations) {
  table <- .pc_read_tsv(path)
  if (identical(kind, "observation")) {
    w <- stats::setNames(as.numeric(table$weight), table$observation_id)[observations]
    return(as.numeric(w))
  }
  m <- as.matrix(table[, observations, drop = FALSE]); storage.mode(m) <- "numeric"; rownames(m) <- table$feature_id
  m[features, , drop = FALSE]
}

fit_limma_model <- function(model, Y, X, contrasts, estimability, observations, weights, block, correlation, ci_level, context) {
  limma_preflight(model$robust)
  trend <- isTRUE(model$trend); robust <- isTRUE(model$robust)
  est <- estimability[estimability$model_id == model$model_id, , drop = FALSE]
  eligible_any <- unique(est$feature_id[est$eligibility == "eligible"])
  universe <- rownames(Y)[rownames(Y) %in% eligible_any]
  failures <- character()
  if (length(universe)) {
    probe <- .pc_lmfit(Y[universe, , drop = FALSE], X, if (is.matrix(weights)) weights[universe, , drop = FALSE] else weights, block, correlation)
    bad <- universe[!is.finite(probe$sigma) | apply(!is.finite(probe$coefficients), 1L, any)]
    failures <- bad; universe <- setdiff(universe, bad)
  }
  if (!length(universe) && length(failures)) stop(sprintf("E_ENGINE_FAILED: every eligible feature of model %s failed numerically", model$model_id), call. = FALSE)
  if (!length(universe)) stop(sprintf("E_NO_ELIGIBLE_FEATURES: model %s has no eligible features", model$model_id), call. = FALSE)
  YU <- Y[universe, , drop = FALSE]
  WU <- if (is.matrix(weights)) weights[universe, , drop = FALSE] else weights
  full <- .pc_lmfit(YU, X, WU, block, correlation)
  eb_full <- limma::eBayes(full, trend = trend, robust = robust)
  rows <- list(); priors <- list(); exactness <- list()
  hypotheses <- unlist(model$hypotheses)
  for (contrast in contrasts) {
    w <- as.numeric(unlist(contrast$weights))
    ex <- exact_contrast_fit(YU, X, w, WU, block, correlation)
    eb <- limma::eBayes(ex$fit, trend = trend, robust = robust)
    priors[[contrast$contrast_id]] <- eb$s2.prior
    exactness[[contrast$contrast_id]] <- list(path = ex$path, weights = !is.null(weights), block = !is.null(block), correlation = correlation)
    ce <- est[est$contrast_id == contrast$contrast_id, , drop = FALSE]; ce <- ce[match(rownames(Y), ce$feature_id), , drop = FALSE]
    for (hypothesis in hypotheses) {
      out <- endpoint_frame(nrow(Y))
      out$model_id <- model$model_id; out$design_id <- model$design_id; out$contrast_id <- contrast$contrast_id; out$feature_id <- rownames(Y)
      out$run_id <- context$run_id; out$plan_hash <- context$plan_hash; out$engine_version <- as.character(utils::packageVersion("limma"))
      out$role <- if (identical(model$role, "primary") && identical(contrast$role, "primary")) "primary" else "secondary"
      out$estimable <- as.logical(ce$estimable); out$n_obs_by_required_group <- ce$n_obs_by_required_group
      out$hypothesis_type <- if (hypothesis == "treat") "protein_treat" else "protein_zero_null"
      out$result_type <- if (hypothesis == "treat") "ProteinTreatResult" else "ProteinZeroNullResult"
      out$eligibility <- ifelse(ce$eligibility == "eligible", "tested", ce$eligibility)
      out$reason_code <- ifelse(ce$eligibility == "eligible", NA_character_, ce$reason)
      out$eligibility[out$feature_id %in% failures & ce$eligibility == "eligible"] <- "numerical_failure"
      out$reason_code[out$feature_id %in% failures & ce$eligibility == "eligible"] <- "E_NUMERICAL_NONFINITE_FIT"
      tested <- out$eligibility == "tested"
      idx <- match(out$feature_id[tested], universe)
      threshold <- if (is.null(model$effect_threshold)) NA_real_ else as.numeric(model$effect_threshold)
      se <- sqrt(eb$s2.post[idx]) * eb$stdev.unscaled[idx, 1]
      out$effect[tested] <- eb$coefficients[idx, 1]
      out$effect_se[tested] <- se
      out$df_residual[tested] <- eb$df.residual[idx]; out$df_inference[tested] <- eb$df.total[idx]
      ci <- ci_moderated_t(eb$coefficients[idx, 1], se, eb$df.total[idx], ci_level)
      out$ci_lower[tested] <- ci[, "lower"]; out$ci_upper[tested] <- ci[, "upper"]; out$ci_level[tested] <- ci_level
      out$ci_method[tested] <- "effect +/- qt(1-(1-level)/2, df_total) * sqrt(s2_post) * stdev_unscaled"
      out$mean_abundance[tested] <- eb$Amean[idx]
      out$exactness_path[tested] <- ex$path
      if (hypothesis == "treat") {
        if (!is.finite(threshold) || threshold <= 0) stop("E_HYPOTHESIS_UNSUPPORTED: treat requires a positive effect_threshold", call. = FALSE)
        tr <- limma::treat(ex$fit, lfc = threshold, trend = trend, robust = robust)
        out$statistic[tested] <- tr$t[idx, 1]; out$statistic_type[tested] <- "moderated_treat_t"; out$p_value[tested] <- tr$p.value[idx, 1]
        out$effect_threshold[tested] <- threshold; out$null_region[tested] <- sprintf("|log2 effect| <= %s", format(threshold))
      } else {
        out$statistic[tested] <- eb$t[idx, 1]; out$statistic_type[tested] <- "moderated_t"; out$p_value[tested] <- eb$p.value[idx, 1]
        out$null_region[tested] <- "log2 effect = 0"
        out$effect_threshold[tested] <- threshold
      }
      out$display_effect_filter[tested] <- if (is.finite(threshold)) abs(out$effect[tested]) >= threshold else NA
      numeric_fields <- c("effect", "effect_se", "statistic", "p_value", "df_inference")
      broken <- tested & !stats::complete.cases(out[, numeric_fields]) | (tested & !apply(is.finite(as.matrix(out[, numeric_fields])), 1L, all))
      if (any(broken)) { out$eligibility[broken] <- "numerical_failure"; out$reason_code[broken] <- "E_NUMERICAL_NONFINITE_STATISTIC"; out[broken, c(numeric_fields, "ci_lower", "ci_upper")] <- NA }
      native <- out$eligibility == "tested"
      out$native_q_value[native] <- adjust_pvalues(out$p_value[native], "BH")
      rows[[length(rows) + 1L]] <- out
    }
  }
  prior_spread <- if (length(priors) > 1L) max(vapply(priors, function(p) max(abs(p - priors[[1]])), numeric(1))) else 0
  list(rows = do.call(rbind, rows), eb_full = eb_full, universe = universe, failures = failures, prior_spread = prior_spread,
       exactness = exactness, settings = moderation_settings(eb_full, trend, robust), diagnostics = moderation_diagnostics(eb_full, model$model_id, trend, robust))
}

limma_stage <- function(request) .pc_run_stage(request, function(out) {
  p <- request$parameters
  outputs <- list(); warnings <- list()
  emit <- function(relative, id, type) outputs[[length(outputs) + 1L]] <<- .pc_output(out, relative, id, type)
  write_tsv <- function(df, relative, id, type) { .pc_write_tsv(df, file.path(out, relative)); emit(relative, id, type) }
  write_json <- function(value, relative, id, type) { .pc_write_json(value, file.path(out, relative)); emit(relative, id, type) }
  plan <- .pc_verify_plan_inputs(request)
  if (!identical(plan$score_testing, "disabled") && identical(as.integer(plan$requested_phase), 1L)) stop("E_PHASE_CAPABILITY: score testing is not available in Phase 1", call. = FALSE)
  context <- list(run_id = request$run_id, plan_hash = request$plan_hash)
  observed <- .pc_matrix_from_tsv(.pc_find_input(request, "primary_observed_mask"), "logical")
  observations <- .pc_read_tsv(.pc_find_input(request, "primary_observations"))
  estimability <- .pc_read_tsv(.pc_find_input(request, "estimability"))
  estimability$estimable <- estimability$estimable == "true"
  ci_level <- if (is.null(p$ci_level)) 0.95 else p$ci_level
  all_rows <- list(); model_status <- list(); settings <- list(); fits <- list()
  for (model in p$models) {
    status <- tryCatch({
      if (identical(model$role, "primary") && (!isTRUE(model$uses_observed_mask) || !identical(model$matrix_artifact, "primary_matrix")))
        stop("E_PRIMARY_NOT_OBSERVED: the primary model must be fitted on observed cells of the primary matrix", call. = FALSE)   # audit 2026-10-02
      values <- .pc_matrix_from_tsv(.pc_find_input(request, model$matrix_artifact), "numeric")
      if (isTRUE(model$uses_observed_mask)) { mask <- observed[rownames(values), colnames(values), drop = FALSE]; values[!mask] <- NA }
      model_estimability <- estimability
      planned_features <- unique(estimability$feature_id[estimability$model_id == model$model_id])
      if (!isTRUE(model$uses_observed_mask)) {
        present <- intersect(planned_features, rownames(values))
        full <- matrix(NA_real_, length(planned_features), ncol(values), dimnames = list(planned_features, colnames(values)))
        full[present, ] <- values[present, , drop = FALSE]; values <- full
        absent <- model_estimability$model_id == model$model_id & !model_estimability$feature_id %in% present & model_estimability$eligibility == "eligible"
        model_estimability$eligibility[absent] <- "excluded"; model_estimability$reason[absent] <- "absent_from_sensitivity_universe"
      }
      X <- .pc_matrix_from_tsv(.pc_find_input(request, paste0("design_", model$design_id)), "numeric")
      if (!identical(rownames(X), colnames(values))) stop("E_ID_ALIGNMENT: design rows do not match matrix observations", call. = FALSE)
      weights <- if (is.null(model$weights_artifact)) NULL else .pc_read_weights(.pc_find_input(request, model$weights_artifact), model$weights_kind, rownames(values), colnames(values))
      block <- NULL; correlation <- NULL
      if (identical(model$blocking$mode, "duplicate_correlation")) {
        block <- observations[[model$blocking$subject_column]][match(colnames(values), observations$observation_id)]
        correlation <- as.numeric(model$blocking$consensus_correlation)
      }
      contrasts <- Filter(function(c) identical(c$design_id, model$design_id), p$contrasts)
      fit <- fit_limma_model(model, values, X, contrasts, model_estimability, observations, weights, block, correlation, ci_level, context)
      fits[[model$model_id]] <- list(fit = fit, values = values, X = X, weights = weights, block = block, correlation = correlation, contrasts = contrasts)
      all_rows[[length(all_rows) + 1L]] <- fit$rows
      settings[[model$model_id]] <- c(fit$settings, list(fitting_universe_size = length(fit$universe), numerical_failures = I(fit$failures),
                                                         prior_identical_across_contrasts = fit$prior_spread <= 1e-12, exactness = fit$exactness,
                                                         input_matrix = model$matrix_artifact, observed_cells_only = isTRUE(model$uses_observed_mask) && identical(model$matrix_artifact, "primary_matrix"),
                                                         input_imputation = if (is.null(model$input_imputation)) "none" else model$input_imputation,
                                                         new_primary_imputation = if (identical(model$role, "primary")) "none" else "not_applicable_non_primary_model",
                                                         weights = model$weights_kind, block_correlation = correlation))
      write_tsv(fit$diagnostics, file.path("diagnostics", paste0(model$model_id, "_moderation.tsv")), paste0("moderation_", model$model_id), "ModelFit")
      list(model_id = model$model_id, state = "COMPLETED", required = identical(model$execution_requirement, "required"), reason_code = NULL)
    }, error = function(e) {
      if (identical(model$execution_requirement, "required")) stop(e)
      code <- regmatches(conditionMessage(e), regexpr("^E_[A-Z0-9_]+", conditionMessage(e))); if (!length(code)) code <- "E_ENGINE_FAILED"
      warnings[[length(warnings) + 1L]] <<- .pc_warning(request, code, sprintf("optional model %s failed: %s", model$model_id, conditionMessage(e)), "dea/model_status.json")
      list(model_id = model$model_id, state = if (.pc_exit_for(code) == 3L) "NOT_RUN" else "FAILED", required = FALSE, reason_code = code, message = conditionMessage(e))
    })
    model_status[[length(model_status) + 1L]] <- status
  }
  rows <- if (length(all_rows)) do.call(rbind, all_rows) else endpoint_frame(0L)

  # omnibus endpoints (exact per-feature GLS moderated F)
  omnibus_rows <- list()
  for (test in p$omnibus_tests) {
    f <- fits[[test$model_id]]; if (is.null(f)) next
    index <- match(unlist(test$coefficient_names), colnames(f$X))
    if (any(is.na(index))) stop(sprintf("E_OMNIBUS_COEFFICIENT: omnibus %s names unknown coefficients", test$id), call. = FALSE)
    U <- f$fit$universe
    Fstat <- omnibus_moderated_f(f$values[U, , drop = FALSE], f$X, index, f$fit$eb_full, if (is.matrix(f$weights)) f$weights[U, , drop = FALSE] else f$weights, f$block, f$correlation)
    o <- endpoint_frame(nrow(f$values)); o$feature_id <- rownames(f$values); o$model_id <- test$model_id; o$design_id <- unique(vapply(f$contrasts, function(c) c$design_id, ""))[1]
    o$contrast_id <- test$id; o$result_type <- "ProteinOmnibusResult"; o$hypothesis_type <- "protein_omnibus"; o$run_id <- context$run_id; o$plan_hash <- context$plan_hash
    o$engine_version <- as.character(utils::packageVersion("limma")); o$role <- "secondary"; o$effect_scale <- NA_character_
    o$eligibility <- ifelse(o$feature_id %in% U, "tested", ifelse(o$feature_id %in% f$fit$failures, "numerical_failure", "excluded"))
    o$reason_code <- ifelse(o$eligibility == "excluded", "not_in_model_fitting_universe", ifelse(o$eligibility == "numerical_failure", "E_NUMERICAL_NONFINITE_FIT", NA_character_))
    o$estimable <- o$eligibility == "tested"
    m <- match(o$feature_id, Fstat$feature_id)
    o$statistic <- Fstat$statistic[m]; o$statistic_type <- ifelse(o$eligibility == "tested", sprintf("moderated_F_df1_%d", length(index)), NA_character_)
    o$df_inference <- Fstat$df2[m]; o$p_value <- Fstat$p_value[m]
    o$null_region <- ifelse(o$eligibility == "tested", paste0(paste(unlist(test$coefficient_names), collapse = " = "), " = 0"), NA_character_)
    tested <- o$eligibility == "tested"; o$native_q_value[tested] <- adjust_pvalues(o$p_value[tested], "BH")
    omnibus_rows[[length(omnibus_rows) + 1L]] <- o
  }

  # declared multiplicity families (pooled once per family)
  family_summaries <- list()
  combined <- rbind(rows, if (length(omnibus_rows)) do.call(rbind, omnibus_rows) else NULL)
  for (family in p$families) {
    adjusted <- adjust_family(combined, family)
    combined <- adjusted$rows; family_summaries[[length(family_summaries) + 1L]] <- adjusted$summary
  }
  planned <- character()
  for (model in p$models) {
    f <- fits[[model$model_id]]; if (is.null(f)) next
    htypes <- ifelse(unlist(model$hypotheses) == "treat", "protein_treat", "protein_zero_null")
    for (contrast in f$contrasts) for (h in htypes) planned <- c(planned, paste(model$model_id, contrast$contrast_id, rownames(f$values), h, sep = "|"))
  }
  for (test in p$omnibus_tests) if (!is.null(fits[[test$model_id]])) planned <- c(planned, paste(test$model_id, test$id, rownames(fits[[test$model_id]]$values), "protein_omnibus", sep = "|"))
  verify_differential_table(combined, planned)
  zero <- combined[combined$hypothesis_type == "protein_zero_null", , drop = FALSE]
  write_tsv(zero, "zero_null.tsv", "dea_zero_null", "ProteinZeroNullResult")
  if (any(combined$hypothesis_type == "protein_treat")) write_tsv(combined[combined$hypothesis_type == "protein_treat", , drop = FALSE], "treat.tsv", "dea_treat", "ProteinTreatResult")
  if (any(combined$hypothesis_type == "protein_omnibus")) write_tsv(combined[combined$hypothesis_type == "protein_omnibus", , drop = FALSE], "omnibus.tsv", "dea_omnibus", "ProteinOmnibusResult")
  families <- if (length(family_summaries)) do.call(rbind, family_summaries) else data.frame(family_id = character())
  write_tsv(families, "families.tsv", "dea_families", "HypothesisFamily")

  # whole-subject influence for the primary model (descriptive)
  if (isTRUE(p$influence$enabled)) {
    for (model in p$models) {
      f <- fits[[model$model_id]]
      if (is.null(f) || !identical(model$role, "primary")) next
      unit_column <- if (!is.null(model$blocking$subject_column)) model$blocking$subject_column else "subject_id"
      ids <- observations[match(colnames(f$values), observations$observation_id), , drop = FALSE]
      subjects <- ifelse(is.na(ids[[unit_column]]), ids$biological_unit_id, ids[[unit_column]])
      groups <- ids[[p$group_column]]
      U <- f$fit$universe
      references <- lapply(f$contrasts, function(c) exact_contrast_fit(f$values[U, , drop = FALSE], f$X, as.numeric(unlist(c$weights)), if (is.matrix(f$weights)) f$weights[U, , drop = FALSE] else f$weights, f$block, f$correlation)$fit$coefficients[, 1])
      names(references) <- vapply(f$contrasts, function(c) c$contrast_id, "")
      inf <- subject_influence(f$values[U, , drop = FALSE], f$X, f$contrasts, subjects, groups, if (is.matrix(f$weights)) f$weights[U, , drop = FALSE] else f$weights, f$block, f$correlation, references)
      inf$model_id <- model$model_id
      write_tsv(inf, file.path("influence", paste0(model$model_id, "_subject_influence.tsv")), paste0("influence_", model$model_id), "InfluenceResult")
    }
  }
  write_json(list(models = model_status, settings = settings, score_p_values_emitted = FALSE,
                  note = "Zero discoveries with all required work executed is a completed analysis, not a failure and not evidence of no biological effect."),
             "model_status.json", "dea_model_status", "ModelFit")
  n_tested <- sum(zero$eligibility == "tested")
  list(outputs = outputs, warnings = warnings, message = sprintf("limma fitted %d model(s); %d zero-null endpoints tested", sum(vapply(model_status, function(s) identical(s$state, "COMPLETED"), logical(1))), n_tested))
})
