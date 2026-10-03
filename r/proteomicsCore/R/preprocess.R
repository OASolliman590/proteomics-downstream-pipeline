# Preprocessing stage handler (packet R03: FR-021 to FR-030).
#
# Reads only a verified R02 canonical bundle.  Applies prespecified
# exclusions and normalization to produce a distinct primary analysis
# matrix artifact (preserve keeps values byte-for-byte equal), then derives
# coverage masks, missingness diagnostics, QC tables, optional imputation
# sensitivities and the optional detection-only endpoint.  No model is fit.

# D-42: every text artifact is written as UTF-8 with LF line endings through a binary connection, so its bytes (and the
# plan hash built from them) do not depend on the platform (a text-mode connection writes CRLF on Windows).
.pc_write_lf <- function(text, path) {
  connection <- file(path, open = "wb")
  on.exit(close(connection), add = TRUE)
  writeBin(charToRaw(paste0(enc2utf8(paste(text, collapse = "\n")), "\n")), connection)
  invisible(path)
}

.pc_write_json <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  .pc_write_lf(jsonlite::toJSON(value, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = NA), path)
}

.pc_matrix_df <- function(values, formatter = .pc_fmt) {
  df <- data.frame(feature_id = rownames(values), stringsAsFactors = FALSE, check.names = FALSE)
  for (j in seq_len(ncol(values))) df[[colnames(values)[j]]] <- formatter(values[, j])
  df
}

.pc_resolve_exclusions <- function(exclusions, ids) {
  rows <- data.frame(observation_id = character(), reason = character(), stringsAsFactors = FALSE)
  for (item in exclusions) {
    reason <- item$reason
    if (is.null(reason) || !nzchar(trimws(reason))) stop(sprintf("E_EXCLUSION_REASON_REQUIRED: exclusion of %s needs a nonblank reason", item$observation_id), call. = FALSE)
    if (!item$observation_id %in% ids) stop(sprintf("E_EXCLUSION_UNKNOWN_OBSERVATION: %s is not a canonical observation", item$observation_id), call. = FALSE)
    if (item$observation_id %in% rows$observation_id) stop(sprintf("E_ID_DUPLICATE: %s is excluded twice", item$observation_id), call. = FALSE)
    rows[nrow(rows) + 1L, ] <- list(item$observation_id, reason)
  }
  rows[order(rows$observation_id), , drop = FALSE]
}

preprocess_stage <- function(request) .pc_run_stage(request, function(out) {
  p <- request$parameters
  outputs <- list(); warnings <- list()
  emit <- function(relative, id, type) outputs[[length(outputs) + 1L]] <<- .pc_output(out, relative, id, type)
  write_tsv <- function(df, relative, id, type) { .pc_write_tsv(df, file.path(out, relative)); emit(relative, id, type) }
  write_json <- function(value, relative, id, type) { .pc_write_json(value, file.path(out, relative)); emit(relative, id, type) }

  manifest_path <- .pc_find_input(request, "canonical_manifest")
  bundle <- read_canonical_bundle(dirname(manifest_path))
  matrix_artifact <- Filter(function(a) identical(a$artifact_id, "matrix"), bundle$manifest$artifacts)[[1L]]
  source_hash <- matrix_artifact$sha256
  group_column <- p$group_column
  values <- bundle$values; observed <- bundle$observed; prior <- bundle$prior_imputed; obs <- bundle$observations
  scope <- if (is.null(p$scope)) "analysis" else p$scope
  output_scale <- bundle$scale$output_scale

  # 1. prespecified exclusions (flags never remove samples; only declared exclusions do)
  exclusions <- .pc_resolve_exclusions(p$exclusions, colnames(values))
  # 2. TMT channel roles
  assay <- p$assay; tmt <- identical(assay, "tmt"); reference_only <- character(); tmt_eligibility <- NULL
  normalization <- p$normalization
  if (identical(normalization, "tmt_loading") && !tmt) stop("E_NORMALIZATION_ASSAY: tmt_loading normalization requires assay=tmt", call. = FALSE)
  if (normalization != "preserve" && !output_scale %in% c("log2")) stop(sprintf("E_NORMALIZATION_SCALE: %s normalization requires log2 abundance, input scale is %s", normalization, output_scale), call. = FALSE)
  keep <- !colnames(values) %in% exclusions$observation_id
  values <- values[, keep, drop = FALSE]; obs <- obs[keep, , drop = FALSE]
  if (!is.null(observed)) { observed <- observed[, keep, drop = FALSE]; prior <- prior[, keep, drop = FALSE] }

  factors <- NULL
  if (tmt) {
    .pc_tmt_columns(obs)
    strategy <- p$tmt_strategy
    if (is.null(strategy)) stop("E_TMT_STRATEGY_REQUIRED: TMT input must declare tmt_strategy", call. = FALSE)
    if (identical(strategy, "upstream_preserved") && normalization != "preserve") stop("E_TMT_STRATEGY_CONFLICT: upstream_preserved requires normalization=preserve", call. = FALSE)
    if (strategy %in% c("bridge", "no_bridge") && normalization != "tmt_loading") stop("E_TMT_STRATEGY_CONFLICT: bridge/no_bridge strategies require normalization=tmt_loading", call. = FALSE)
    if (normalization == "tmt_loading") {
      loading <- tmt_loading(values, obs, observed); values <- loading$values; factors <- loading$factors
    }
    bridge_flag <- obs$channel_role == "bridge"
    if (!is.null(p$bridge_channel_column)) {
      column <- p$bridge_channel_column
      if (!column %in% names(obs)) stop(sprintf("E_TMT_METADATA: bridge_channel_column %s is absent", column), call. = FALSE)
      bridge_flag <- obs[[column]] %in% "true"
    }
    if (identical(strategy, "bridge")) {
      bridged <- tmt_bridge(values, obs, bridge_flag)
      values <- bridged$values; tmt_eligibility <- bridged$eligibility
      write_tsv(bridged$bridge, "tmt/bridge_medians.tsv", "tmt_bridge_medians", "PreprocessingResult")
      write_tsv(tmt_eligibility, "tmt/bridge_eligibility.tsv", "tmt_bridge_eligibility", "CoverageResult")
    }
    if (identical(strategy, "no_bridge")) tmt_no_bridge_check(obs, group_column)
    reference_only <- colnames(values)[obs$channel_role != "sample" | bridge_flag]
    analytic <- !colnames(values) %in% reference_only
    write_tsv(data.frame(observation_id = reference_only, role = "reference_only_not_a_biological_unit", stringsAsFactors = FALSE), "tmt/reference_channels.tsv", "tmt_reference_channels", "ObservationHierarchy")
    values <- values[, analytic, drop = FALSE]; obs <- obs[analytic, , drop = FALSE]
    if (!is.null(observed)) { observed <- observed[, analytic, drop = FALSE]; prior <- prior[, analytic, drop = FALSE] }
  } else {
    if (identical(normalization, "median")) { n <- normalize_median(values, observed); values <- n$values; factors <- n$factors }
    else if (identical(normalization, "reference")) { n <- normalize_reference(values, observed, p$reference_features, if (is.null(p$reference_minimum_fraction)) 1 else p$reference_minimum_fraction); values <- n$values; factors <- n$factors }
    else factors <- normalize_preserve(values)$factors
  }
  if (is.null(factors)) factors <- normalize_preserve(values)$factors
  if (!ncol(values)) stop("E_NO_OBSERVATIONS: no analytical observations remain after exclusions", call. = FALSE)

  # 3. primary analysis matrix artifact and lineage
  write_tsv(.pc_matrix_df(values), "primary/matrix.tsv", "primary_matrix", "MatrixArtifact")
  write_tsv(.pc_matrix_df(!is.na(values), .pc_fmt_bool), "primary/numeric_mask.tsv", "primary_numeric_mask", "NumericAvailabilityMask")
  if (!is.null(observed)) {
    write_tsv(.pc_matrix_df(observed, .pc_fmt_bool), "primary/observed_mask.tsv", "primary_observed_mask", "OriginalObservedMask")
    write_tsv(.pc_matrix_df(prior, .pc_fmt_bool), "primary/prior_imputed_mask.tsv", "primary_prior_imputed_mask", "PriorImputedMask")
  }
  write_tsv(obs, "primary/observations.tsv", "primary_observations", "ObservationMetadata")
  write_tsv(bundle$features, "primary/features.tsv", "primary_features", "FeatureMetadata")
  write_tsv(factors, "normalization_factors.tsv", "normalization_factors", "PreprocessingResult")
  values_changed <- normalization != "preserve" || tmt
  lineage <- list(source_artifact_id = "canonical.matrix", source_sha256 = source_hash, output_artifact_id = "primary_matrix",
                  output_sha256 = sha256_file(file.path(out, "primary/matrix.tsv")), input_scale = output_scale, output_scale = output_scale,
                  operations = list(list(step = "exclusion", excluded = exclusions$observation_id), list(step = "normalization", method = normalization, values_changed = values_changed),
                                    list(step = "reference_only_channels", removed = reference_only)),
                  new_primary_imputation = "none", display_matrix_used = FALSE)
  write_json(lineage, "primary/lineage.json", "primary_lineage", "TransformLineage")
  write_tsv(exclusions, "exclusions/exclusions.tsv", "exclusions", "ExclusionDecision")

  # 4. coverage masks for every planned model x contrast
  coverage_state <- list(state = "NOT_REQUESTED", reason_code = NULL)
  if (scope == "analysis" && length(p$coverage_rules)) {
    extra <- NULL
    if (!is.null(tmt_eligibility)) { bad <- tmt_eligibility$feature_id[!tmt_eligibility$bridge_eligible]; extra <- stats::setNames(as.list(rep("tmt_bridge_measurement_missing_in_plex", length(bad))), bad) }
    if (is.null(observed)) {
      coverage_state <- list(state = "NOT_RUN", reason_code = "E_ORIGINAL_MASK_REQUIRED", message = "original-observed mask unknown; primary observed-coverage inference is blocked")
    } else {
      cov <- coverage_tables(observed, obs, group_column, p$coverage_rules, extra)
      write_tsv(cov$summary, "coverage/coverage.tsv", "coverage", "CoverageResult")
      write_tsv(cov$by_group, "coverage/coverage_by_group.tsv", "coverage_by_group", "CoverageResult")
      coverage_state <- list(state = "COMPLETED", reason_code = NULL, n_eligible = sum(cov$summary$eligibility == "eligible"), n_rows = nrow(cov$summary))
    }
  }

  # 5. missingness diagnostics
  miss <- missingness_tables(values, observed, prior, obs, group_column)
  write_tsv(miss$by_observation, "missingness/by_observation.tsv", "missingness_by_observation", "MissingnessSummary")
  write_tsv(miss$by_group, "missingness/by_group.tsv", "missingness_by_group", "MissingnessSummary")
  write_tsv(miss$by_feature, "missingness/by_feature.tsv", "missingness_by_feature", "MissingnessSummary")
  if (nrow(miss$abundance)) write_tsv(miss$abundance, "missingness/abundance_dependence.tsv", "missingness_abundance", "MissingnessSummary")
  write_json(miss$statement, "missingness/statement.json", "missingness_statement", "MissingnessSummary")

  # 6. QC diagnostics (display-only PCA input is a separate artifact)
  n_injections <- if ("n_injections" %in% names(obs)) as.integer(obs$n_injections) else rep(1L, nrow(obs))
  sample_n <- do.call(rbind, lapply(unique(obs[[group_column]]), function(g) { sel <- obs[[group_column]] == g
    data.frame(group = g, n_biological_units = sum(sel), n_injections = sum(n_injections[sel]), n_subjects = length(unique(ifelse(is.na(obs$subject_id[sel]), obs$observation_id[sel], obs$subject_id[sel]))), stringsAsFactors = FALSE) }))
  write_tsv(sample_n, "qc/sample_n.tsv", "qc_sample_n", "QCResult")
  write_tsv(qc_distributions(values), "qc/distributions.tsv", "qc_distributions", "QCResult")
  usable <- if (is.null(observed)) !is.na(values) else observed & !is.na(values)
  correlations <- qc_correlations(values, usable)
  write_tsv(correlations, "qc/correlations.tsv", "qc_correlations", "QCResult")
  pca <- qc_pca(values)
  write_tsv(.pc_matrix_df(pca$display$centered), "qc/pca_display_input.tsv", "qc_pca_display_input", "DisplayOnlyMatrix")
  pca_state <- list(state = pca$state, reason_code = pca$reason_code, message = pca$message, display_rule = "features with >=1 value; NA filled with feature observed median; centered by feature; display only, never an inference input",
                    dropped_all_missing_features = pca$display$dropped_all_missing, primary_matrix_sha256 = sha256_file(file.path(out, "primary/matrix.tsv")))
  if (identical(pca$state, "COMPLETED")) {
    scores <- data.frame(observation_id = rownames(pca$scores), pca$scores, check.names = FALSE, stringsAsFactors = FALSE)
    write_tsv(scores, "qc/pca_scores.tsv", "qc_pca_scores", "QCResult")
    write_tsv(pca$variance, "qc/pca_variance.tsv", "qc_pca_variance", "QCResult")
  }
  write_json(pca_state, "qc/pca_state.json", "qc_pca_state", "QCResult")
  watch <- qc_watchlist(pca, correlations, colnames(values), if (is.null(p$watchlist_threshold)) 3.5 else p$watchlist_threshold)
  write_tsv(watch, "qc/watchlist.tsv", "qc_watchlist", "ExclusionDecision")
  for (id in watch$observation_id[watch$flagged]) warnings[[length(warnings) + 1L]] <- .pc_warning(request, "W_QC_WATCHLIST", sprintf("observation %s flagged by QC watchlist (not removed)", id), "qc/watchlist.tsv")
  fragment <- list(exclusions = lapply(seq_len(nrow(exclusions)), function(i) list(observation_id = exclusions$observation_id[i], reason = exclusions$reason[i])),
                   retained_observation_ids = I(colnames(values)), watchlist_policy = list(method = "pca_robust_z", threshold = if (is.null(p$watchlist_threshold)) 3.5 else p$watchlist_threshold, action = "flag_only"),
                   watchlist_flagged = I(watch$observation_id[watch$flagged]))
  write_json(fragment, "exclusions/fragment.json", "exclusion_fragment", "ExclusionDecision")

  # 7. imputation and normalization sensitivities (primary values never change)
  primary_before <- sha256_file(file.path(out, "primary/matrix.tsv"))
  sensitivity_records <- list()
  for (spec in p$sensitivities) {
    result <- run_sensitivity(values, spec, request$rng$seed)
    base <- file.path("sensitivity", spec$id)
    write_tsv(.pc_matrix_df(result$values), file.path(base, "matrix.tsv"), paste0("sensitivity_", spec$id, "_matrix"), "SensitivityMatrix")
    write_tsv(.pc_matrix_df(result$imputed, .pc_fmt_bool), file.path(base, "imputed_mask.tsv"), paste0("sensitivity_", spec$id, "_imputed_mask"), "SensitivityMask")
    write_json(c(result$parameters, list(id = spec$id, model_id = spec$model_id, excluded_features = I(result$excluded), role = "sensitivity_only_secondary_model_input",
                                         limitation = "imputed values are not additional independent information")),
               file.path(base, "parameters.json"), paste0("sensitivity_", spec$id, "_parameters"), "SensitivityResult")
    sensitivity_records[[spec$id]] <- list(method = spec$method, model_id = spec$model_id, n_features = nrow(result$values))
  }
  if ("quantile" %in% unlist(p$normalization_sensitivities)) {
    q <- normalize_quantile_sensitivity(values)
    write_tsv(.pc_matrix_df(q$values), "sensitivity/quantile/matrix.tsv", "sensitivity_quantile_matrix", "SensitivityMatrix")
    write_json(list(method = "quantile", function_called = "limma::normalizeBetweenArrays(method='quantile')", universe = "complete_case", lost_features = I(q$lost_features), role = "named_normalization_sensitivity_only"),
               "sensitivity/quantile/parameters.json", "sensitivity_quantile_parameters", "SensitivityResult")
  }
  if (!identical(sha256_file(file.path(out, "primary/matrix.tsv")), primary_before)) stop("E_INTEGRITY: primary matrix changed during sensitivity generation", call. = FALSE)

  # 8. detection-only exploratory endpoint
  detection_state <- list(state = "NOT_REQUESTED")
  if (isTRUE(p$detection$enabled)) {
    detection <- tryCatch(detection_tests(observed, obs, group_column, p$detection$comparisons, p$blocking_mode, if (is.null(p$detection$family_id)) "detection-exploratory" else p$detection$family_id),
                          error = function(e) e)
    if (inherits(detection, "error")) {
      code <- regmatches(conditionMessage(detection), regexpr("^E_[A-Z0-9_]+", conditionMessage(detection)))
      if (!length(code) || !code %in% c("E_DETECTION_DESIGN_UNSUPPORTED", "E_ORIGINAL_MASK_REQUIRED")) stop(detection)
      detection_state <- list(state = "INAPPLICABLE", reason_code = code, message = conditionMessage(detection), p_values_emitted = FALSE)
      warnings[[length(warnings) + 1L]] <- .pc_warning(request, code, conditionMessage(detection), "detection/status.json")
    } else {
      write_tsv(detection, "detection/detection.tsv", "detection", "DetectionResult")
      detection_state <- list(state = "COMPLETED", family_id = unique(detection$family_id), n_tests = nrow(detection), p_values_emitted = TRUE,
                              note = "exploratory detection family; never merged with protein abundance families")
    }
    write_json(detection_state, "detection/status.json", "detection_status", "DetectionResult")
  }

  summary <- list(result_type = "PreprocessingResult", scope = scope, assay = assay, scale = output_scale, normalization = normalization,
                  source_matrix_sha256 = source_hash, primary_matrix_sha256 = primary_before, primary_values_changed = values_changed,
                  original_mask_state = bundle$mask_state, n_features = nrow(values), n_observations = ncol(values),
                  excluded_observations = I(exclusions$observation_id), reference_only_channels = I(reference_only),
                  coverage = coverage_state, pca = list(state = pca$state, reason_code = pca$reason_code), detection = detection_state,
                  sensitivities = sensitivity_records, model_fit_performed = FALSE)
  write_json(summary, "preprocessing_result.json", "preprocessing_result", "PreprocessingResult")
  list(outputs = outputs, warnings = warnings, message = sprintf("preprocessing completed: %d features x %d observations (%s)", nrow(values), ncol(values), normalization))
})
