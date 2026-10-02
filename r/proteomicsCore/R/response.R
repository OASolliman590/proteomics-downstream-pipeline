# Treatment-response stage (packet R09: FR-081 to FR-090, SM19-SM24).
#
# Descriptive same-model axes d = U - C, t = T - U, r = T - C with full joint
# covariance (no P/q fields); optional TOST equivalence on r; optional formal
# rescue conjunction with independently frozen directions; fixed independent
# scores with exact or Monte Carlo randomization only when every training
# and selection participant is disjoint from the tested units.

.rs_family_adjust <- function(p, adjustment) {
  out <- rep(NA_real_, length(p)); ok <- is.finite(p)
  if (any(ok)) out[ok] <- adjust_pvalues(p[ok], if (identical(adjustment, "BY")) "BY" else "BH")
  out
}

.rs_envelope <- function(n, context, result_type, hypothesis_type, family_id, model_id) {
  data.frame(schema_version = rep("1.2.0", n), run_id = context$run_id, plan_hash = context$plan_hash, result_type = result_type, hypothesis_type = hypothesis_type,
             family_id = family_id, engine = "limma", engine_version = as.character(utils::packageVersion("limma")), model_id = model_id, stringsAsFactors = FALSE)
}

response_stage <- function(request) .pc_run_stage(request, function(out) {
  p <- request$parameters
  outputs <- list(); warnings <- list()
  emit <- function(relative, id, type) outputs[[length(outputs) + 1L]] <<- .pc_output(out, relative, id, type)
  write_tsv <- function(df, relative, id, type) { .pc_write_tsv(df, file.path(out, relative)); emit(relative, id, type) }
  write_json <- function(value, relative, id, type) { .pc_write_json(value, file.path(out, relative)); emit(relative, id, type) }
  plan <- .pm_verify_inputs_response(request)
  context <- list(run_id = request$run_id, plan_hash = request$plan_hash)
  values <- .pc_matrix_from_tsv(.pc_find_input(request, "primary_matrix"), "numeric")
  observed <- .pc_matrix_from_tsv(.pc_find_input(request, "primary_observed_mask"), "logical")
  values[!observed[rownames(values), colnames(values)]] <- NA
  obs <- .pc_read_tsv(.pc_find_input(request, "primary_observations")); obs <- obs[match(colnames(values), obs$observation_id), , drop = FALSE]
  dea <- .pc_read_tsv(.pc_find_input(request, "dea_zero_null"))
  summary <- list(axes = list(), equivalence = NULL, formal_rescue = NULL, scores = list())

  descriptive_rows <- list(); cov_rows <- list(); eq_rows <- list(); rescue_rows <- list()
  for (axis in p$axes) {
    X <- .pc_matrix_from_tsv(.pc_find_input(request, paste0("design_", axis$design_id)), "numeric")
    moderation <- .pc_read_tsv(.pc_find_input(request, paste0("moderation_", axis$model_id)))
    s2post <- stats::setNames(as.numeric(moderation$s2_posterior), moderation$feature_id)
    wd <- as.numeric(unlist(axis$disease_weights)); wt <- as.numeric(unlist(axis$treatment_weights)); wr <- as.numeric(unlist(axis$residual_weights))
    if (max(abs(wr - (wd + wt))) > 1e-12) stop(sprintf("E_AXIS_INCOHERENT: axis %s residual weights differ from disease + treatment", axis$id), call. = FALSE)
    C <- cbind(d = wd, t = wt, r = wr)
    weights <- if (is.null(axis$weights_artifact)) NULL else .pc_read_weights(.pc_find_input(request, axis$weights_artifact), "observation", rownames(values), colnames(values))
    block <- if (identical(axis$blocking_mode, "duplicate_correlation")) obs[[axis$subject_column]] else NULL
    rho <- if (is.null(block)) NULL else as.numeric(axis$consensus_correlation)
    rows_of <- function(cid) { r <- dea[dea$model_id == axis$model_id & dea$contrast_id == cid, , drop = FALSE]; r[match(rownames(values), r$feature_id), , drop = FALSE] }
    D <- rows_of(axis$disease_contrast); Tt <- rows_of(axis$treatment_contrast); R <- rows_of(axis$residual_contrast)
    tested <- D$eligibility == "tested" & Tt$eligibility == "tested" & R$eligibility == "tested"
    tested[is.na(tested)] <- FALSE
    n <- nrow(values)
    cv <- matrix(NA_real_, n, 6, dimnames = list(rownames(values), c("u_dd", "u_tt", "u_rr", "u_dt", "u_dr", "u_tr")))
    for (i in which(tested)) {
      U <- feature_contrast_covariance(values[i, ], X, C, weights, block, rho)
      if (!is.null(U)) cv[i, ] <- c(U["d", "d"], U["t", "t"], U["r", "r"], U["d", "t"], U["d", "r"], U["t", "r"])
    }
    s2 <- s2post[rownames(values)]
    d <- as.numeric(D$effect); t <- as.numeric(Tt$effect); r <- as.numeric(R$effect)
    algebra <- abs(r - (d + t))
    if (any(algebra[tested] > 1e-8, na.rm = TRUE)) stop(sprintf("E_AXIS_INCOHERENT: axis %s exported effects violate r = d + t", axis$id), call. = FALSE)
    ri <- reversal_index(d, t, p$dmin); ri[!tested] <- NA
    classes <- vapply(ri, response_class, "")
    shared <- length(intersect(unlist(axis$disease_groups), unlist(axis$treatment_groups))) > 0L
    df_total <- as.numeric(R$df_inference)
    out_rows <- cbind(.rs_envelope(n, context, "DescriptiveResponse", "descriptive_reversal", NA_character_, axis$model_id),
      data.frame(axis_id = axis$id, disease_contrast_id = axis$disease_contrast, treatment_contrast_id = axis$treatment_contrast, residual_contrast_id = axis$residual_contrast,
                 contrast_id = NA_character_, feature_id = rownames(values), estimable = tested, eligibility = ifelse(tested, "tested", "excluded"),
                 reason_code = ifelse(tested, NA_character_, "axis_contrast_not_tested"), d = ifelse(tested, d, NA), t = ifelse(tested, t, NA), r = ifelse(tested, r, NA),
                 var_d = cv[, "u_dd"] * s2, var_t = cv[, "u_tt"] * s2, var_r = cv[, "u_rr"] * s2, cov_dt = cv[, "u_dt"] * s2,
                 unscaled_var_d = cv[, "u_dd"], unscaled_var_t = cv[, "u_tt"], unscaled_var_r = cv[, "u_rr"], unscaled_cov_dt = cv[, "u_dt"],
                 algebra_abs_error = ifelse(tested, algebra, NA), dmin = p$dmin, reversal_index = ri, descriptive_class = classes,
                 directional_opposition = ifelse(tested & d != 0 & t != 0, sign(t) == -sign(d), ifelse(tested, FALSE, NA)),
                 residual = ifelse(is.na(ri), NA, d * (1 - ri)), crossed_control = ifelse(is.na(ri), NA, ri > 1), residual_exceeds_disease = ifelse(is.na(ri), NA, ri > 2),
                 shared_control_warning = shared, stringsAsFactors = FALSE))
    if (identical(p$ratio_uncertainty, "fieller")) {
      q <- stats::qt(1 - p$alpha / 2, df_total)
      kinds <- character(n); intervals <- character(n)
      for (i in seq_len(n)) {
        if (!tested[i] || is.na(ri[i]) || !is.finite(cv[i, "u_dd"])) { kinds[i] <- "unavailable"; intervals[i] <- "[]"; next }
        f <- fieller_interval(-t[i], d[i], cv[i, "u_tt"] * s2[i], cv[i, "u_dd"] * s2[i], -cv[i, "u_dt"] * s2[i], q[i])
        kinds[i] <- f$kind; intervals[i] <- as.character(jsonlite::toJSON(lapply(f$intervals, function(v) ifelse(is.finite(v), v, ifelse(v > 0, "Inf", "-Inf"))), auto_unbox = FALSE, digits = NA))
      }
      out_rows$ratio_interval_kind <- kinds; out_rows$ratio_interval <- intervals; out_rows$ratio_interval_level <- 1 - p$alpha
    } else {
      out_rows$ratio_interval_kind <- "not_requested"
    }
    descriptive_rows[[length(descriptive_rows) + 1L]] <- out_rows
    if (shared) warnings[[length(warnings) + 1L]] <- .pc_warning(request, "W_SHARED_CONTROL", sprintf("axis %s: disease and treatment contrasts share a group; their estimates are negatively correlated and opposition is descriptive, not independent confirmation", axis$id), "response/descriptive.tsv")
    if (p$response_mode %in% c("equivalence", "formal_rescue")) {
      se_r <- as.numeric(R$effect_se)
      eq <- tost(r, se_r, df_total, p$equivalence_margin, p$alpha)
      eq_df <- cbind(.rs_envelope(n, context, "EquivalenceResult", "equivalence", axis$equivalence_family, axis$model_id),
                     data.frame(axis_id = axis$id, contrast_id = axis$residual_contrast, feature_id = rownames(values), eligibility = ifelse(tested, "tested", "excluded"),
                                reason_code = ifelse(tested, NA_character_, "axis_contrast_not_tested"), residual_estimate = r, residual_se = se_r, df = df_total,
                                margin = p$equivalence_margin, alpha = p$alpha, eq[, c("p_lower", "p_upper", "p_value", "ci_lower", "ci_upper", "ci_level")], stringsAsFactors = FALSE))
      eq_df[!tested, c("p_lower", "p_upper", "p_value", "ci_lower", "ci_upper")] <- NA
      eq_df$q_value <- .rs_family_adjust(eq_df$p_value, axis$equivalence_adjustment)
      eq_rows[[length(eq_rows) + 1L]] <- eq_df
    }
    if (identical(p$response_mode, "formal_rescue")) {
      direction <- jsonlite::fromJSON(.pc_find_input(request, "direction_resource"), simplifyVector = FALSE)
      tested_subjects <- unique(ifelse(is.na(obs$subject_id), obs$biological_unit_id, obs$subject_id))
      overlap <- intersect(unique(c(unlist(direction$selection_subject_ids), unlist(direction$selection_unit_ids))), c(tested_subjects, obs$biological_unit_id))
      if (!length(unlist(direction$selection_subject_ids)) || length(overlap)) stop(sprintf("E_RESCUE_DIRECTION_NOT_INDEPENDENT: direction resource selection participants are unknown or overlap tested units (%s)", paste(overlap, collapse = ", ")), call. = FALSE)
      dir_map <- stats::setNames(vapply(direction$directions, function(x) as.numeric(x$direction), 0), vapply(direction$directions, function(x) x$feature_id, ""))
      s <- dir_map[rownames(values)]
      ok <- tested & is.finite(s) & s %in% c(-1, 1)
      comp <- rescue_components(s, d, as.numeric(D$effect_se), t, as.numeric(Tt$effect_se), r, as.numeric(R$effect_se), df_total, p$disease_margin, p$treatment_margin, p$equivalence_margin)
      comp[!ok, ] <- NA
      rr <- cbind(.rs_envelope(n, context, "FormalRescueResult", "formal_rescue", axis$rescue_family, axis$model_id),
                  data.frame(axis_id = axis$id, feature_id = rownames(values), eligibility = ifelse(ok, "tested", "excluded"),
                             reason_code = ifelse(ok, NA_character_, ifelse(tested, "no_independent_direction", "axis_contrast_not_tested")), direction = s,
                             direction_resource_sha256 = .pc_input_sha(request, "direction_resource"), comp, stringsAsFactors = FALSE))
      rr$q_value <- .rs_family_adjust(rr$p_value, axis$rescue_adjustment)
      rescue_rows[[length(rescue_rows) + 1L]] <- rr
    }
    summary$axes[[axis$id]] <- list(n_tested = sum(tested), classes = as.list(table(classes[tested])), shared_control = shared)
  }
  if (length(descriptive_rows)) write_tsv(do.call(rbind, descriptive_rows), "descriptive.tsv", "response_descriptive", "DescriptiveResponse")
  if (length(eq_rows)) write_tsv(do.call(rbind, eq_rows), "equivalence.tsv", "response_equivalence", "EquivalenceResult")
  if (length(rescue_rows)) write_tsv(do.call(rbind, rescue_rows), "formal_rescue.tsv", "response_formal_rescue", "FormalRescueResult")

  # scores
  score_tests <- list(); descriptive_scores <- list()
  for (score in p$scores) {
    manifest_path <- .pc_find_input(request, score$resource_artifact)
    manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
    groups <- unlist(score$required_groups)
    keep <- obs[[p$group_column]] %in% groups
    subjects <- unique(ifelse(is.na(obs$subject_id[keep]), obs$biological_unit_id[keep], obs$subject_id[keep]))
    eligibility <- score_eligibility(manifest, subjects, obs$biological_unit_id[keep], p$score_test)
    sc <- apply_fixed_score(values[, keep, drop = FALSE], manifest)
    sc$group <- obs[[p$group_column]][keep]
    transform_hash <- .pc_input_sha(request, score$resource_artifact)
    if (!isTRUE(eligibility$inferential)) {
      descriptive_scores[[length(descriptive_scores) + 1L]] <- data.frame(schema_version = "1.2.0", result_type = "DescriptiveScore", hypothesis_type = "descriptive_score",
        family_id = NA_character_, score_id = score$id, observation_id = sc$observation_id, group = sc$group, value = sc$value, observed_feature_fraction = sc$observed_feature_fraction,
        transform_sha256 = transform_hash, score_eligibility_reason = eligibility$reason, stringsAsFactors = FALSE)
      next
    }
    ok <- is.finite(sc$value)
    a <- sc$group[ok] == groups[1]
    if (identical(score$scheme, "independent_labels")) {
      res <- if (identical(p$score_test, "independent_exact")) exact_label_randomization(sc$value[ok], a) else monte_carlo_label_randomization(sc$value[ok], a, as.integer(score$draws), as.integer(request$rng$seed))
    } else if (identical(score$scheme, "paired_sign_flip")) {
      ids <- obs[keep, , drop = FALSE]
      subj <- ids$subject_id
      if (any(is.na(subj))) stop("E_SCORE_EXCHANGEABILITY_UNSUPPORTED: paired sign flips need subject IDs", call. = FALSE)
      first <- stats::setNames(sc$value[sc$group == groups[1]], subj[sc$group == groups[1]]); second <- stats::setNames(sc$value[sc$group == groups[2]], subj[sc$group == groups[2]])
      common <- intersect(names(first), names(second))
      res <- exact_sign_flip(first[common] - second[common])
    } else stop(sprintf("E_SCORE_EXCHANGEABILITY_UNSUPPORTED: scheme %s is not exposed in this version", score$scheme), call. = FALSE)
    score_tests[[length(score_tests) + 1L]] <- data.frame(schema_version = "1.2.0", run_id = context$run_id, plan_hash = context$plan_hash, result_type = "IndependentScoreTest",
      hypothesis_type = "independent_score_test", family_id = score$family_id, score_id = score$id, contrast_id = score$contrast_id, scheme = score$scheme,
      exchangeability_evidence = score$exchangeability_evidence, method = res$method, statistic = res$statistic, k = res$k,
      N = if (is.null(res$N)) NA_integer_ else res$N, B = if (is.null(res$B)) NA_integer_ else res$B,
      observed_allocation_count = if (is.null(res$observed_allocation_count)) NA_integer_ else res$observed_allocation_count,
      forced_observed = if (is.null(res$forced_observed)) NA_integer_ else res$forced_observed, p_value = res$p_value, q_value = NA_real_,
      tie_rule = "|T| >= |T_obs| - (1e-12 + 1e-10*|T_obs|)", seed = if (is.null(res$seed)) NA_integer_ else res$seed,
      sampled_tail_ci_lower = if (is.null(res$sampled_tail_ci_lower)) NA_real_ else res$sampled_tail_ci_lower,
      sampled_tail_ci_upper = if (is.null(res$sampled_tail_ci_upper)) NA_real_ else res$sampled_tail_ci_upper, transform_sha256 = transform_hash, stringsAsFactors = FALSE)
  }
  if (length(score_tests)) {
    tests <- do.call(rbind, score_tests)
    for (fam in unique(tests$family_id)) { sel <- tests$family_id == fam; tests$q_value[sel] <- adjust_pvalues(tests$p_value[sel], "BH") }
    write_tsv(tests, "independent_score_tests.tsv", "response_independent_score_tests", "IndependentScoreTest")
  }
  if (length(descriptive_scores)) write_tsv(do.call(rbind, descriptive_scores), "descriptive_scores.tsv", "response_descriptive_scores", "DescriptiveScore")
  summary$scores <- lapply(p$scores, function(s) s$id)
  write_json(c(summary, list(result_type = "ResponseSummary", response_mode = p$response_mode, dmin = p$dmin, ratio_uncertainty = p$ratio_uncertainty,
                             vocabulary = "SM20 frozen vocabulary; classes are descriptive; near_restoration is not equivalence; no causal rescue claim")),
             "response_result.json", "response_result", "ResponseSummary")
  list(outputs = outputs, warnings = warnings, message = sprintf("response: %d axis/axes, %d score(s)", length(p$axes), length(p$scores)))
})

.pc_input_sha <- function(request, artifact_id) {
  item <- Filter(function(x) identical(x$artifact_id, artifact_id), request$inputs)
  if (!length(item)) NA_character_ else item[[1L]]$sha256
}

# Plan artifacts must match the frozen plan; post-plan inputs (dea tables, resources) carry their own verified hashes.
.pm_verify_inputs_response <- function(request) {
  plan_path <- .pc_find_input(request, "plan")
  plan <- jsonlite::fromJSON(plan_path, simplifyVector = FALSE)
  if (is.null(request$plan_hash) || !identical(plan$plan_hash, request$plan_hash)) stop("E_PLAN_CHANGED: request plan_hash does not match the frozen plan", call. = FALSE)
  planned <- list(); for (a in plan$artifacts) planned[[a$artifact_id]] <- a
  resources <- list(); for (r in plan$config$resources) resources[[r$id]] <- r
  for (item in request$inputs) {
    # dea_/moderation_ inputs come from completed post-plan stages (hash-checked by the orchestrator and .pc_find_input);
    # resfile_ inputs are verified against their hash-verified snapshot manifest by the consuming handler.
    if (item$artifact_id == "plan" || startsWith(item$artifact_id, "dea_") || startsWith(item$artifact_id, "moderation_") || startsWith(item$artifact_id, "resfile_") ||
        startsWith(item$artifact_id, "gene_")) next
    if (startsWith(item$artifact_id, "counts_")) {
      model <- Filter(function(m) identical(m$model_id, sub("^counts_", "", item$artifact_id)), plan$models)
      if (!length(model) || !identical(model[[1]]$count_evidence_sha256, item$sha256)) stop(sprintf("E_PLAN_CHANGED: count evidence %s differs from the frozen plan", item$artifact_id), call. = FALSE)
      next
    }
    if (startsWith(item$artifact_id, "resource_") || identical(item$artifact_id, "direction_resource")) {
      id <- sub("^resource_", "", item$artifact_id)
      if (identical(item$artifact_id, "direction_resource")) id <- plan$config$response$direction_resource_id
      ref <- resources[[id]]
      if (is.null(ref) || !identical(ref$sha256, item$sha256)) stop(sprintf("E_RESOURCE_HASH: resource %s differs from the frozen plan", id), call. = FALSE)
      next
    }
    ref <- planned[[item$artifact_id]]
    if (is.null(ref) || !identical(ref$sha256, item$sha256)) stop(sprintf("E_PLAN_CHANGED: input %s differs from the frozen plan", item$artifact_id), call. = FALSE)
  }
  plan
}
