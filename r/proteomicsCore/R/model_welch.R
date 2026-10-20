# Welch engine stage (packet R15a: FR-168 to FR-173, SM42/SM43, V168, V170-V173).
#
# Per feature and per declared pairwise contrast: a two-sample Welch t on the
# OBSERVED log2 cells of the primary matrix (SM04, no new imputation), with
# Welch-Satterthwaite df and effect = mean(numerator) - mean(denominator).
# A design that pairs the observations (subject blocking) uses the paired t on
# within-subject differences with n - 1 df.  The Welch result never uses
# limma moderation and carries its own family.  BH is applied inside each
# declared family through the shared family machinery (adjust_family), and
# the endpoint rows use the shared differential envelope, so verification
# and downstream readers treat them like limma rows with engine = "welch".

# Two-sample Welch statistics for one feature.  x and y are the observed cells
# of the numerator and denominator groups (NA = not observed).
welch_two_sample_row <- function(x, y, ci_level) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  nx <- length(x); ny <- length(y)
  if (nx < 2L || ny < 2L) return(welch_empty_row(nx, ny))
  vx <- stats::var(x) / nx; vy <- stats::var(y) / ny
  se <- sqrt(vx + vy)
  effect <- mean(x) - mean(y)
  df <- (vx + vy)^2 / (vx^2 / (nx - 1) + vy^2 / (ny - 1))
  welch_result_row(effect, se, df, ci_level, nx, ny, "welch_t")
}

# Paired statistics on complete pairs (within-subject differences x - y, n - 1 df).
welch_paired_row <- function(x, y, ci_level) {
  keep <- !is.na(x) & !is.na(y)
  d <- x[keep] - y[keep]
  n <- length(d)
  if (n < 2L) return(welch_empty_row(n, n))
  effect <- mean(d)
  se <- stats::sd(d) / sqrt(n)
  welch_result_row(effect, se, n - 1, ci_level, n, n, "paired_t")
}

welch_empty_row <- function(nx, ny) c(effect = NA_real_, se = NA_real_, t = NA_real_, df = NA_real_, p = NA_real_, lower = NA_real_, upper = NA_real_,
                                      n_num = nx, n_den = ny, zero_variance = FALSE)

welch_result_row <- function(effect, se, df, ci_level, nx, ny, kind) {
  if (!is.finite(se) || se <= 0) {
    out <- welch_empty_row(nx, ny); out["zero_variance"] <- TRUE; return(out)   # numerical failure, never P = 1
  }
  t <- effect / se
  p <- 2 * stats::pt(-abs(t), df)
  q <- stats::qt(1 - (1 - ci_level) / 2, df)
  c(effect = effect, se = se, t = t, df = df, p = p, lower = effect - q * se, upper = effect + q * se, n_num = nx, n_den = ny, zero_variance = FALSE)
}

# Subject pairs for a paired contrast: one numerator and one denominator observation per subject.
# Returns NULL when the design cannot be honoured as a within-subject comparison.
welch_pairs <- function(obs, group, numerator, denominator, subject_column) {
  if (is.null(subject_column) || !subject_column %in% names(obs)) return(NULL)
  subject <- obs[[subject_column]]
  idx_num <- which(group == numerator); idx_den <- which(group == denominator)
  s_num <- subject[idx_num]; s_den <- subject[idx_den]
  if (anyNA(s_num) || anyNA(s_den) || anyDuplicated(s_num) || anyDuplicated(s_den)) return(NULL)
  common <- intersect(s_num, s_den)
  if (length(common) < 2L) return(NULL)
  list(num = idx_num[match(common, s_num)], den = idx_den[match(common, s_den)], n_unpaired = length(s_num) + length(s_den) - 2L * length(common))
}

welch_contrast_matrix <- function(Y, num_idx, den_idx, pairs, ci_level) {
  stats_rows <- lapply(seq_len(nrow(Y)), function(i) {
    if (is.null(pairs)) welch_two_sample_row(Y[i, num_idx], Y[i, den_idx], ci_level)
    else welch_paired_row(Y[i, pairs$num], Y[i, pairs$den], ci_level)
  })
  out <- do.call(rbind, stats_rows)
  rownames(out) <- rownames(Y)
  out
}

# Decides the variant (two-sample, paired, or an adapted unadjusted two-sample) for one model.  The refusals are
# the typed validity errors of ADR 0010; adaptations are recorded, never silent.
welch_variant <- function(model, obs, group, numerator, denominator) {
  blocked <- isTRUE(model$blocked)
  declared <- if (is.null(model$paired)) "auto" else model$paired
  if (identical(declared, "false") && blocked) stop("E_WELCH_PAIRED_UNIT: a declared paired or subject-blocked design cannot run the two-sample Welch form", call. = FALSE)
  if (identical(declared, "true") && !blocked) stop("E_WELCH_PAIRED_UNIT: paired = true requires a subject-blocked design", call. = FALSE)
  if (!blocked) return(list(variant = "two_sample", pairs = NULL, adaptation = NULL))
  pairs <- welch_pairs(obs, group, numerator, denominator, model$subject_column)
  if (!is.null(pairs)) return(list(variant = "paired", pairs = pairs, adaptation = NULL))
  if (identical(declared, "true")) stop("E_WELCH_PAIRED_UNIT: the subject structure does not pair the two groups", call. = FALSE)
  list(variant = "two_sample", pairs = NULL, adaptation = list(analysis = "welch", item = paste0("block:", model$subject_column), requested = "within_subject",
                                                              used = "unadjusted", reason = "welch_block_not_pairable", model_id = model$model_id, contrast_id = NA_character_))
}

welch_stage <- function(request) .pc_run_stage(request, function(out) {
  p <- request$parameters
  outputs <- list(); warnings <- list()
  emit <- function(relative, id, type) outputs[[length(outputs) + 1L]] <<- .pc_output(out, relative, id, type)
  write_tsv <- function(df, relative, id, type) { .pc_write_tsv(df, file.path(out, relative)); emit(relative, id, type) }
  write_json <- function(value, relative, id, type) { .pc_write_json(value, file.path(out, relative)); emit(relative, id, type) }
  .pc_verify_plan_inputs(request)
  context <- list(run_id = request$run_id, plan_hash = request$plan_hash)
  observed <- .pc_matrix_from_tsv(.pc_find_input(request, "primary_observed_mask"), "logical")
  values_all <- .pc_matrix_from_tsv(.pc_find_input(request, "primary_matrix"), "numeric")
  observations <- .pc_read_tsv(.pc_find_input(request, "primary_observations"))
  estimability <- .pc_read_tsv(.pc_find_input(request, "estimability"))
  estimability$estimable <- estimability$estimable == "true"
  ci_level <- if (is.null(p$ci_level)) 0.95 else as.numeric(p$ci_level)
  obs_index <- match(colnames(values_all), observations$observation_id)
  if (anyNA(obs_index)) stop("E_ID_ALIGNMENT: primary matrix columns are not all listed in primary_observations", call. = FALSE)
  obs <- observations[obs_index, , drop = FALSE]
  group <- obs[[p$group_column]]
  values_all[!observed[rownames(values_all), colnames(values_all), drop = FALSE]] <- NA_real_
  if ("biological_unit_id" %in% names(obs) && anyDuplicated(obs$biological_unit_id[!is.na(obs$biological_unit_id)]) && !any(vapply(p$models, function(m) isTRUE(m$blocked), logical(1))))
    stop("E_WELCH_BIOLOGICAL_UNIT: unblocked Welch requires one observation per biological unit (SM03); technical replicates are not independent", call. = FALSE)

  rows <- list(); adaptations <- list(); decisions <- list()
  for (model in p$models) {
    if (length(model$refused)) for (r in model$refused) adaptations[[length(adaptations) + 1L]] <- list(analysis = "welch", item = r$item, requested = r$requested, used = "not_run", reason = r$reason, model_id = model$model_id, contrast_id = r$contrast_id %||% NA_character_, refusal = TRUE)
    for (a in model$adaptations) {
      adaptations[[length(adaptations) + 1L]] <- c(list(model_id = model$model_id, contrast_id = NA_character_), a)
      warnings[[length(warnings) + 1L]] <- .pc_warning(request, "W_WELCH_ADAPTED", sprintf("welch model %s: %s requested %s, used %s (%s)", model$model_id, a$item, a$requested, a$used, a$reason),
                                                       "methods_extension/welch/eligibility.json")
    }
    for (contrast in model$contrasts) {
      num_idx <- which(group == contrast$numerator); den_idx <- which(group == contrast$denominator)
      variant <- welch_variant(model, obs, group, contrast$numerator, contrast$denominator)
      if (!is.null(variant$adaptation)) {
        variant$adaptation$contrast_id <- contrast$contrast_id
        adaptations[[length(adaptations) + 1L]] <- variant$adaptation
        warnings[[length(warnings) + 1L]] <- .pc_warning(request, "W_WELCH_ADAPTED", sprintf("welch model %s contrast %s: %s requested, %s used (%s)", model$model_id, contrast$contrast_id,
                                                                                         variant$adaptation$requested, variant$adaptation$used, variant$adaptation$reason), "methods_extension/welch/eligibility.json")
      }
      n_pairs <- if (is.null(variant$pairs)) NA_integer_ else length(variant$pairs$num)
      decisions[[length(decisions) + 1L]] <- list(model_id = model$model_id, contrast_id = contrast$contrast_id, numerator = contrast$numerator, denominator = contrast$denominator,
                                                  variant = variant$variant, n_pairs = n_pairs, n_unpaired = if (is.null(variant$pairs)) NA_integer_ else variant$pairs$n_unpaired)
      est <- estimability[estimability$model_id == model$model_id & estimability$contrast_id == contrast$contrast_id, , drop = FALSE]
      est <- est[match(rownames(values_all), est$feature_id), , drop = FALSE]
      absent <- is.na(est$feature_id)
      est$eligibility[absent] <- "excluded"; est$reason[absent] <- "absent_from_estimability"; est$estimable[absent] <- FALSE; est$n_obs_by_required_group[absent] <- NA_character_
      est$feature_id <- rownames(values_all)
      stats_m <- welch_contrast_matrix(values_all, num_idx, den_idx, variant$pairs, ci_level)
      out_rows <- endpoint_frame(nrow(values_all))
      out_rows$schema_version <- "1.2.0"; out_rows$run_id <- context$run_id; out_rows$plan_hash <- context$plan_hash
      out_rows$result_type <- "WelchResult"; out_rows$hypothesis_type <- "protein_zero_null"; out_rows$engine <- "welch"
      out_rows$engine_version <- as.character(getRversion())
      out_rows$model_id <- model$model_id; out_rows$design_id <- model$design_id; out_rows$contrast_id <- contrast$contrast_id; out_rows$feature_id <- rownames(values_all)
      out_rows$role <- if (identical(model$role, "primary")) "primary" else "secondary"
      out_rows$estimable <- as.logical(est$estimable); out_rows$n_obs_by_required_group <- est$n_obs_by_required_group
      out_rows$eligibility <- ifelse(est$eligibility == "eligible", "tested", est$eligibility)
      out_rows$reason_code <- ifelse(est$eligibility == "eligible", NA_character_, est$reason)
      tested <- out_rows$eligibility == "tested"
      small <- tested & (stats_m[, "n_num"] < 2 | stats_m[, "n_den"] < 2)
      out_rows$eligibility[small] <- "nonestimable"; out_rows$reason_code[small] <- "insufficient_units_per_group"
      tested <- out_rows$eligibility == "tested"
      broken <- tested & !(is.finite(stats_m[, "effect"]) & is.finite(stats_m[, "t"]) & is.finite(stats_m[, "p"]) & !as.logical(stats_m[, "zero_variance"]))
      out_rows$eligibility[broken] <- "numerical_failure"; out_rows$reason_code[broken] <- "E_NUMERICAL_NONFINITE_STATISTIC"
      tested <- out_rows$eligibility == "tested"
      out_rows$effect[tested] <- stats_m[tested, "effect"]; out_rows$effect_se[tested] <- stats_m[tested, "se"]
      out_rows$statistic[tested] <- stats_m[tested, "t"]; out_rows$statistic_type[tested] <- if (is.null(variant$pairs)) "Welch t" else "paired t"
      out_rows$df_residual[tested] <- stats_m[tested, "df"]; out_rows$df_inference[tested] <- stats_m[tested, "df"]
      out_rows$p_value[tested] <- stats_m[tested, "p"]
      out_rows$ci_lower[tested] <- stats_m[tested, "lower"]; out_rows$ci_upper[tested] <- stats_m[tested, "upper"]; out_rows$ci_level[tested] <- ci_level
      out_rows$ci_method[tested] <- if (is.null(variant$pairs)) "effect +/- qt(1-(1-level)/2, df) * sqrt(var_num/n_num + var_den/n_den)" else "paired: mean difference +/- qt(1-(1-level)/2, n-1) * sd(difference)/sqrt(n)"
      out_rows$exactness_path[tested] <- if (is.null(variant$pairs)) "welch_closed_form" else "paired_closed_form"
      out_rows$effect_scale <- "log2_difference"; out_rows$null_region <- ifelse(tested, "log2 effect = 0", NA_character_)
      out_rows$mean_abundance[tested] <- vapply(which(tested), function(i) mean(values_all[i, c(num_idx, den_idx)], na.rm = TRUE), numeric(1))
      out_rows$effect_threshold <- NA_real_; out_rows$display_effect_filter <- NA
      rows[[length(rows) + 1L]] <- out_rows
    }
    # A declared family member that cannot be oriented as a pairwise contrast stays a planned endpoint: not_run, with its reason (ADR 0010).
    for (refusal in model$refused) {
      nr <- endpoint_frame(nrow(values_all))
      nr$schema_version <- "1.2.0"; nr$run_id <- context$run_id; nr$plan_hash <- context$plan_hash
      nr$result_type <- "WelchResult"; nr$hypothesis_type <- "protein_zero_null"; nr$engine <- "welch"; nr$engine_version <- as.character(getRversion())
      nr$model_id <- model$model_id; nr$design_id <- model$design_id; nr$contrast_id <- refusal$contrast_id; nr$feature_id <- rownames(values_all)
      nr$role <- if (identical(model$role, "primary")) "primary" else "secondary"
      nr$eligibility <- "not_run"; nr$reason_code <- refusal$reason; nr$estimable <- FALSE
      rows[[length(rows) + 1L]] <- nr
    }
  }
  if (!length(rows)) stop("E_NO_ELIGIBLE_FEATURES: no Welch model or contrast was requested", call. = FALSE)
  combined <- do.call(rbind, rows)
  family_summaries <- list()
  for (family in p$families) {
    spec <- list(family_id = family$family_id, hypothesis_type = "protein_zero_null", role = family$role, adjustment = family$adjustment,
                 dependence_assumption = family$dependence_assumption, q_cutoff = as.numeric(family$q_cutoff),
                 members = lapply(family$members, function(m) list(model_id = m$model_id, contrast_id = m$contrast_id)))
    adjusted <- adjust_family(combined, spec)
    combined <- adjusted$rows; family_summaries[[length(family_summaries) + 1L]] <- adjusted$summary
  }
  combined$native_q_value <- combined$q_value
  planned <- unlist(lapply(rows, function(r) paste(r$model_id, r$contrast_id, r$feature_id, r$hypothesis_type, sep = "|")))
  verify_differential_table(combined, planned)
  write_tsv(combined, "welch_results.tsv", "welch_results", "WelchResult")
  write_tsv(do.call(rbind, family_summaries), "families.tsv", "welch_families", "HypothesisFamily")
  # Welch as the primary engine writes the dea/ contract (amendment A-2026-10-01-46): the same zero-null table name, families and model status limma writes
  if (isTRUE(p$primary_dea)) {
    write_tsv(combined[combined$hypothesis_type == "protein_zero_null", , drop = FALSE], "zero_null.tsv", "dea_zero_null", "ProteinZeroNullResult")
    write_json(list(models = lapply(p$models, function(m) list(model_id = m$model_id, state = "COMPLETED", required = identical(m$execution_requirement, "required"), reason_code = NULL)),
                    settings = setNames(lapply(p$models, function(m) list(input_matrix = "primary_matrix", observed_cells_only = TRUE, input_imputation = "none",
                                                                          new_primary_imputation = "none", engine = "welch", moderation = "none",
                                                                          paired = m$paired, blocked = isTRUE(m$blocked))), vapply(p$models, function(m) m$model_id, "")),
                    score_p_values_emitted = FALSE, note = "Welch primary engine: dea/ holds the Welch tables; no limma moderation is applied."),
               "model_status.json", "dea_model_status", "ModelFit")
  }
  pc <- p$paper_compatible
  paper_rows <- 0L
  if (isTRUE(pc$enabled)) {
    cut <- as.numeric(pc$raw_p_cutoff)
    sel <- combined$eligibility == "tested" & !is.na(combined$p_value) & combined$p_value < cut
    paper <- data.frame(model_id = combined$model_id[sel], contrast_id = combined$contrast_id[sel], feature_id = combined$feature_id[sel],
                        effect = combined$effect[sel], statistic = combined$statistic[sel], raw_p_value = combined$p_value[sel],
                        bh_q_value = combined$q_value[sel], family_id = combined$family_id[sel], raw_p_cutoff = cut,
                        selection_rule = "raw_p_below_cutoff", claim_label = "not_multiplicity_controlled", stringsAsFactors = FALSE)
    paper_rows <- nrow(paper)
    write_tsv(paper, "paper_compatible.tsv", "welch_paper_compatible", "PaperCompatibleSelection")
  }
  write_json(list(result_type = "WelchEligibility", engine = "welch", claim_labels = I(c("welch", if (isTRUE(pc$enabled)) "not_multiplicity_controlled")),
                  method = "two-sample Welch t (Welch-Satterthwaite df) on observed log2 cells; paired t on within-subject differences when the design pairs observations",
                  decisions = decisions, adaptations = adaptations, paper_compatible = list(enabled = isTRUE(pc$enabled), raw_p_cutoff = if (isTRUE(pc$enabled)) as.numeric(pc$raw_p_cutoff) else NULL, selected = paper_rows),
                  families = lapply(p$families, function(f) f$family_id), no_new_imputation = TRUE, limma_moderation_applied = FALSE,
                  versions = list(R = as.character(getRversion()))),
             "eligibility.json", "welch_eligibility", "WelchEligibility")
  list(outputs = outputs, warnings = warnings,
       message = sprintf("welch: %d endpoint rows across %d model-contrast pairs; %d adaptation record(s)", nrow(combined), length(decisions), length(adaptations)))
})
