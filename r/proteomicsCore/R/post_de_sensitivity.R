# Post-DE robustness and sensitivity (packet R14b, SM32, FR-137 to FR-143).
#
# Covariate-adjusted and subgroup models are named secondary sensitivity models: their designs were re-planned in
# Python through the R04 grammar with exact rank/estimability (design_service.replan); here the rank is re-verified
# with design_rank, SM04 coverage and featurewise estimability are re-applied on the analysed observations, and the
# model is fitted with the primary engine settings (fit_limma_model).  Every model has its own families; the primary
# families and tables are only read.  Matched-n resampling and influence act on biological units (subjects in
# subject-blocked designs) with all their observations.  Robustness summaries are descriptive.

PD_INTERACTION_WARNING <- "Differing significance between the primary and a sensitivity or subgroup analysis is not evidence of an interaction; only the formal interaction contrast (where estimable) tests effect modification."

# Fit one named model on the given observations: coverage (SM04) -> featurewise estimability (SM07) -> limma (primary settings) -> families.
pd_fit_model <- function(model_id, values, observed, obs, X, contrasts, settings, families, context, coverage_rule) {
  keep <- colnames(values)
  rules <- lapply(contrasts, function(c) list(model_id = model_id, contrast_id = c$contrast_id, required_groups = c$required_groups,
                                              policy = coverage_rule$policy, minimum_observed_per_group = coverage_rule$minimum_observed_per_group,
                                              minimum_fraction = coverage_rule$minimum_fraction))
  cov <- coverage_tables(observed[, keep, drop = FALSE], obs, settings$group_column, rules)
  groups <- obs[[settings$group_column]]
  est <- featurewise_estimability(values, observed[, keep, drop = FALSE], X, contrasts, cov$summary, model_id, groups)
  est$estimable <- as.logical(est$estimable)
  model <- list(model_id = model_id, design_id = model_id, role = "sensitivity", hypotheses = list("zero_null"), trend = settings$trend, robust = settings$robust,
                effect_threshold = NULL)
  Y <- values; Y[!observed[, keep, drop = FALSE]] <- NA
  block <- NULL; correlation <- NULL
  if (identical(settings$blocking_mode, "duplicate_correlation")) { block <- obs[[settings$subject_column]]; correlation <- settings$consensus_correlation }
  fit <- fit_limma_model(model, Y, X, contrasts, est, obs, NULL, block, correlation, settings$ci_level, context)
  rows <- fit$rows; summaries <- list()
  for (family in families) { adj <- adjust_family(rows, family); rows <- adj$rows; summaries[[length(summaries) + 1L]] <- adj$summary }
  list(rows = rows, families = if (length(summaries)) do.call(rbind, summaries) else data.frame(family_id = character()), fit = fit)
}

# Families of a sensitivity model mirror the primary model's zero-null families (same adjustment and cutoff), restricted to its contrasts.
pd_mirror_families <- function(primary_families, model_id, contrast_ids, suffix) {
  out <- list()
  for (f in primary_families) {
    keep <- Filter(function(m) m$contrast_id %in% contrast_ids, f$members)
    if (!length(keep)) next
    out[[length(out) + 1L]] <- list(family_id = paste0(f$family_id, "__", suffix), hypothesis_type = "protein_zero_null", role = "sensitivity",
      adjustment = f$adjustment, q_cutoff = f$q_cutoff, dependence_assumption = f$dependence_assumption,
      members = lapply(keep, function(m) list(model_id = model_id, contrast_id = m$contrast_id)))
  }
  out
}

.pd_design_from <- function(spec) {
  X <- matrix(as.numeric(unlist(spec$matrix)), nrow = length(spec$observation_ids), byrow = TRUE, dimnames = list(unlist(spec$observation_ids), unlist(spec$coefficients)))
  X
}

.pd_members <- function(rows, contrast_id, criterion) {
  r <- rows[rows$contrast_id == contrast_id & rows$eligibility == "tested", , drop = FALSE]
  value <- if (identical(criterion$criterion, "exploratory_raw_p")) r$p_value else r$q_value
  m <- !is.na(value) & value < criterion$threshold
  d <- if (is.null(criterion$direction)) "any" else criterion$direction
  if (identical(d, "up")) m <- m & r$effect > 0
  if (identical(d, "down")) m <- m & r$effect < 0
  r$feature_id[m]
}

.pd_effects <- function(rows, contrast_id) { r <- rows[rows$contrast_id == contrast_id & rows$eligibility == "tested", , drop = FALSE]; stats::setNames(r$effect, r$feature_id) }

# Primary-versus-sensitivity comparison (FR-140).
pd_compare <- function(primary_rows, sens_rows, contrast_id, criterion) {
  a <- .pd_members(primary_rows, contrast_id, criterion); b <- .pd_members(sens_rows, contrast_id, criterion)
  ea <- .pd_effects(primary_rows, contrast_id); eb <- .pd_effects(sens_rows, contrast_id)
  common <- intersect(names(ea), names(eb))
  x <- ea[common]; y <- eb[common]
  list(retained = intersect(a, b), lost = setdiff(a, b), gained = setdiff(b, a), n_common_tested = length(common),
       effect_correlation = if (length(common) > 2L) stats::cor(x, y) else NA_real_,
       attenuation_slope = if (length(common)) sum(x * y) / sum(x * x) else NA_real_,
       jaccard = if (length(union(a, b))) length(intersect(a, b)) / length(union(a, b)) else NA_real_)
}

# Matched-n draws: biological units of each targeted group, sampled without replacement in declared group order.
pd_matched_draws <- function(units, groups, targets, draws, seed) {
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  unit_group <- tapply(groups, units, function(g) g[1]); unit_order <- unique(units)
  lapply(seq_len(draws), function(d) {
    chosen <- character()
    for (g in names(targets)) {
      pool <- unit_order[unit_group[unit_order] == g]
      chosen <- c(chosen, sample(pool, as.integer(targets[[g]])))
    }
    untouched <- unit_order[!unit_group[unit_order] %in% names(targets)]
    c(chosen, untouched)
  })
}

.pd_fig_scatter <- function(out, stem, x, y, main, formats, xl, yl) .pm_devices(out, stem, formats, 5.5, 5.5, function() {
  graphics::par(mar = c(4.5, 4.5, 3, 1)); graphics::plot(x, y, pch = 19, col = "grey35", xlab = xl, ylab = yl, main = main, cex.main = 0.85)
  graphics::abline(0, 1, lty = 3); graphics::abline(h = 0, v = 0, col = "grey80") })

post_de_sensitivity_stage <- function(request) .pc_run_stage(request, function(out) {
  p <- request$parameters
  outputs <- list(); warnings <- list(); refusals <- list(); figures <- list()
  emit <- function(relative, id, type) outputs[[length(outputs) + 1L]] <<- .pc_output(out, relative, id, type)
  write_tsv <- function(df, relative, id, type) { .pc_write_tsv(df, file.path(out, relative)); emit(relative, id, type) }
  write_json <- function(value, relative, id, type) { .pc_write_json(value, file.path(out, relative)); emit(relative, id, type) }
  add_fig <- function(x) { for (f in x$files) emit(f$relative_path, f$artifact_id, "Figure"); figures <<- c(figures, x$records) }
  refuse <- function(analysis, item, code, reason) {
    refusals[[length(refusals) + 1L]] <<- data.frame(analysis = analysis, item = item, reason_code = code, reason = reason, stringsAsFactors = FALSE)
    warnings[[length(warnings) + 1L]] <<- .pc_warning(request, code, sprintf("%s %s refused: %s", analysis, item, reason), "post_de/sensitivity/refusals.tsv")
  }
  plan <- .pd_verify_inputs(request)
  warnings <- c(warnings, .pd_adaptation_warnings(request))
  prim <- .pd_primary(request)
  values <- .pc_matrix_from_tsv(.pc_find_input(request, "primary_matrix"), "numeric")
  observed <- .pc_matrix_from_tsv(.pc_find_input(request, "primary_observed_mask"), "logical")[rownames(values), colnames(values), drop = FALSE]
  obs <- prim$obs
  zero <- .pd_read_stage_table(request, "stage__dea_zero_null")
  primary_rows <- zero[zero$model_id == p$primary_model_id & zero$hypothesis_type == "protein_zero_null", , drop = FALSE]
  context <- list(run_id = request$run_id, plan_hash = request$plan_hash)
  settings <- list(group_column = p$group_column, subject_column = p$subject_column, blocking_mode = p$blocking_mode, consensus_correlation = p$consensus_correlation,
                   trend = isTRUE(p$trend), robust = isTRUE(p$robust), ci_level = p$ci_level)
  criteria <- p$criteria
  primary_families <- Filter(function(f) identical(f$hypothesis_type, "protein_zero_null") && any(vapply(f$members, function(m) identical(m$model_id, p$primary_model_id), logical(1))), plan$families)
  lineage_path <- .pc_find_input(request, "observation_lineage", required = FALSE)
  lineage <- if (is.null(lineage_path)) NULL else .pc_read_tsv(lineage_path)
  injections_of <- function(ids) {
    if (is.null(lineage)) return(ids)
    col <- intersect(c("canonical_observation_id", "observation_id"), names(lineage))[1]
    src <- intersect(c("source_observation_id", "source_column", "injection_id"), names(lineage))[1]
    if (is.na(col) || is.na(src)) return(ids)
    unique(lineage[[src]][lineage[[col]] %in% ids])
  }

  # FR-137 imbalance diagnostics and estimability classes (decided by the planner; re-verified here)
  write_tsv(.pd_rows(p$imbalance, c("covariate", "level", "group", "n_units")), "imbalance.tsv", "post_de_imbalance", "SensitivityComparison")
  write_tsv(.pd_rows(p$classification, c("model_id", "kind", "class", "reason_code", "detail")), "estimability_classes.tsv", "post_de_estimability_classes", "SensitivityComparison")

  model_rows <- list(); comparison_rows <- list(); fitted <- list(); sens_tables <- list()
  for (m in p$models) {
    if (!identical(m$state, "ELIGIBLE")) { refuse(m$kind, m$model_id, m$reason_code, m$reason); next }
    X <- .pd_design_from(m$design)
    rk <- design_rank(X)
    if (!identical(as.integer(rk$rank), as.integer(m$rank))) stop(sprintf("E_INTEGRITY: model %s rank %d differs from the planned exact rank %d", m$model_id, rk$rank, m$rank), call. = FALSE)
    ids <- rownames(X); o <- obs[match(ids, obs$observation_id), , drop = FALSE]
    interaction <- identical(m$kind, "interaction")
    fams <- if (interaction) list(list(family_id = paste0("interaction__", m$model_id), hypothesis_type = "protein_zero_null", role = "sensitivity",
      adjustment = "BH", q_cutoff = 0.05, dependence_assumption = "BH: independence or positive regression dependence",
      members = lapply(m$contrasts, function(c) list(model_id = m$model_id, contrast_id = c$contrast_id))))
      else pd_mirror_families(primary_families, m$model_id, vapply(m$contrasts, function(c) c$contrast_id, ""), m$model_id)
    res <- pd_fit_model(m$model_id, values[, ids, drop = FALSE], observed[, ids, drop = FALSE], o, X, m$contrasts, settings, fams, context, p$coverage)
    res$rows$role <- "sensitivity"; res$rows$result_type <- if (interaction) "InteractionContrastResult" else "SensitivityModelResult"
    base <- file.path("models", m$model_id)
    write_tsv(res$rows, file.path(base, "zero_null.tsv"), paste0("sensitivity_", m$model_id), "SensitivityModelResult")
    write_tsv(res$families, file.path(base, "families.tsv"), paste0("sensitivity_families_", m$model_id), "HypothesisFamily")
    fitted[[m$model_id]] <- res$rows; sens_tables[[m$model_id]] <- m
    model_rows[[length(model_rows) + 1L]] <- data.frame(model_id = m$model_id, kind = m$kind, n_observations = length(ids), n_units = length(unique(.pd_units(o, settings$subject_column))),
      design_rank = rk$rank, coefficients = paste(colnames(X), collapse = ";"), families = paste(vapply(fams, function(f) f$family_id, ""), collapse = ";"),
      description = m$description, claim_label = "descriptive", stringsAsFactors = FALSE)
    if (interaction) next
    imodel <- paste0(m$model_id, "__interaction")
    ifamily <- if (any(vapply(p$models, function(x) identical(x$model_id, imodel) && identical(x$state, "ELIGIBLE"), logical(1)))) paste0("interaction__", imodel) else NA_character_
    for (c in m$contrasts) for (crit in criteria) {
      cmp <- pd_compare(primary_rows, res$rows, c$contrast_id, crit)
      comparison_rows[[length(comparison_rows) + 1L]] <- data.frame(model_id = m$model_id, kind = m$kind, contrast_id = c$contrast_id,
        criterion = sprintf("%s < %s (%s)", crit$criterion, format(crit$threshold), if (is.null(crit$direction)) "any" else crit$direction),
        n_primary = length(cmp$retained) + length(cmp$lost), n_sensitivity = length(cmp$retained) + length(cmp$gained),
        n_retained = length(cmp$retained), n_lost = length(cmp$lost), n_gained = length(cmp$gained), jaccard = cmp$jaccard,
        n_common_tested = cmp$n_common_tested, effect_correlation = cmp$effect_correlation, attenuation_slope = cmp$attenuation_slope,
        retained = paste(sort(cmp$retained), collapse = ";"), lost = paste(sort(cmp$lost), collapse = ";"), gained = paste(sort(cmp$gained), collapse = ";"),
        interaction_family = ifamily,
        warning = PD_INTERACTION_WARNING, claim_label = "descriptive", stringsAsFactors = FALSE)
      ea <- .pd_effects(primary_rows, c$contrast_id); eb <- .pd_effects(res$rows, c$contrast_id); common <- intersect(names(ea), names(eb))
      stem <- paste0("comparison_", .pd_safe(m$model_id), "_", .pd_safe(c$contrast_id))
      src <- data.frame(feature_id = common, primary_effect = ea[common], sensitivity_effect = eb[common], stringsAsFactors = FALSE)
      if (!identical(crit, criteria[[1]])) next
      write_tsv(src, file.path("figure_sources", paste0(stem, ".tsv")), paste0(stem, "_source"), "FigureSource")
      add_fig(.pd_fig_scatter(out, stem, src$primary_effect, src$sensitivity_effect, sprintf("%s vs primary: %s (descriptive)", m$model_id, c$contrast_id), unlist(p$figure_formats), "primary effect", "sensitivity effect"))
    }
  }
  write_tsv(.pd_rows(model_rows, c("model_id", "kind")), "sensitivity_models.tsv", "post_de_sensitivity_models", "SensitivityModelResult")
  comparison <- .pd_rows(comparison_rows, c("model_id", "contrast_id", "warning", "claim_label"))
  write_tsv(comparison, "comparison.tsv", "post_de_sensitivity_comparison", "SensitivityComparison")

  # FR-141 matched-n resampling on biological units
  X0 <- .pd_design_from(p$primary_design)
  units_all <- .pd_units(obs, settings$subject_column); groups_all <- obs[[settings$group_column]]
  primary_contrasts <- p$primary_contrasts
  draws_df <- NULL; freq_df <- NULL
  if (!is.null(p$matched_n) && identical(p$matched_n$state, "ELIGIBLE")) {
    targets <- p$matched_n$targets
    draws <- pd_matched_draws(units_all, groups_all, targets, as.integer(p$matched_n$draws), as.integer(p$matched_n$seed))
    counts <- list(); hits <- list()
    for (d in seq_along(draws)) {
      keep <- units_all %in% draws[[d]]
      ids <- obs$observation_id[keep]
      fams <- pd_mirror_families(primary_families, "matched_n", vapply(primary_contrasts, function(c) c$contrast_id, ""), "matched_n")
      res <- pd_fit_model("matched_n", values[, ids, drop = FALSE], observed[, ids, drop = FALSE], obs[keep, , drop = FALSE], X0[ids, , drop = FALSE], primary_contrasts, settings, fams, context, p$coverage)
      for (c in primary_contrasts) {
        mem <- .pd_members(res$rows, c$contrast_id, criteria[[1]])
        counts[[length(counts) + 1L]] <- data.frame(draw = d, contrast_id = c$contrast_id, n_units = length(unique(units_all[keep])), n_discoveries = length(mem),
                                                    units = paste(sort(unique(units_all[keep])), collapse = ";"), stringsAsFactors = FALSE)
        hits[[length(hits) + 1L]] <- data.frame(contrast_id = c$contrast_id, feature_id = mem, stringsAsFactors = FALSE)
      }
    }
    draws_df <- do.call(rbind, counts)
    h <- do.call(rbind, hits)
    freq_df <- do.call(rbind, lapply(primary_contrasts, function(c) {
      tested <- primary_rows$feature_id[primary_rows$contrast_id == c$contrast_id & primary_rows$eligibility == "tested"]
      k <- if (is.null(h) || !nrow(h)) integer(length(tested)) else vapply(tested, function(f) sum(h$contrast_id == c$contrast_id & h$feature_id == f), integer(1))
      data.frame(contrast_id = c$contrast_id, feature_id = tested, selection_count = as.integer(k), draws = length(draws), selection_frequency = as.numeric(k) / length(draws),
                 claim_label = "descriptive", stringsAsFactors = FALSE)
    }))
    summary <- data.frame(seed = as.integer(p$matched_n$seed), draws = length(draws), targets = as.character(jsonlite::toJSON(targets, auto_unbox = TRUE)),
                          unit_level = if (is.null(settings$subject_column)) "biological_unit" else paste0("subject:", settings$subject_column),
                          rng_kind = "L'Ecuyer-CMRG", sampling = "without replacement within each targeted group, declared group order", claim_label = "descriptive", stringsAsFactors = FALSE)
    write_tsv(cbind(draws_df, claim_label = "descriptive"), "matched_n.tsv", "post_de_matched_n", "ResamplingResult")
    write_tsv(freq_df, "matched_n_selection.tsv", "post_de_matched_n_selection", "ResamplingResult")
    write_tsv(summary, "matched_n_settings.tsv", "post_de_matched_n_settings", "ResamplingResult")
    stem <- "matched_n_discoveries"
    write_tsv(draws_df[, c("draw", "contrast_id", "n_discoveries")], file.path("figure_sources", paste0(stem, ".tsv")), paste0(stem, "_source"), "FigureSource")
    add_fig(.pm_devices(out, stem, unlist(p$figure_formats), 6, 4.5, function() {
      graphics::par(mar = c(4.5, 4.5, 3, 1))
      graphics::hist(draws_df$n_discoveries, breaks = 20, col = "grey70", border = "white", xlab = "discoveries per draw", main = "Matched-n resampling (descriptive)") }))
  } else if (!is.null(p$matched_n)) refuse("matched_n", "matched_n", p$matched_n$reason_code, p$matched_n$reason)

  # FR-142 unit-level influence (omit each biological unit with all its observations)
  infl_df <- NULL; infl_members <- list()
  if (isTRUE(p$influence)) {
    rows <- list()
    primary_members <- lapply(primary_contrasts, function(c) .pd_members(primary_rows, c$contrast_id, criteria[[1]])); names(primary_members) <- vapply(primary_contrasts, function(c) c$contrast_id, "")
    for (u in unique(units_all)) {
      keep <- units_all != u
      ids <- obs$observation_id[keep]; Xs <- X0[ids, , drop = FALSE]
      nonzero <- colSums(abs(Xs)) > 0; Xs <- Xs[, nonzero, drop = FALSE]
      omitted <- obs$observation_id[!keep]
      state <- "refit"; reason <- NA_character_
      cs <- lapply(primary_contrasts, function(c) { c2 <- c; c2$weights <- as.list(as.numeric(unlist(c$weights))[nonzero]); c2 })
      if (any(as.numeric(unlist(lapply(primary_contrasts, function(c) unlist(c$weights)[!nonzero]))) != 0)) { state <- "ineligible_refit"; reason <- "contrast_uses_omitted_column" }
      else if (qr(Xs)$rank < ncol(Xs)) { state <- "ineligible_refit"; reason <- "design_rank_deficient_after_omission" }
      else if (any(vapply(primary_contrasts, function(c) any(vapply(unlist(c$required_groups), function(g) length(unique(units_all[keep & groups_all == g])) < 2L, logical(1))), logical(1)))) { state <- "ineligible_refit"; reason <- "fewer_than_two_units_in_required_group" }
      res <- NULL
      if (state == "refit") {
        fams <- pd_mirror_families(primary_families, "influence", names(primary_members), "influence")
        res <- pd_fit_model("influence", values[, ids, drop = FALSE], observed[, ids, drop = FALSE], obs[keep, , drop = FALSE], Xs, cs, settings, fams, context, p$coverage)
      }
      for (c in primary_contrasts) {
        cid <- c$contrast_id
        mem <- if (is.null(res)) character() else .pd_members(res$rows, cid, criteria[[1]])
        ea <- .pd_effects(primary_rows, cid); eb <- if (is.null(res)) numeric() else .pd_effects(res$rows, cid)
        common <- intersect(names(ea), names(eb)); delta <- abs(eb[common] - ea[common])
        infl_members[[length(infl_members) + 1L]] <- list(unit = u, contrast_id = cid, members = mem, state = state, effects = eb)
        rows[[length(rows) + 1L]] <- data.frame(unit_id = u, contrast_id = cid, omitted_observation_ids = paste(omitted, collapse = ";"),
          omitted_source_observations = paste(injections_of(omitted), collapse = ";"), n_omitted_observations = length(omitted), state = state, reason = reason,
          n_lost = if (is.null(res)) NA_integer_ else length(setdiff(primary_members[[cid]], mem)), n_gained = if (is.null(res)) NA_integer_ else length(setdiff(mem, primary_members[[cid]])),
          lost = if (is.null(res)) NA_character_ else paste(sort(setdiff(primary_members[[cid]], mem)), collapse = ";"),
          gained = if (is.null(res)) NA_character_ else paste(sort(setdiff(mem, primary_members[[cid]])), collapse = ";"),
          max_abs_effect_change = if (length(delta)) max(delta) else NA_real_, feature_of_max_change = if (length(delta)) names(delta)[which.max(delta)] else NA_character_,
          claim_label = "descriptive", stringsAsFactors = FALSE)
      }
    }
    infl_df <- do.call(rbind, rows)
    infl_df$rank <- stats::ave(-ifelse(is.na(infl_df$max_abs_effect_change), -Inf, infl_df$max_abs_effect_change), infl_df$contrast_id, FUN = function(v) rank(v, ties.method = "min"))
    infl_df <- infl_df[order(infl_df$contrast_id, infl_df$rank), , drop = FALSE]
    write_tsv(infl_df, "influence.tsv", "post_de_influence", "InfluenceResult")
    stem <- "influence_max_effect_change"
    write_tsv(infl_df[, c("unit_id", "contrast_id", "max_abs_effect_change", "rank")], file.path("figure_sources", paste0(stem, ".tsv")), paste0(stem, "_source"), "FigureSource")
    add_fig(.pm_devices(out, stem, unlist(p$figure_formats), 7, 4.5, function() {
      graphics::par(mar = c(6, 4.5, 3, 1)); top <- infl_df[order(-infl_df$max_abs_effect_change), ][seq_len(min(30L, nrow(infl_df))), ]
      graphics::barplot(top$max_abs_effect_change, names.arg = paste(top$unit_id, top$contrast_id), las = 2, cex.names = 0.6, col = "grey40",
                        ylab = "max |effect change| (log2)", main = "Leave-one-unit-out influence (descriptive)") }))
  }

  # FR-143 per-feature robustness summary (descriptive; never re-adjusts primary P values)
  summary_rows <- list()
  for (c in primary_contrasts) {
    cid <- c$contrast_id
    pr <- primary_rows[primary_rows$contrast_id == cid & primary_rows$eligibility == "tested", , drop = FALSE]
    pmem <- .pd_members(primary_rows, cid, criteria[[1]])
    analyses <- list()
    for (mid in names(fitted)) if (cid %in% vapply(sens_tables[[mid]]$contrasts, function(x) x$contrast_id, "")) analyses[[length(analyses) + 1L]] <- list(kind = "sensitivity_model", members = .pd_members(fitted[[mid]], cid, criteria[[1]]), effects = .pd_effects(fitted[[mid]], cid))
    for (im in infl_members) if (identical(im$contrast_id, cid) && identical(im$state, "refit")) analyses[[length(analyses) + 1L]] <- list(kind = "influence", members = im$members, effects = im$effects)
    n_sens <- sum(vapply(analyses, function(a) identical(a$kind, "sensitivity_model"), logical(1))); n_infl <- sum(vapply(analyses, function(a) identical(a$kind, "influence"), logical(1)))
    for (i in seq_len(nrow(pr))) {
      f <- pr$feature_id[i]
      meets <- vapply(analyses, function(a) f %in% a$members, logical(1))
      tested <- vapply(analyses, function(a) f %in% names(a$effects), logical(1))
      same_sign <- vapply(analyses, function(a) f %in% names(a$effects) && sign(a$effects[[f]]) == sign(pr$effect[i]), logical(1))
      kinds <- vapply(analyses, function(a) a$kind, "")
      mfreq <- if (!is.null(freq_df)) freq_df$selection_frequency[freq_df$contrast_id == cid & freq_df$feature_id == f] else numeric()
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(contrast_id = cid, feature_id = f, primary_member = f %in% pmem, primary_effect = pr$effect[i],
        n_analyses = length(analyses), n_meets_criterion = sum(meets), robustness_fraction = if (length(analyses)) mean(meets) else NA_real_,
        sensitivity_model_fraction = if (n_sens) mean(meets[kinds == "sensitivity_model"]) else NA_real_,
        influence_retained_fraction = if (n_infl) mean(meets[kinds == "influence"]) else NA_real_,
        matched_n_selection_frequency = if (length(mfreq)) mfreq else NA_real_,
        sign_stability = if (any(tested)) sum(same_sign) / sum(tested) else NA_real_,
        influence_flag = (f %in% pmem) && n_infl > 0 && any(!meets[kinds == "influence"]),
        criterion = sprintf("%s < %s", criteria[[1]]$criterion, format(criteria[[1]]$threshold)), claim_label = "descriptive", stringsAsFactors = FALSE)
    }
  }
  robust <- .pd_rows(summary_rows, c("contrast_id", "feature_id", "claim_label"))
  write_tsv(robust, "robustness_summary.tsv", "post_de_robustness_summary", "RobustnessSummary")

  refusals_df <- .pd_refusal_frame(refusals)
  write_tsv(refusals_df, "refusals.tsv", "post_de_sensitivity_refusals", "PostDeEligibility")
  for (table in list(comparison, robust)) .pd_check_claims(table)
  write_json(list(module = "sensitivity", state = "COMPLETED", reason_code = NULL, claim_label = "descriptive", eligibility = p$eligibility, refusals = refusals_df,
                  primary_engine = list(engine = "limma", trend = settings$trend, robust = settings$robust, consensus_correlation = settings$consensus_correlation),
                  criteria = criteria, warning = PD_INTERACTION_WARNING, input_hashes = .pd_input_hashes(request), figures = figures,
                  rule = "sensitivity models are named secondary models re-planned through the R04 grammar with exact rank/estimability, SM04 coverage re-applied, primary engine settings, own families; resampling and influence act on biological units"),
             "eligibility.json", "post_de_sensitivity_eligibility", "PostDeEligibility")
  list(outputs = outputs, warnings = warnings, message = sprintf("post-DE sensitivity: %d model(s) fitted", length(fitted)))
})

.pd_rows <- function(rows, columns) {
  if (is.data.frame(rows)) return(rows)
  if (!length(rows)) return(stats::setNames(data.frame(matrix(character(), 0, length(columns)), stringsAsFactors = FALSE), columns))
  if (is.data.frame(rows[[1]])) return(do.call(rbind, rows))
  do.call(rbind, lapply(rows, function(r) as.data.frame(lapply(r, function(v) if (is.null(v)) NA else if (length(v) > 1L) paste(unlist(v), collapse = ";") else unlist(v)), stringsAsFactors = FALSE)))
}
