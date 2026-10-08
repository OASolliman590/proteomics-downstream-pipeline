# Post-DE protein-phenotype association (packet R14c, SM33, FR-144 to FR-148).
#
# Model method: proteins are the response and the phenotype is a design term of a re-planned design (adjusted for the
# declared covariates and, by default, group), fitted with the primary engine settings (pd_fit_model); one family per
# phenotype (and per group for within-group analyses).  Correlation method: Pearson, Spearman or partial (Pearson on
# residuals after the adjustment terms) on features observed in every analysed unit, with permutation P = (k+1)/(B+1)
# over biological units (whole subjects, or within subjects, in subject-blocked designs) and BH within the family.
# Phenotypes are never imputed; complete-case n is reported.  Pooled/within-group sign disagreement is flagged.

PD_SIMPSON_WARNING <- "Pooled and within-group directions disagree for flagged features (Simpson pattern): the pooled association may reflect group differences rather than a within-group relationship."

# Permutation vectors over units: free, whole-subject (value constant within subject) or within-subject.
pd_unit_permutations <- function(n, B, seed, subjects = NULL) {
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  if (is.null(subjects)) return(lapply(seq_len(B), function(b) sample.int(n)))
  lapply(seq_len(B), function(b) {
    idx <- seq_len(n)
    for (s in unique(subjects)) { w <- which(subjects == s); if (length(w) > 1L) idx[w] <- w[sample.int(length(w))] }
    idx
  })
}

pd_subject_permutations <- function(subjects, B, seed) {
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  u <- unique(subjects)
  lapply(seq_len(B), function(b) { perm <- stats::setNames(u[sample.int(length(u))], u); match(perm[subjects], subjects) })
}

.pd_residualize <- function(v, Z) if (is.null(Z) || !ncol(Z)) v - mean(v) else as.numeric(stats::lm.fit(Z, v)$residuals)

# Correlation of one feature with the phenotype and its permutation P (k counts |r_perm| >= |r_obs|, the observed included via +1).
# method: pearson | spearman (plain), partial or pearson_partial (residuals on Z), spearman_partial (ranks of y and x
# residualised on Z, then Pearson; D-59). Partial forms permute the residualised phenotype (D-47).
pd_correlation <- function(y, x, method, Z, perms) {
  if (identical(method, "spearman_partial")) { y <- rank(y); x <- rank(x); method <- "partial" }
  if (method %in% c("partial", "pearson_partial")) { ry <- .pd_residualize(y, Z); rx <- .pd_residualize(x, Z); r <- stats::cor(ry, rx); null <- vapply(perms, function(p) stats::cor(ry, rx[p]), numeric(1)) }
  else { r <- stats::cor(y, x, method = method); null <- vapply(perms, function(p) stats::cor(y, x[p], method = method), numeric(1)) }
  k <- sum(abs(null) >= abs(r) - 1e-12)
  list(r = r, k = k, p = (k + 1) / (length(perms) + 1))
}

.pd_within_slopes <- function(y, x, groups) {
  out <- c()
  for (g in unique(groups)) { w <- groups == g & is.finite(y); if (sum(w) >= 3L && stats::var(x[w]) > 0) out[g] <- unname(stats::coef(stats::lm(y[w] ~ x[w]))[2]) }
  out
}

post_de_association_stage <- function(request) .pc_run_stage(request, function(out) {
  p <- request$parameters
  outputs <- list(); warnings <- list(); refusals <- list(); figures <- list()
  emit <- function(relative, id, type) outputs[[length(outputs) + 1L]] <<- .pc_output(out, relative, id, type)
  write_tsv <- function(df, relative, id, type) { .pc_write_tsv(df, file.path(out, relative)); emit(relative, id, type) }
  write_json <- function(value, relative, id, type) { .pc_write_json(value, file.path(out, relative)); emit(relative, id, type) }
  add_fig <- function(x) { for (f in x$files) emit(f$relative_path, f$artifact_id, "Figure"); figures <<- c(figures, x$records) }
  plan <- .pd_verify_inputs(request)
  warnings <- c(warnings, .pd_adaptation_warnings(request))
  prim <- .pd_primary(request); obs <- prim$obs
  values <- .pc_matrix_from_tsv(.pc_find_input(request, "primary_matrix"), "numeric")
  observed <- .pc_matrix_from_tsv(.pc_find_input(request, "primary_observed_mask"), "logical")[rownames(values), colnames(values), drop = FALSE]
  context <- list(run_id = request$run_id, plan_hash = request$plan_hash)
  settings <- list(group_column = p$group_column, subject_column = p$subject_column, blocking_mode = p$blocking_mode, consensus_correlation = p$consensus_correlation,
                   trend = isTRUE(p$trend), robust = isTRUE(p$robust), ci_level = p$ci_level)
  for (s in p$eligibility$subanalyses) refusals[[length(refusals) + 1L]] <- data.frame(analysis = s$analysis, item = s$item, reason_code = s$reason_code, reason = s$reason, stringsAsFactors = FALSE)
  for (r in refusals) warnings[[length(warnings) + 1L]] <- .pc_warning(request, r$reason_code, sprintf("association %s %s refused: %s", r$analysis, r$item, r$reason), "post_de/association/refusals.tsv")
  results <- list(); summaries <- list(); simpson_any <- FALSE
  for (a in p$analyses) {
    X <- .pd_design_from(a$design); ids <- rownames(X); o <- obs[match(ids, obs$observation_id), , drop = FALSE]
    rk <- design_rank(X); if (!identical(as.integer(rk$rank), as.integer(a$rank))) stop(sprintf("E_INTEGRITY: association design %s rank differs from the plan", a$analysis_id), call. = FALSE)
    x <- as.numeric(unlist(a$values)[ids])
    family <- list(family_id = a$family_id, hypothesis_type = "protein_zero_null", role = "secondary", adjustment = "BH", q_cutoff = 0.05,
                   dependence_assumption = "BH: independence or positive regression dependence", members = list(list(model_id = a$analysis_id, contrast_id = a$analysis_id)))
    Y <- values[, ids, drop = FALSE]; Y[!observed[, ids, drop = FALSE]] <- NA
    if (identical(p$method, "model")) {
      res <- pd_fit_model(a$analysis_id, values[, ids, drop = FALSE], observed[, ids, drop = FALSE], o, X, a$contrasts, settings, list(family), context, p$coverage)
      rows <- res$rows
      tab <- data.frame(analysis_id = a$analysis_id, phenotype = a$phenotype, scope = a$scope, group = if (is.null(a$group)) NA_character_ else a$group,
        feature_id = rows$feature_id, eligibility = rows$eligibility, reason_code = rows$reason_code, effect = rows$effect, effect_scale = "log2 abundance per phenotype unit",
        ci_lower = rows$ci_lower, ci_upper = rows$ci_upper, statistic = rows$statistic, statistic_type = rows$statistic_type, p_value = rows$p_value, q_value = rows$q_value,
        family_id = rows$family_id, method = "limma (primary engine settings)", stringsAsFactors = FALSE)
    } else {
      complete <- rowSums(!is.na(Y)) == ncol(Y)
      used_method <- if (is.null(a$correlation_used)) p$correlation else a$correlation_used
      Z <- NULL
      if (used_method %in% c("partial", "pearson_partial", "spearman_partial")) { keepcol <- setdiff(colnames(X), c("(Intercept)", a$coefficient)); Z <- cbind(1, X[, keepcol, drop = FALSE]) }
      subjects <- if (!is.null(settings$subject_column)) o[[settings$subject_column]] else NULL
      constant <- !is.null(subjects) && all(tapply(x, subjects, function(v) length(unique(v))) == 1L)
      perms <- if (is.null(subjects)) pd_unit_permutations(length(ids), as.integer(p$permutations), as.integer(p$seed))
               else if (constant) pd_subject_permutations(subjects, as.integer(p$permutations), as.integer(p$seed))
               else pd_unit_permutations(length(ids), as.integer(p$permutations), as.integer(p$seed), subjects)
      scheme <- if (is.null(subjects)) "biological_units" else if (constant) paste0("whole_subjects:", settings$subject_column) else paste0("within:", settings$subject_column)
      stats <- lapply(rownames(Y), function(f) if (complete[f]) pd_correlation(Y[f, ], x, used_method, Z, perms) else NULL)
      r <- vapply(stats, function(s) if (is.null(s)) NA_real_ else s$r, numeric(1)); pv <- vapply(stats, function(s) if (is.null(s)) NA_real_ else s$p, numeric(1))
      kk <- vapply(stats, function(s) if (is.null(s)) NA_integer_ else as.integer(s$k), integer(1))
      q <- rep(NA_real_, length(pv)); q[complete] <- stats::p.adjust(pv[complete], "BH")
      tab <- data.frame(analysis_id = a$analysis_id, phenotype = a$phenotype, scope = a$scope, group = if (is.null(a$group)) NA_character_ else a$group,
        feature_id = rownames(Y), eligibility = ifelse(complete, "tested", "excluded"), reason_code = ifelse(complete, NA_character_, "incomplete_for_correlation"),
        effect = r, effect_scale = paste(used_method, "correlation"), ci_lower = NA_real_, ci_upper = NA_real_, statistic = r, statistic_type = paste0(used_method, "_r"),
        p_value = pv, q_value = q, family_id = a$family_id, method = sprintf("%s correlation (requested %s), permutation P = (k+1)/(B+1), B = %d, seed %d, scheme %s", used_method, p$correlation, as.integer(p$permutations), as.integer(p$seed), scheme),
        correlation_requested = p$correlation, correlation_used = used_method,
        permutation_k = kk, permutations = as.integer(p$permutations), permutation_scheme = scheme, stringsAsFactors = FALSE)
    }
    # review 2026-10-05 M1: report only the adjustment actually applied (the model and partial correlation use the design's group/covariate columns;
    # Pearson/Spearman use none, and the planner refuses them when an adjustment is declared)
    partial <- !is.null(a$correlation_used) && grepl("partial$", a$correlation_used) || identical(p$correlation, "partial")
    used <- if (identical(p$method, "model") || partial) c(if (isTRUE(a$adjust_for_group)) "group", unlist(a$adjust_for)) else character()
    tab$n_units <- a$n_units; tab$n_missing_phenotype <- a$n_missing; tab$adjusted_for <- if (length(used)) paste(used, collapse = ";") else "none"
    tab$simpson_flag <- FALSE; tab$within_group_slopes <- NA_character_; tab$pooled_unadjusted_slope <- NA_real_
    if (identical(a$scope, "pooled")) {
      groups <- o[[settings$group_column]]
      for (i in seq_len(nrow(tab))) {
        if (!identical(tab$eligibility[i], "tested")) next
        yi <- Y[tab$feature_id[i], ]; ok <- is.finite(yi) & is.finite(x)
        # review 2026-10-05 minor 3: Simpson's pattern compares the within-group slopes with the *unadjusted* pooled slope
        # (a group-adjusted pooled slope is a weighted average of the within-group slopes and cannot disagree with all of them)
        pooled <- if (sum(ok) >= 3L && stats::var(x[ok]) > 0) unname(stats::coef(stats::lm.fit(cbind(1, x[ok]), yi[ok]))[2]) else NA_real_
        tab$pooled_unadjusted_slope[i] <- pooled
        sl <- .pd_within_slopes(yi, x, groups)
        if (length(sl)) {
          tab$within_group_slopes[i] <- paste(sprintf("%s:%.6g", names(sl), sl), collapse = ";")
          tab$simpson_flag[i] <- length(sl) >= 2L && is.finite(pooled) && all(sign(sl) == sign(sl[1])) && sign(sl[1]) != 0 && sign(pooled) != 0 && sign(sl[1]) != sign(pooled)
        }
      }
      if (any(tab$simpson_flag)) { simpson_any <- TRUE; warnings[[length(warnings) + 1L]] <- .pc_warning(request, "W_PHENOTYPE_SIMPSON", sprintf("%s: %d feature(s) with pooled/within-group sign disagreement", a$phenotype, sum(tab$simpson_flag)), "post_de/association") }
      if (isTRUE(a$separated_by_group)) warnings[[length(warnings) + 1L]] <- .pc_warning(request, "W_PHENOTYPE_GROUP_SEPARATED", sprintf("%s: phenotype ranges do not overlap between groups; the pooled association is confounded with group differences", a$phenotype), "post_de/association")
    }
    tab$result_type <- "PhenotypeAssociationResult"
    results[[a$analysis_id]] <- tab
    summaries[[length(summaries) + 1L]] <- data.frame(analysis_id = a$analysis_id, phenotype = a$phenotype, scope = a$scope, group = if (is.null(a$group)) NA_character_ else a$group,
      family_id = a$family_id, method = p$method, n_units = a$n_units, n_missing_phenotype = a$n_missing, adjusted_for = tab$adjusted_for[1], n_tested = sum(tab$eligibility == "tested"),
      n_q_below_0_05 = sum(tab$q_value < 0.05, na.rm = TRUE), n_simpson_flag = sum(tab$simpson_flag), stringsAsFactors = FALSE)
  }
  for (ph in unique(vapply(p$analyses, function(a) a$phenotype, ""))) {
    tabs <- do.call(rbind, results[vapply(p$analyses, function(a) identical(a$phenotype, ph), logical(1))])
    write_tsv(tabs, paste0("association_", .pd_safe(ph), ".tsv"), paste0("association_", .pd_safe(ph)), "PhenotypeAssociationResult")
  }
  summary_df <- .pd_rows(summaries, c("analysis_id", "phenotype"))
  write_tsv(summary_df, "association_summary.tsv", "post_de_association_summary", "PhenotypeAssociationResult")
  # heatmap: features x analyses of the test statistic (moderated t or r), rows clustered; every plotted cell is in the source table
  all_tab <- do.call(rbind, results)
  stat <- tapply(all_tab$statistic, list(all_tab$feature_id, all_tab$analysis_id), function(v) v[1])
  stat <- stat[rowSums(is.finite(stat)) > 0, , drop = FALSE]
  if (nrow(stat)) {
    ord <- if (nrow(stat) > 2L) stats::hclust(stats::dist(ifelse(is.finite(stat), stat, 0)))$order else seq_len(nrow(stat))
    src <- expand.grid(row = seq_len(nrow(stat)), column = seq_len(ncol(stat)))
    heat <- data.frame(feature_id = rownames(stat)[src$row], analysis_id = colnames(stat)[src$column], value = stat[cbind(src$row, src$column)],
                       statistic = all_tab$statistic_type[1], row_order = match(src$row, ord), column_order = src$column, stringsAsFactors = FALSE)
    write_tsv(heat, "heatmap_source.tsv", "post_de_association_heatmap_source", "FigureSource")
    M <- stat[ord, , drop = FALSE]
    add_fig(.pm_devices(out, "association_heatmap", unlist(p$figure_formats), 6.5, max(4, 0.12 * nrow(M) + 2), function() {
      graphics::par(mar = c(8, 6, 3, 1)); lim <- max(abs(M), na.rm = TRUE); if (!is.finite(lim) || lim == 0) lim <- 1
      graphics::image(seq_len(ncol(M)), seq_len(nrow(M)), t(M), col = grDevices::hcl.colors(41, "Blue-Red 3"), zlim = c(-lim, lim), axes = FALSE, xlab = "", ylab = "",
                      main = sprintf("Phenotype association (%s)", all_tab$statistic_type[1]), cex.main = 0.9)
      graphics::axis(1, at = seq_len(ncol(M)), labels = colnames(M), las = 2, cex.axis = 0.6); graphics::axis(2, at = seq_len(nrow(M)), labels = rownames(M), las = 1, cex.axis = 0.5) }))
  }
  refusals_df <- .pd_refusal_frame(refusals)
  write_tsv(refusals_df, "refusals.tsv", "post_de_association_refusals", "PostDeEligibility")
  write_json(list(module = "association", state = "COMPLETED", reason_code = NULL, claim_label = NULL, inference = "one secondary family per phenotype (and per group for within-group analyses); BH",
                  method = p$method, correlation = if (identical(p$method, "correlation")) p$correlation else NULL, permutations = if (identical(p$method, "correlation")) p$permutations else NULL,
                  seed = p$seed, simpson_warning = if (simpson_any) PD_SIMPSON_WARNING else NULL, eligibility = p$eligibility, refusals = refusals_df, input_hashes = .pd_input_hashes(request),
                  figures = figures, missing_values = "complete case per phenotype; phenotypes are never imputed",
                  rule = "proteins are the response; the phenotype is a design term adjusted for declared covariates and, by default, group; pooled analysis refused when the phenotype is aliased with group"),
             "eligibility.json", "post_de_association_eligibility", "PostDeEligibility")
  list(outputs = outputs, warnings = warnings, message = sprintf("post-DE association: %d analysis/analyses", length(results)))
})
