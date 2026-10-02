# Multivariate PERMANOVA and PERMDISP (packet R13, SM26-SM30, FR-121 to FR-130).
#
# Distance: declared metric (default Euclidean) on per-feature scaled
# (default z-scored) genuinely observed log2 abundance, computed within each
# test's observations; features must be observed in every analysed
# observation and non-constant.  vegan::adonis2 with marginal SS; the group
# term comes only from unrestricted (or subject-blocked) permutations and
# covariate/interaction terms only from permutations restricted within the
# primary group.  PERMDISP accompanies every group/covariate row.

PERMANOVA_INTERPRETATION_RULE <- paste(
  "A significant PERMANOVA with a non-significant PERMDISP indicates a location (centroid) shift;",
  "when both are significant, location and/or dispersion differences are possible;",
  "a significant PERMDISP alone indicates a dispersion difference without established location shift;",
  "non-significance is not evidence of no difference.")

.pm_require <- function() {
  for (pkg in c("vegan", "permute")) if (!requireNamespace(pkg, quietly = TRUE)) stop(sprintf("E_ENGINE_NOT_AVAILABLE: R package %s is required for PERMANOVA", pkg), call. = FALSE)
  invisible(TRUE)
}

.pm_how <- function(nperm, blocks = NULL) {
  if (is.null(blocks)) permute::how(nperm = nperm) else permute::how(nperm = nperm, blocks = factor(blocks))
}

# Usable features of a test: observed in every analysed observation and non-constant.
pm_usable_features <- function(Y, features, keep) {
  sub <- Y[features, keep, drop = FALSE]
  complete <- stats::complete.cases(sub)
  spread <- apply(sub, 1L, function(v) if (all(is.finite(v))) stats::sd(v) else NA_real_)
  constant <- complete & (is.na(spread) | spread <= sqrt(.Machine$double.eps) * pmax(1, abs(rowMeans(sub))))
  data.frame(feature_id = features, used = complete & !constant,
             reason = ifelse(!complete, "incomplete_in_test_observations", ifelse(constant, "constant_within_test", "used")), stringsAsFactors = FALSE)
}

pm_scaled <- function(Y, features, keep, scaling) {
  X <- t(Y[features, keep, drop = FALSE])
  if (identical(scaling, "zscore")) X <- scale(X, center = TRUE, scale = TRUE)
  else if (identical(scaling, "center")) X <- scale(X, center = TRUE, scale = FALSE)
  X <- unclass(X); attr(X, "scaled:center") <- NULL; attr(X, "scaled:scale") <- NULL
  X
}

pm_distance <- function(X, metric) {
  if (identical(metric, "euclidean")) stats::dist(X, method = "euclidean") else vegan::vegdist(X, method = metric)
}

.pm_term <- function(a, term) {
  row <- a[term, , drop = FALSE]
  list(df = row$Df, ss = row$SumOfSqs, r2 = row$R2, f = row$F, p = row[["Pr(>F)"]],
       residual_df = a["Residual", "Df"], residual_ss = a["Residual", "SumOfSqs"], total_ss = a["Total", "SumOfSqs"])
}

pm_display <- function(p, nperm) {
  floor <- 1 / (nperm + 1)
  if (is.na(p)) return(NA_character_)
  if (p <= floor * (1 + 1e-9)) paste0("< ", formatC(floor, format = "g", digits = 3)) else formatC(p, format = "g", digits = 3)
}

# by is always explicit: vegan >= 2.7 defaults to a single "Model" row when by = NULL.
pm_adonis <- function(D, data, formula_rhs, permutations, seed, by = "terms") {
  set.seed(seed)
  f <- stats::as.formula(paste("D ~", formula_rhs))
  environment(f) <- environment()
  vegan::adonis2(f, data = data, permutations = permutations, by = by)
}

pm_permdisp <- function(D, groups, permutations, seed) {
  bd <- vegan::betadisper(D, factor(groups))
  set.seed(seed)
  pt <- vegan::permutest(bd, permutations = permutations)
  list(betadisper = bd, f = pt$tab[1, "F"], p = pt$tab[1, "Pr(>F)"], nperm = NROW(pt$perm))
}

pm_interpret <- function(p, p_disp, alpha) {
  if (is.na(p_disp)) return("permdisp_not_applicable")
  if (p < alpha && p_disp >= alpha) "location_shift"
  else if (p < alpha && p_disp < alpha) "location_and_or_dispersion"
  else if (p >= alpha && p_disp < alpha) "dispersion_difference_location_not_established"
  else "no_separation_detected"
}

pm_univariate_r2 <- function(X, groups) {
  g <- factor(groups)
  vapply(seq_len(ncol(X)), function(j) { a <- stats::anova(stats::lm(X[, j] ~ g)); a[["Sum Sq"]][1] / sum(a[["Sum Sq"]]) }, numeric(1))
}

pm_check_interaction <- function(groups, covariate, column) {
  cells <- table(factor(groups), factor(covariate))
  empty <- which(cells == 0, arr.ind = TRUE)
  if (nrow(empty)) stop(sprintf("E_PERMANOVA_INTERACTION_NONESTIMABLE: group x %s has empty cell(s): %s", column,
                                paste(sprintf("%s/%s", rownames(cells)[empty[, 1]], colnames(cells)[empty[, 2]]), collapse = ", ")), call. = FALSE)
  invisible(cells)
}

.pm_adjust <- function(p, method) if (identical(method, "none")) p else stats::p.adjust(p, method = method)

# One feature set: global, pairwise, covariate-adjusted, interaction, PERMDISP, contributions.
pm_feature_set <- function(Y, set, obs, settings) {
  group_col <- settings$group_column; groups_all <- obs[[group_col]]
  levels_present <- settings$group_levels[settings$group_levels %in% groups_all]
  blocks <- if (!is.null(settings$subject_column)) obs[[settings$subject_column]] else NULL
  metric <- settings$metric; scaling <- settings$scaling; nperm <- settings$permutations; seed <- settings$seed; alpha <- settings$alpha
  rows <- list(); usage <- list(); nulls <- list(); ordination <- NULL; feature_r2 <- NULL; identity <- NULL
  run_test <- function(keep, analysis, comparison) {
    use <- pm_usable_features(Y, set$features, keep)
    usage[[length(usage) + 1L]] <<- cbind(feature_set_id = set$id, test = paste(analysis, comparison, sep = ":"), use)
    feats <- use$feature_id[use$used]
    if (length(feats) < 2L) return(list(refused = TRUE, features = feats))
    X <- pm_scaled(Y, feats, keep, scaling)
    list(refused = FALSE, features = feats, X = X, D = pm_distance(X, metric))
  }
  add_row <- function(analysis, comparison, term, n, n_features, stats, actual_nperm, scheme, disp = NULL, disp_scheme = NA_character_) {
    rows[[length(rows) + 1L]] <<- data.frame(feature_set_id = set$id, analysis = analysis, comparison = comparison, term = term, n = n, n_features = n_features,
      df = stats$df, ss = stats$ss, r2 = stats$r2, pseudo_f = stats$f, residual_df = stats$residual_df, residual_ss = stats$residual_ss, total_ss = stats$total_ss,
      p_value = stats$p, p_floor = 1 / (actual_nperm + 1), p_display = pm_display(stats$p, actual_nperm), p_adjusted = NA_real_, p_adjusted_display = NA_character_,
      adjustment = NA_character_, nperm = actual_nperm, seed = seed, permutation_scheme = scheme,
      permdisp_f = if (is.null(disp)) NA_real_ else disp$f, permdisp_p = if (is.null(disp)) NA_real_ else disp$p,
      permdisp_p_display = if (is.null(disp)) NA_character_ else pm_display(disp$p, disp$nperm), permdisp_scheme = disp_scheme,
      interpretation = NA_character_, metric = metric, scaling = scaling, stringsAsFactors = FALSE)
  }
  record_null <- function(a, term, test_id) {
    ps <- vegan::permustats(a)
    j <- match(term, colnames(ps$permutations)); if (is.na(j)) j <- 1L
    nulls[[length(nulls) + 1L]] <<- data.frame(feature_set_id = set$id, test_id = test_id, term = term, permutation = seq_len(nrow(ps$permutations)),
                                                pseudo_f = unname(ps$permutations[, j]), observed_pseudo_f = unname(ps$statistic[j]), stringsAsFactors = FALSE, row.names = NULL)
    nrow(ps$permutations)
  }
  group_scheme <- if (is.null(blocks)) "free" else paste0("within:", settings$subject_column)

  # global
  keep <- groups_all %in% levels_present
  g <- run_test(keep, "global", "all_groups")
  if (g$refused) return(list(refused = TRUE, reason = "fewer_than_two_usable_features", usage = do.call(rbind, usage)))
  data <- data.frame(group = factor(groups_all[keep], levels = levels_present))
  perm <- .pm_how(nperm, if (is.null(blocks)) NULL else blocks[keep])
  a <- pm_adonis(g$D, data, "group", perm, seed)
  actual <- record_null(a, "group", "global:all_groups")
  disp <- pm_permdisp(g$D, data$group, .pm_how(nperm, if (is.null(blocks)) NULL else blocks[keep]), seed)
  add_row("global", "all_groups", "group", sum(keep), length(g$features), .pm_term(a, "group"), actual, group_scheme, disp, group_scheme)
  bd <- disp$betadisper
  eig <- bd$eig[bd$eig > 0]
  coords <- bd$vectors[, seq_len(min(2L, ncol(bd$vectors))), drop = FALSE]
  if (ncol(coords) < 2L) coords <- cbind(coords, 0)
  centroid <- bd$centroids[, seq_len(min(2L, ncol(bd$centroids))), drop = FALSE]
  if (ncol(centroid) < 2L) centroid <- cbind(centroid, 0)
  grp <- as.character(data$group)
  ordination <- data.frame(feature_set_id = set$id, observation_id = obs$observation_id[keep], group = grp, axis1 = coords[, 1], axis2 = coords[, 2],
                           centroid_axis1 = centroid[grp, 1], centroid_axis2 = centroid[grp, 2], distance_to_centroid = bd$distances, stringsAsFactors = FALSE)
  axes <- data.frame(feature_set_id = set$id, axis = c("axis1", "axis2"), percent_positive_eigenvalue = 100 * c(eig[1], if (length(eig) > 1) eig[2] else 0) / sum(eig), stringsAsFactors = FALSE)
  centroids <- data.frame(feature_set_id = set$id, group = rownames(centroid), axis1 = centroid[, 1], axis2 = centroid[, 2], stringsAsFactors = FALSE)
  # per-feature contribution and identity (global groups, same scaling)
  r2_uni <- pm_univariate_r2(g$X, data$group)
  means <- sapply(levels_present, function(l) colMeans(t(Y[g$features, keep, drop = FALSE])[data$group == l, , drop = FALSE]))
  means <- matrix(means, nrow = length(g$features), dimnames = list(g$features, levels_present))
  feature_r2 <- data.frame(feature_set_id = set$id, feature_id = g$features, univariate_group_r2 = r2_uni, stringsAsFactors = FALSE)
  for (l in levels_present) feature_r2[[paste0("mean_", l)]] <- means[, l]
  global_r2 <- .pm_term(a, "group")$r2
  identity_applicable <- identical(metric, "euclidean") && identical(scaling, "zscore")
  identity <- list(feature_set_id = set$id, multivariate_r2 = global_r2, mean_univariate_r2 = mean(r2_uni),
                   abs_difference = abs(mean(r2_uni) - global_r2), tolerance = 1e-8,
                   state = if (!identity_applicable) "not_applicable" else if (abs(mean(r2_uni) - global_r2) <= 1e-8) "pass" else "fail")
  if (identical(identity$state, "fail")) stop(sprintf("E_PERMANOVA_IDENTITY: feature set %s: mean univariate R2 %.12g differs from multivariate R2 %.12g", set$id, mean(r2_uni), global_r2), call. = FALSE)

  # pairwise
  if (isTRUE(settings$pairwise) && length(levels_present) > 2L) {
    pairs <- utils::combn(levels_present, 2L, simplify = FALSE)
    first <- length(rows) + 1L
    for (pair in pairs) {
      keep_pair <- groups_all %in% pair
      t <- run_test(keep_pair, "pairwise", paste(pair, collapse = "_vs_"))
      if (t$refused) next
      dpair <- data.frame(group = factor(groups_all[keep_pair], levels = pair))
      blk <- if (is.null(blocks)) NULL else blocks[keep_pair]
      ap <- pm_adonis(t$D, dpair, "group", .pm_how(nperm, blk), seed)
      act <- record_null(ap, "group", paste0("pairwise:", paste(pair, collapse = "_vs_")))
      dp <- pm_permdisp(t$D, dpair$group, .pm_how(nperm, blk), seed)
      add_row("pairwise", paste(pair[2], "vs", pair[1]), "group", sum(keep_pair), length(t$features), .pm_term(ap, "group"), act, group_scheme, dp, group_scheme)
    }
    idx <- seq.int(first, length(rows))
    if (length(rows) >= first) {
      adjusted <- .pm_adjust(vapply(rows[idx], function(r) r$p_value, numeric(1)), settings$adjustment)
      for (k in seq_along(idx)) { rows[[idx[k]]]$p_adjusted <- adjusted[k]; rows[[idx[k]]]$adjustment <- settings$adjustment
        rows[[idx[k]]]$p_adjusted_display <- pm_display(adjusted[k], rows[[idx[k]]]$nperm) }
    }
  }

  # covariates: marginal SS; group term unrestricted, covariate terms blocked within group
  covariates <- settings$covariates
  if (length(covariates)) {
    if (!is.null(blocks)) stop("E_PERMANOVA_DESIGN_UNSUPPORTED: covariate PERMANOVA is not supported for subject-blocked designs", call. = FALSE)
    cdata <- data.frame(group = factor(groups_all[keep], levels = levels_present))
    names_cov <- character()
    for (cv in covariates) {
      column <- cv$column
      if (!column %in% names(obs)) stop(sprintf("E_PERMANOVA_COVARIATE: covariate column %s is absent", column), call. = FALSE)
      value <- obs[[column]][keep]
      if (identical(cv$type, "continuous")) {
        num <- suppressWarnings(as.numeric(value))
        if (any(!is.finite(num))) stop(sprintf("E_PERMANOVA_COVARIATE: continuous covariate %s has missing or nonfinite values", column), call. = FALSE)
        cdata[[paste0("cov_", column)]] <- num
      } else {
        if (any(is.na(value))) stop(sprintf("E_PERMANOVA_COVARIATE: categorical covariate %s has missing values", column), call. = FALSE)
        cdata[[paste0("cov_", column)]] <- factor(value)
      }
      names_cov <- c(names_cov, paste0("cov_", column))
    }
    rhs <- paste(c("group", names_cov), collapse = " + ")
    free <- pm_adonis(g$D, cdata, rhs, .pm_how(nperm), seed, by = "margin")
    act <- record_null(free, "group", "group_adjusted:all_groups")
    add_row("group_adjusted", paste("all groups adjusted for", paste(vapply(covariates, function(cv) cv$column, ""), collapse = ", ")), "group",
            sum(keep), length(g$features), .pm_term(free, "group"), act, "free", disp, "free")
    blocked <- pm_adonis(g$D, cdata, rhs, .pm_how(nperm, cdata$group), seed, by = "margin")
    for (k in seq_along(covariates)) {
      term <- names_cov[k]; cv <- covariates[[k]]
      act_b <- record_null(blocked, term, paste0("covariate:", cv$column))
      cdisp <- if (identical(cv$type, "categorical")) pm_permdisp(g$D, cdata[[term]], .pm_how(nperm, cdata$group), seed) else NULL
      add_row("covariate", paste(cv$column, "adjusted for group"), cv$column, sum(keep), length(g$features), .pm_term(blocked, term), act_b,
              paste0("within:", group_col), cdisp, if (is.null(cdisp)) NA_character_ else paste0("within:", group_col))
    }
    if (isTRUE(settings$interaction)) {
      for (k in seq_along(covariates)) {
        cv <- covariates[[k]]; term <- names_cov[k]
        if (!identical(cv$type, "categorical")) next
        pm_check_interaction(cdata$group, cdata[[term]], cv$column)
        rhs_i <- paste(c(paste0("group * ", term), setdiff(names_cov, term)), collapse = " + ")
        ai <- pm_adonis(g$D, cdata, rhs_i, .pm_how(nperm, cdata$group), seed, by = "margin")
        iterm <- paste0("group:", term)
        act_i <- record_null(ai, iterm, paste0("interaction:", cv$column))
        add_row("interaction", paste("group x", cv$column), paste0("group:", cv$column), sum(keep), length(g$features), .pm_term(ai, iterm), act_i, paste0("within:", group_col))
      }
    }
  }
  tests <- do.call(rbind, rows)
  for (i in seq_len(nrow(tests))) {
    p_use <- if (!is.na(tests$p_adjusted[i])) tests$p_adjusted[i] else tests$p_value[i]
    tests$interpretation[i] <- if (tests$analysis[i] %in% c("interaction")) "permdisp_not_applicable" else pm_interpret(p_use, tests$permdisp_p[i], alpha)
  }
  list(refused = FALSE, tests = tests, usage = do.call(rbind, usage), nulls = do.call(rbind, nulls), ordination = ordination, axes = axes, centroids = centroids,
       feature_r2 = feature_r2, identity = identity, global_r2 = global_r2, r2_universe = stats::setNames(r2_uni, g$features), features = g$features)
}

# Random equal-size feature-set null and best-possible set (SM30).
pm_selection_context <- function(set_id, set_r2, k, universe_r2, draws, seed, identity_applicable, direct_r2 = NULL) {
  set.seed(seed)
  null <- if (identity_applicable) replicate(draws, mean(sample(universe_r2, k))) else replicate(draws, direct_r2(sample(names(universe_r2), k)))
  best <- if (identity_applicable) mean(sort(universe_r2, decreasing = TRUE)[seq_len(k)]) else NA_real_
  list(null = data.frame(feature_set_id = set_id, draw = seq_len(draws), r2 = null, stringsAsFactors = FALSE),
       summary = data.frame(feature_set_id = set_id, k = k, set_r2 = set_r2, null_draws = draws, seed = seed, null_median = stats::median(null),
                            null_q95 = unname(stats::quantile(null, 0.95)), fraction_null_at_or_above_set = mean(null >= set_r2 - 1e-12),
                            best_possible_r2 = best, best_possible_note = if (identity_applicable) "mean of the top-k univariate R2 (exact under z-scored Euclidean distance)" else "not computable without the mean-R2 identity",
                            all_feature_mean_r2 = mean(universe_r2), in_sample_optimistic = TRUE, stringsAsFactors = FALSE))
}

.pm_verify_inputs <- function(request) {
  plan_path <- .pc_find_input(request, "plan")
  plan <- jsonlite::fromJSON(plan_path, simplifyVector = FALSE)
  if (is.null(request$plan_hash) || !identical(plan$plan_hash, request$plan_hash)) stop("E_PLAN_CHANGED: request plan_hash does not match the frozen plan", call. = FALSE)
  planned <- list(); for (a in plan$artifacts) planned[[a$artifact_id]] <- a
  for (item in request$inputs) {
    if (item$artifact_id %in% c("plan", "dea_zero_null")) next
    ref <- planned[[item$artifact_id]]
    if (is.null(ref) || !identical(ref$sha256, item$sha256)) stop(sprintf("E_PLAN_CHANGED: input %s differs from the frozen plan", item$artifact_id), call. = FALSE)
    if (identical(ref$result_type, "DisplayOnlyMatrix")) stop("E_DISPLAY_MATRIX_REJECTED: display-only matrices cannot feed PERMANOVA", call. = FALSE)
  }
  plan
}

pm_resolve_sets <- function(sets, features, dea_path) {
  out <- list()
  dea <- if (!is.null(dea_path)) .pc_read_tsv(dea_path) else NULL
  for (s in sets) {
    kind <- s$kind
    if (identical(kind, "all_complete")) members <- features
    else if (identical(kind, "declared_panel")) {
      members <- unlist(s$feature_ids)
      unknown <- setdiff(members, features)
      if (length(unknown)) stop(sprintf("E_PERMANOVA_FEATURE_SET: panel %s names features not in the primary matrix: %s", s$id, paste(unknown, collapse = ", ")), call. = FALSE)
    } else if (identical(kind, "dep_derived")) {
      if (is.null(dea)) stop(sprintf("E_PERMANOVA_FEATURE_SET: DEP-derived set %s needs a completed differential table", s$id), call. = FALSE)
      rows <- dea[dea$model_id == s$model_id & dea$contrast_id == s$contrast_id & dea$eligibility == "tested", , drop = FALSE]
      if (!nrow(rows)) stop(sprintf("E_PERMANOVA_FEATURE_SET: no tested rows for %s/%s", s$model_id, s$contrast_id), call. = FALSE)
      value <- suppressWarnings(as.numeric(if (identical(s$criterion, "family_q")) rows$q_value else rows$p_value))
      members <- rows$feature_id[!is.na(value) & value <= s$threshold]
    } else stop(sprintf("E_PERMANOVA_FEATURE_SET: unknown kind %s", kind), call. = FALSE)
    provenance <- if (identical(kind, "dep_derived")) "same_data" else if (identical(kind, "all_complete")) "not_selected" else s$selection_provenance
    out[[length(out) + 1L]] <- list(id = s$id, kind = kind, features = members, selection_provenance = provenance)
  }
  out
}

permanova_stage <- function(request) .pc_run_stage(request, function(out) {
  .pm_require()
  p <- request$parameters
  outputs <- list(); warnings <- list()
  emit <- function(relative, id, type) outputs[[length(outputs) + 1L]] <<- .pc_output(out, relative, id, type)
  write_tsv <- function(df, relative, id, type) { .pc_write_tsv(df, file.path(out, relative)); emit(relative, id, type) }
  write_json <- function(value, relative, id, type) { .pc_write_json(value, file.path(out, relative)); emit(relative, id, type) }
  .pm_verify_inputs(request)
  if (identical(p$execution_profile, "production") && p$permutations < 999) stop("E_PERMANOVA_RESOLUTION: production PERMANOVA needs at least 999 permutations", call. = FALSE)
  Y <- .pc_matrix_from_tsv(.pc_find_input(request, "primary_matrix"), "numeric")
  observed <- .pc_matrix_from_tsv(.pc_find_input(request, "primary_observed_mask"), "logical")
  Y[!observed[rownames(Y), colnames(Y)]] <- NA
  obs <- .pc_read_tsv(.pc_find_input(request, "primary_observations"))
  obs <- obs[match(colnames(Y), obs$observation_id), , drop = FALSE]
  dea_path <- .pc_find_input(request, "dea_zero_null", required = FALSE)
  settings <- list(group_column = p$group_column, group_levels = unlist(p$group_levels), subject_column = p$subject_column, metric = p$metric, scaling = p$scaling,
                   permutations = as.integer(p$permutations), seed = as.integer(p$seed), alpha = p$alpha, pairwise = isTRUE(p$pairwise), adjustment = p$adjustment,
                   covariates = p$covariates, interaction = isTRUE(p$interaction))
  sets <- pm_resolve_sets(p$feature_sets, rownames(Y), dea_path)
  tests <- list(); usage <- list(); nulls <- list(); ord <- list(); axes <- list(); cents <- list(); fr2 <- list(); identity <- list(); context <- list(); random <- list(); states <- list()
  universe <- NULL
  for (set in sets) {
    res <- pm_feature_set(Y, set, obs, settings)
    usage[[length(usage) + 1L]] <- res$usage
    if (isTRUE(res$refused)) {
      states[[set$id]] <- list(state = "REFUSED", reason_code = "E_PERMANOVA_FEATURE_SET", reason = res$reason, n_candidate_features = length(set$features))
      warnings[[length(warnings) + 1L]] <- .pc_warning(request, "E_PERMANOVA_FEATURE_SET", sprintf("feature set %s has fewer than two usable features and was not tested", set$id), "permanova/feature_usage.tsv")
      next
    }
    states[[set$id]] <- list(state = "COMPLETED", kind = set$kind, selection_provenance = set$selection_provenance, n_features_used = length(res$features))
    tests[[length(tests) + 1L]] <- res$tests; nulls[[length(nulls) + 1L]] <- res$nulls; ord[[length(ord) + 1L]] <- res$ordination
    axes[[length(axes) + 1L]] <- res$axes; cents[[length(cents) + 1L]] <- res$centroids; fr2[[length(fr2) + 1L]] <- res$feature_r2; identity[[set$id]] <- res$identity
    if (identical(set$kind, "all_complete") && is.null(universe)) universe <- res$r2_universe
    set$result <- res
    sets[[which(vapply(sets, function(s) identical(s$id, set$id), logical(1)))]] <- set
  }
  if (!length(tests)) stop("E_PERMANOVA_FEATURE_SET: no feature set had at least two usable features", call. = FALSE)
  # selection-circularity context for same-data / unknown-provenance / DEP-derived sets
  identity_applicable <- identical(settings$metric, "euclidean") && identical(settings$scaling, "zscore")
  if (is.null(universe)) {
    groups_all <- obs[[settings$group_column]]
    keep_all <- groups_all %in% settings$group_levels
    base <- pm_usable_features(Y, rownames(Y), keep_all)
    feats <- base$feature_id[base$used]
    if (length(feats) >= 2L) universe <- stats::setNames(pm_univariate_r2(pm_scaled(Y, feats, keep_all, settings$scaling), groups_all[keep_all]), feats)
  }
  for (set in sets) {
    if (is.null(set$result) || !set$selection_provenance %in% c("same_data", "unknown")) next
    k <- length(set$result$features)
    if (is.null(universe) || length(universe) <= k) { warnings[[length(warnings) + 1L]] <- .pc_warning(request, "W_SELECTION_CONTEXT_UNAVAILABLE", sprintf("random-set null for %s needs a larger complete universe", set$id)); next }
    direct <- function(members) {
      groups_all <- obs[[settings$group_column]]; keep <- groups_all %in% settings$group_levels
      X <- pm_scaled(Y, members, keep, settings$scaling); D <- pm_distance(X, settings$metric)
      vegan::adonis2(D ~ g, data = data.frame(g = factor(groups_all[keep])), permutations = 0, by = "terms")$R2[1]
    }
    sc <- pm_selection_context(set$id, set$result$global_r2, k, universe, as.integer(p$random_sets), settings$seed, identity_applicable, direct)
    context[[length(context) + 1L]] <- sc$summary; random[[length(random) + 1L]] <- sc$null
  }
  tests_df <- do.call(rbind, tests)
  write_tsv(tests_df, "tests.tsv", "permanova_tests", "PermanovaResult")
  write_tsv(do.call(rbind, usage), "feature_usage.tsv", "permanova_feature_usage", "MultivariateFeatureSet")
  write_tsv(do.call(rbind, nulls), "permutation_null.tsv", "permanova_permutation_null", "FigureSource")
  write_tsv(do.call(rbind, ord), "ordination.tsv", "permanova_ordination", "FigureSource")
  write_tsv(do.call(rbind, axes), "ordination_axes.tsv", "permanova_ordination_axes", "FigureSource")
  write_tsv(do.call(rbind, cents), "centroids.tsv", "permanova_centroids", "FigureSource")
  write_tsv(do.call(rbind, fr2), "feature_r2.tsv", "permanova_feature_r2", "PermanovaContribution")
  write_json(unname(identity), "identity_check.json", "permanova_identity", "PermanovaContribution")
  if (length(context)) {
    write_tsv(do.call(rbind, context), "selection_context.tsv", "permanova_selection_context", "SelectionContext")
    write_tsv(do.call(rbind, random), "random_set_null.tsv", "permanova_random_set_null", "FigureSource")
  }
  figures <- permanova_figures(out, tests_df, do.call(rbind, ord), do.call(rbind, axes), do.call(rbind, cents), do.call(rbind, nulls), do.call(rbind, fr2),
                               if (length(random)) do.call(rbind, random) else NULL, if (length(context)) do.call(rbind, context) else NULL, unlist(p$figure_formats))
  for (f in figures$files) emit(f$relative_path, f$artifact_id, "Figure")
  write_json(list(result_type = "PermanovaResult", settings = list(metric = settings$metric, scaling = settings$scaling, permutations_requested = settings$permutations,
                    seed = settings$seed, rng_kind = "L'Ecuyer-CMRG", adjustment = settings$adjustment, alpha = settings$alpha, pairwise = settings$pairwise,
                    covariates = settings$covariates, interaction = settings$interaction, group_column = settings$group_column,
                    group_permutation_scheme = if (is.null(settings$subject_column)) "free" else paste0("within:", settings$subject_column),
                    feature_universe = "genuinely observed cells; complete and non-constant within each test's observations; no imputation"),
                  versions = list(vegan = as.character(utils::packageVersion("vegan")), permute = as.character(utils::packageVersion("permute"))),
                  feature_sets = states, identity = unname(identity), figures = figures$records, interpretation_rule = PERMANOVA_INTERPRETATION_RULE,
                  limitations = c("PERMANOVA describes multivariate separation of the analysed biological units; it is not a classifier or held-out validation.",
                                  "R2 of a feature set chosen on the same data is optimistic by construction (see selection_context.tsv).",
                                  "Permutation P values are Monte Carlo estimates with resolution 1/(nperm+1).")),
             "permanova_result.json", "permanova_result", "PermanovaResult")
  list(outputs = outputs, warnings = warnings, message = sprintf("PERMANOVA: %d test row(s) over %d feature set(s)", nrow(tests_df), length(tests)))
})
