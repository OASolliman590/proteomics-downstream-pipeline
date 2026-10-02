# Rank, aliasing and featurewise estimability (packet R04, SM07, V033/V036).

design_rank <- function(X, tol = 1e-7) {
  q <- qr(X, tol = tol)
  p <- ncol(X); rank <- q$rank
  aliased <- if (rank < p) colnames(X)[q$pivot[(rank + 1L):p]] else character()
  relations <- list()
  if (length(aliased)) {
    kept <- colnames(X)[q$pivot[seq_len(rank)]]
    for (name in aliased) {
      coef <- qr.coef(qr(X[, kept, drop = FALSE]), X[, name])
      coef <- coef[abs(coef) > 1e-9]
      relations[[name]] <- as.list(round(coef, 12))
    }
  }
  sv <- svd(X)$d
  list(rank = rank, full_rank = rank == p, aliased = aliased, alias_relations = relations,
       condition_number = if (rank == p && min(sv) > 0) max(sv) / min(sv) else NA_real_)
}

contrast_estimable <- function(X, contrast, tol = 1e-8) {
  if (!nrow(X)) return(FALSE)
  s <- svd(X)
  keep <- s$d > max(s$d) * 1e-10
  if (!any(keep)) return(FALSE)
  V <- s$v[, keep, drop = FALSE]
  residual <- contrast - V %*% (t(V) %*% contrast)
  max(abs(residual)) <= tol * max(1, max(abs(contrast)))
}

# Featurewise n/rank/df/estimability for every model x contrast, on observed rows.
featurewise_estimability <- function(values, observed, X, contrasts, coverage, model_id, groups, tol = 1e-8) {
  rows <- list()
  usable <- !is.na(values) & (if (is.null(observed)) TRUE else observed)
  for (contrast in contrasts) {
    cov <- coverage[coverage$model_id == model_id & coverage$contrast_id == contrast$contrast_id, , drop = FALSE]
    cov_by_feature <- stats::setNames(seq_len(nrow(cov)), cov$feature_id)
    required <- unlist(contrast$required_groups)
    for (i in seq_len(nrow(values))) {
      feature <- rownames(values)[i]
      use <- usable[i, ]
      Xi <- X[use, , drop = FALSE]
      rank <- if (any(use)) qr(Xi, tol = 1e-7)$rank else 0L
      n <- sum(use); df <- n - rank
      estimable <- any(use) && contrast_estimable(Xi, unlist(contrast$weights), tol)
      n_req <- vapply(required, function(g) sum(use & groups == g), integer(1))
      prior_state <- if (length(cov_by_feature) && !is.na(cov_by_feature[feature])) cov$eligibility[cov_by_feature[feature]] else "eligible"
      prior_reason <- if (length(cov_by_feature) && !is.na(cov_by_feature[feature])) cov$reason[cov_by_feature[feature]] else "eligible"
      if (prior_state != "eligible") { eligibility <- prior_state; reason <- prior_reason }
      else if (!estimable) { eligibility <- "nonestimable"; reason <- "contrast_not_estimable_on_observed_rows" }
      else if (df < 1L) { eligibility <- "nonestimable"; reason <- "no_residual_df" }
      else { eligibility <- "eligible"; reason <- "eligible" }
      rows[[length(rows) + 1L]] <- data.frame(model_id = model_id, contrast_id = contrast$contrast_id, feature_id = feature, n_obs = n, rank = rank,
        df_residual = df, estimable = estimable, eligibility = eligibility, reason = reason,
        n_obs_by_required_group = as.character(jsonlite::toJSON(as.list(stats::setNames(as.integer(n_req), required)), auto_unbox = TRUE)), stringsAsFactors = FALSE)
    }
  }
  if (length(rows)) do.call(rbind, rows) else data.frame()
}

# Express a contrast as weights on group cell means when it is a pure group contrast.
implied_group_weights <- function(X, contrast, term_map, group_levels) {
  M <- matrix(0, nrow = length(group_levels), ncol = ncol(X), dimnames = list(group_levels, colnames(X)))
  for (j in seq_along(term_map)) {
    term <- term_map[[j]]
    if (identical(term$term, "(Intercept)")) M[, j] <- 1
    if (identical(term$term, "group")) M[term$level, j] <- 1
  }
  w <- qr.solve(t(M), contrast)
  residual <- max(abs(t(M) %*% w - contrast))
  list(weights = stats::setNames(as.numeric(w), group_levels), representable = residual <= 1e-8 * max(1, max(abs(contrast))))
}
