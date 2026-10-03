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
# Rank and estimability depend only on the observed-row pattern, so they are computed once per
# distinct pattern; columns are built as vectors (one data frame per contrast) so the cost is
# linear in the number of features (audit 2026-10-02; outputs unchanged).
featurewise_estimability <- function(values, observed, X, contrasts, coverage, model_id, groups, tol = 1e-8) {
  rows <- list()
  usable <- !is.na(values) & (if (is.null(observed)) TRUE else observed)
  features <- rownames(values); n_feat <- nrow(values)
  pattern <- apply(usable, 1L, function(u) paste(as.integer(u), collapse = ""))
  for (contrast in contrasts) {
    cov <- coverage[coverage$model_id == model_id & coverage$contrast_id == contrast$contrast_id, , drop = FALSE]
    cov_by_feature <- stats::setNames(seq_len(nrow(cov)), cov$feature_id)
    required <- unlist(contrast$required_groups)
    w <- unlist(contrast$weights)
    rank <- integer(n_feat); estimable <- logical(n_feat)
    for (key in unique(pattern)) {
      members <- which(pattern == key); use <- usable[members[1L], ]
      Xi <- X[use, , drop = FALSE]
      rank[members] <- if (any(use)) qr(Xi, tol = 1e-7)$rank else 0L
      estimable[members] <- any(use) && contrast_estimable(Xi, w, tol)
    }
    n <- as.integer(rowSums(usable)); df <- n - rank
    n_req <- vapply(required, function(g) as.integer(rowSums(usable[, groups == g, drop = FALSE])), integer(n_feat))
    n_req <- matrix(n_req, nrow = n_feat, dimnames = list(features, required))
    idx <- if (length(cov_by_feature)) cov_by_feature[features] else rep(NA_integer_, n_feat)
    prior_state <- ifelse(is.na(idx), "eligible", cov$eligibility[idx]); prior_reason <- ifelse(is.na(idx), "eligible", cov$reason[idx])
    eligibility <- ifelse(prior_state != "eligible", prior_state, ifelse(!estimable, "nonestimable", ifelse(df < 1L, "nonestimable", "eligible")))
    reason <- ifelse(prior_state != "eligible", prior_reason, ifelse(!estimable, "contrast_not_estimable_on_observed_rows", ifelse(df < 1L, "no_residual_df", "eligible")))
    counts_json <- vapply(seq_len(n_feat), function(i) as.character(jsonlite::toJSON(as.list(stats::setNames(n_req[i, ], required)), auto_unbox = TRUE)), "")
    rows[[length(rows) + 1L]] <- data.frame(model_id = rep(model_id, n_feat), contrast_id = rep(contrast$contrast_id, n_feat), feature_id = features, n_obs = n, rank = rank,
      df_residual = df, estimable = estimable, eligibility = eligibility, reason = reason, n_obs_by_required_group = counts_json, stringsAsFactors = FALSE, row.names = NULL)
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
