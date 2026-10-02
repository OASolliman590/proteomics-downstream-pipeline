# Whole-subject influence diagnostics (packet R05, SM13, V049).
#
# Omits every record of one subject (all its visits; technical injections
# were already aggregated) and refits the same fixed model.  The output is a
# descriptive effect-change diagnostic, not held-out validation.

subject_influence <- function(Y, X, contrasts, subjects, groups, weights = NULL, block = NULL, correlation = NULL, reference_effects) {
  rows <- list()
  for (subject in unique(subjects)) {
    keep <- subjects != subject
    Xs <- X[keep, , drop = FALSE]
    nonzero <- colSums(abs(Xs)) > 0
    # omitting a reference subject aliases the remaining subject columns; drop aliased nuisance columns
    q <- qr(Xs[, nonzero, drop = FALSE])
    if (q$rank < sum(nonzero)) {
      kept <- which(nonzero)
      aliased <- kept[q$pivot[(q$rank + 1L):length(kept)]]
      nonzero[aliased] <- FALSE
    }
    Xs <- Xs[, nonzero, drop = FALSE]
    Ys <- Y[, keep, drop = FALSE]
    ws <- if (is.null(weights)) NULL else if (is.matrix(weights)) weights[, keep, drop = FALSE] else weights[keep]
    bs <- if (is.null(block)) NULL else block[keep]
    for (contrast in contrasts) {
      w <- as.numeric(unlist(contrast$weights))
      required <- unlist(contrast$required_groups)
      relevant <- sum(!keep & groups %in% required)
      units_left <- vapply(required, function(g) sum(keep & groups == g), integer(1))
      state <- "refit"; reason <- NA_character_
      if (any(w[!nonzero] != 0)) { state <- "ineligible_refit"; reason <- "contrast_uses_omitted_or_aliased_column" }
      else if (!contrast_estimable(Xs, w[nonzero])) { state <- "ineligible_refit"; reason <- "contrast_not_estimable_after_omission" }
      else if (any(units_left < 2L)) { state <- "ineligible_refit"; reason <- "fewer_than_two_units_in_required_group" }
      else if (qr(Xs)$rank < ncol(Xs)) { state <- "ineligible_refit"; reason <- "design_rank_deficient_after_omission" }
      changes <- rep(NA_real_, nrow(Y))
      if (state == "refit") {
        fit <- exact_contrast_fit(Ys, Xs, w[nonzero], ws, bs, correlation)$fit
        changes <- fit$coefficients[, 1] - reference_effects[[contrast$contrast_id]]
      }
      finite <- is.finite(changes)
      rows[[length(rows) + 1L]] <- data.frame(subject_id = subject, contrast_id = contrast$contrast_id, n_omitted_observations = sum(!keep),
        n_omitted_in_required_groups = relevant, state = state, reason = reason, n_features_compared = sum(finite),
        median_abs_effect_change = if (any(finite)) stats::median(abs(changes[finite])) else NA_real_,
        max_abs_effect_change = if (any(finite)) max(abs(changes[finite])) else NA_real_,
        interpretation = "descriptive whole-subject influence; not held-out validation", stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}
