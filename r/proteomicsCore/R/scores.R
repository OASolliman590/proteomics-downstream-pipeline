# Fixed independent scores and descriptive scores (packet R09, SM23, V082/V087).
#
# A score is inferential only when its feature set, signs/weights, centers,
# scales and missing-feature policy were frozen on participants disjoint
# from every tested biological unit and subject.  Otherwise the result is a
# DescriptiveScore that physically has no P/q/statistic columns.

score_eligibility <- function(manifest, tested_subjects, tested_units, score_test) {
  training <- unique(c(unlist(manifest$training_subject_ids), unlist(manifest$selection_subject_ids)))
  training_units <- unique(unlist(manifest$training_unit_ids))
  if (identical(score_test, "off")) return(list(inferential = FALSE, reason = "score_test_off"))
  if (!length(training) || isTRUE(manifest$selection_provenance_unknown)) return(list(inferential = FALSE, reason = "E_SCORE_INDEPENDENCE_UNVERIFIED"))
  overlap <- intersect(training, tested_subjects)
  unit_overlap <- intersect(training_units, tested_units)
  if (length(overlap) || length(unit_overlap)) return(list(inferential = FALSE, reason = "E_SCORE_SELECTION_OVERLAP", overlap = c(overlap, unit_overlap)))
  if (isTRUE(manifest$transform_refit_on_validation)) return(list(inferential = FALSE, reason = "E_SCORE_INDEPENDENCE_UNVERIFIED"))
  list(inferential = TRUE, reason = "independent")
}

apply_fixed_score <- function(Y, manifest) {
  features <- vapply(manifest$features, function(f) f$feature_id, "")
  weight <- vapply(manifest$features, function(f) as.numeric(f$weight), 0)
  center <- vapply(manifest$features, function(f) as.numeric(f$center), 0)
  scale <- vapply(manifest$features, function(f) as.numeric(f$scale), 0)
  if (any(!is.finite(scale) | scale <= 0)) stop("E_SCORE_MANIFEST: score scales must be finite and positive", call. = FALSE)
  policy <- if (is.null(manifest$missing_feature_policy)) "require_all" else manifest$missing_feature_policy
  present <- features %in% rownames(Y)
  values <- vapply(seq_len(ncol(Y)), function(j) {
    y <- rep(NA_real_, length(features)); y[present] <- Y[features[present], j]
    ok <- is.finite(y)
    if (identical(policy, "require_all") && !all(ok)) return(NA_real_)
    if (!any(ok)) return(NA_real_)
    contribution <- weight[ok] * (y[ok] - center[ok]) / scale[ok]
    if (identical(policy, "renormalize")) sum(contribution) * sum(abs(weight)) / sum(abs(weight[ok])) else sum(contribution)
  }, numeric(1))
  fraction <- vapply(seq_len(ncol(Y)), function(j) { y <- rep(NA_real_, length(features)); y[present] <- Y[features[present], j]; mean(is.finite(y)) }, numeric(1))
  data.frame(observation_id = colnames(Y), value = values, observed_feature_fraction = fraction, stringsAsFactors = FALSE)
}
