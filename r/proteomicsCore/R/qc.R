# QC diagnostics (packet R03, SM06, V026/V027/V030).
#
# PCA uses a separate display-only matrix: features with at least one value,
# missing cells filled with the feature's observed median, then centered.
# That display artifact never feeds inference.  Pairwise correlations use
# only cells observed in both observations and report the shared n.
# Watchlist flags never remove an observation.

qc_distributions <- function(values) {
  do.call(rbind, lapply(seq_len(ncol(values)), function(j) {
    v <- values[!is.na(values[, j]), j]
    if (!length(v)) return(data.frame(observation_id = colnames(values)[j], n = 0L, min = NA_real_, q25 = NA_real_, median = NA_real_, q75 = NA_real_, max = NA_real_, mean = NA_real_, sd = NA_real_))
    q <- stats::quantile(v, c(0, .25, .5, .75, 1), names = FALSE, type = 7)
    finite_or_na <- function(x) if (is.finite(x)) x else NA_real_
    data.frame(observation_id = colnames(values)[j], n = length(v), min = q[1], q25 = q[2], median = q[3], q75 = q[4], max = q[5],
               mean = finite_or_na(mean(v)), sd = if (length(v) > 1) finite_or_na(stats::sd(v)) else NA_real_)
  }))
}

qc_correlations <- function(values, usable = !is.na(values), minimum_shared = 3L) {
  ids <- colnames(values); rows <- list()
  for (a in seq_along(ids)) for (b in seq_along(ids)) if (a < b) {
    shared <- usable[, a] & usable[, b]
    n <- sum(shared)
    r <- NA_real_; reason <- "computed"
    if (n < minimum_shared) reason <- "too_few_shared_observed_cells"
    else {
      sa <- stats::sd(values[shared, a]); sb <- stats::sd(values[shared, b])
      if (!is.finite(sa) || !is.finite(sb)) reason <- "nonfinite_variance"
      else if (sa == 0 || sb == 0) reason <- "constant_on_shared_cells"
      else { r <- stats::cor(values[shared, a], values[shared, b]); if (!is.finite(r)) { r <- NA_real_; reason <- "nonfinite_correlation" } }
    }
    rows[[length(rows) + 1L]] <- data.frame(observation_a = ids[a], observation_b = ids[b], n_shared = n, pearson_r = r, state = reason, stringsAsFactors = FALSE)
  }
  if (length(rows)) do.call(rbind, rows) else data.frame(observation_a = character(), observation_b = character(), n_shared = integer(), pearson_r = numeric(), state = character())
}

pca_display_matrix <- function(values) {
  keep <- rowSums(!is.na(values)) > 0
  display <- values[keep, , drop = FALSE]
  if (nrow(display)) for (i in seq_len(nrow(display))) { miss <- is.na(display[i, ]); if (any(miss)) display[i, miss] <- stats::median(display[i, !miss]) }
  list(filled = display, centered = display - rowMeans(display), dropped_all_missing = rownames(values)[!keep])
}

qc_pca <- function(values) {
  display <- pca_display_matrix(values)
  centered <- display$centered
  if (nrow(centered) == 0L) return(list(state = "INAPPLICABLE", reason_code = "E_QC_NO_FEATURES", message = "no retained features; PCA not computed", display = display))
  if (ncol(centered) < 2L) return(list(state = "INAPPLICABLE", reason_code = "E_QC_TOO_FEW_OBSERVATIONS", message = "fewer than two observations; PCA not computed", display = display))
  total <- sum(centered^2)
  if (!is.finite(total)) return(list(state = "INAPPLICABLE", reason_code = "E_QC_NUMERIC_OVERFLOW", message = "display matrix variance overflows double precision; PCA not computed", display = display))
  if (total <= .Machine$double.eps * max(1, sum(display$filled^2))) return(list(state = "INAPPLICABLE", reason_code = "E_QC_CONSTANT", message = "all retained features are constant; no variance to decompose", display = display))
  decomposition <- svd(t(centered))
  d <- decomposition$d
  k <- sum(d > max(d) * 1e-10)
  scores <- decomposition$u[, seq_len(k), drop = FALSE] %*% diag(d[seq_len(k)], k, k)
  colnames(scores) <- paste0("PC", seq_len(k)); rownames(scores) <- colnames(centered)
  variance <- data.frame(component = colnames(scores), eigenvalue = d[seq_len(k)]^2 / (ncol(centered) - 1), singular_value = d[seq_len(k)],
                         variance_explained = d[seq_len(k)]^2 / sum(d^2), stringsAsFactors = FALSE)
  list(state = "COMPLETED", reason_code = NULL, scores = scores, variance = variance, display = display)
}

qc_watchlist <- function(pca, correlations, ids, threshold = 3.5) {
  flags <- data.frame(observation_id = ids, pca_distance = NA_real_, pca_robust_z = NA_real_, median_correlation = NA_real_, flagged = FALSE, reasons = "", action = "none_flag_only", stringsAsFactors = FALSE)
  robust_z <- function(x) { m <- stats::median(x); s <- stats::mad(x); if (!is.finite(s) || s == 0) s <- stats::sd(x); if (!is.finite(s) || s == 0) return(rep(0, length(x))); (x - m) / s }
  if (identical(pca$state, "COMPLETED")) {
    k <- min(2L, ncol(pca$scores))
    centre <- apply(pca$scores[, seq_len(k), drop = FALSE], 2L, stats::median)
    distance <- sqrt(rowSums(sweep(pca$scores[, seq_len(k), drop = FALSE], 2L, centre)^2))
    flags$pca_distance <- distance[ids]; flags$pca_robust_z <- robust_z(distance)[ids]
  }
  if (nrow(correlations)) {
    med <- vapply(ids, function(id) { r <- correlations$pearson_r[(correlations$observation_a == id | correlations$observation_b == id) & !is.na(correlations$pearson_r)]; if (length(r)) stats::median(r) else NA_real_ }, numeric(1))
    flags$median_correlation <- med
  }
  for (i in seq_along(ids)) {
    reasons <- character()
    if (!is.na(flags$pca_robust_z[i]) && flags$pca_robust_z[i] > threshold) reasons <- c(reasons, sprintf("pca_robust_z>%.2f", threshold))
    flags$flagged[i] <- length(reasons) > 0L; flags$reasons[i] <- paste(reasons, collapse = ";")
  }
  flags
}
