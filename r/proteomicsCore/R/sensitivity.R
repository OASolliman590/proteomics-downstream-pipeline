# Explicit imputation sensitivities (packet R03, SM05, V028).
#
# Every algorithm writes a separate sensitivity matrix plus an imputed-cell
# mask; the primary matrix and masks are never modified.  No algorithm reads
# group labels when choosing donors or fill parameters.

.pc_sensitivity_params <- function(spec, allowed) {
  present <- intersect(setdiff(names(spec), c("id", "method", "model_id")), c("shift_sd", "scale_sd", "k", "rowmax", "colmax", "maxp"))
  extra <- setdiff(present, allowed)
  if (length(extra)) stop(sprintf("E_SENSITIVITY_PARAMETERS: method %s does not accept parameters %s", spec$method, paste(extra, collapse = ", ")), call. = FALSE)
  invisible(TRUE)
}

impute_min_deterministic <- function(values, spec) {
  .pc_sensitivity_params(spec, character())
  out <- values; mask <- matrix(FALSE, nrow(values), ncol(values), dimnames = dimnames(values)); excluded <- character()
  for (i in seq_len(nrow(values))) {
    miss <- is.na(values[i, ])
    if (all(miss)) { excluded <- c(excluded, rownames(values)[i]); next }
    if (any(miss)) { out[i, miss] <- min(values[i, !miss]); mask[i, miss] <- TRUE }
  }
  keep <- !rownames(out) %in% excluded
  list(values = out[keep, , drop = FALSE], imputed = mask[keep, , drop = FALSE], excluded = excluded, parameters = list(method = "min_deterministic", rule = "feature minimum observed log2 value across all groups"))
}

impute_left_shifted_gaussian <- function(values, spec, seed) {
  .pc_sensitivity_params(spec, c("shift_sd", "scale_sd"))
  shift <- spec$shift_sd; scale <- spec$scale_sd
  if (is.null(shift) || is.null(scale) || !is.finite(shift) || !is.finite(scale) || shift <= 0 || scale <= 0) stop("E_SENSITIVITY_PARAMETERS: left_shifted_gaussian requires positive shift_sd and scale_sd", call. = FALSE)
  mu <- apply(values, 2L, function(v) mean(v, na.rm = TRUE)); sigma <- apply(values, 2L, function(v) if (sum(!is.na(v)) > 1) stats::sd(v, na.rm = TRUE) else NA_real_)
  bad <- names(sigma)[!is.finite(sigma) | sigma == 0 | !is.finite(mu)]
  if (length(bad)) stop(sprintf("E_IMPUTATION_DISTRIBUTION: undefined or zero standard deviation in observation(s) %s", paste(bad, collapse = ", ")), call. = FALSE)
  out <- values; mask <- is.na(values)
  set.seed(as.integer(seed), kind = "L'Ecuyer-CMRG")
  feature_order <- order(rownames(values)); observation_order <- order(colnames(values))
  for (i in feature_order) for (j in observation_order) if (mask[i, j]) out[i, j] <- stats::rnorm(1L, mean = mu[j] - shift * sigma[j], sd = scale * sigma[j])
  list(values = out, imputed = mask, excluded = character(), parameters = list(method = "left_shifted_gaussian", shift_sd = shift, scale_sd = scale, seed = as.integer(seed), rng_kind = "L'Ecuyer-CMRG", draw_order = "feature_id then observation_id"))
}

impute_knn <- function(values, spec, seed) {
  .pc_sensitivity_params(spec, c("k", "rowmax", "colmax", "maxp"))
  if (!requireNamespace("impute", quietly = TRUE)) stop("E_DEPENDENCY_UNAVAILABLE: Bioconductor impute is required for the knn sensitivity", call. = FALSE)
  if (any(rowSums(!is.na(values)) == 0)) stop("E_IMPUTATION_ALL_MISSING: knn rejects all-missing feature rows", call. = FALSE)
  if (any(colSums(!is.na(values)) == 0)) stop("E_IMPUTATION_ALL_MISSING: knn rejects all-missing observation columns", call. = FALSE)
  k <- if (is.null(spec$k)) 10L else as.integer(spec$k); rowmax <- if (is.null(spec$rowmax)) 0.5 else spec$rowmax; colmax <- if (is.null(spec$colmax)) 0.8 else spec$colmax; maxp <- if (is.null(spec$maxp)) 1500L else as.integer(spec$maxp)
  row_fraction <- rowMeans(is.na(values))
  fallback <- rownames(values)[row_fraction > rowmax]
  result <- utils::capture.output(fit <- impute::impute.knn(values, k = k, rowmax = rowmax, colmax = colmax, maxp = maxp, rng.seed = as.integer(seed)))
  list(values = fit$data, imputed = is.na(values), excluded = character(),
       parameters = list(method = "knn", package = "impute", package_version = as.character(utils::packageVersion("impute")), k = k, rowmax = rowmax, colmax = colmax, maxp = maxp, rng_seed = as.integer(seed),
                         rows_using_package_mean_fallback = fallback, fallback_note = "impute.knn replaces rows with more than rowmax missing by column means"))
}

impute_complete_case <- function(values, spec) {
  .pc_sensitivity_params(spec, character())
  keep <- stats::complete.cases(values)
  list(values = values[keep, , drop = FALSE], imputed = matrix(FALSE, sum(keep), ncol(values), dimnames = list(rownames(values)[keep], colnames(values))),
       excluded = rownames(values)[!keep], parameters = list(method = "complete_case", rule = "features observed in every retained observation; no fill"))
}

run_sensitivity <- function(values, spec, seed) {
  switch(spec$method,
         min_deterministic = impute_min_deterministic(values, spec),
         left_shifted_gaussian = impute_left_shifted_gaussian(values, spec, seed),
         knn = impute_knn(values, spec, seed),
         complete_case = impute_complete_case(values, spec),
         stop(sprintf("E_SENSITIVITY_METHOD: unknown sensitivity method %s", spec$method), call. = FALSE))
}
