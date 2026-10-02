# Prespecified LFQ normalization policies (packet R03, SM04, FR-021/FR-022).
#
# preserve leaves values untouched.  median centering subtracts, for each
# observation, the difference between its observed median and the median of
# those medians.  reference centering does the same using only the declared
# reference proteins and fails when an observation lacks coverage.
# Quantile normalization is a named complete-case sensitivity only.

.pc_obs_for_centering <- function(values, observed) {
  usable <- !is.na(values)
  if (!is.null(observed)) usable <- usable & observed
  usable
}

normalize_preserve <- function(values) {
  list(values = values, factors = data.frame(observation_id = colnames(values), method = "preserve", center = NA_real_, factor = 0, stringsAsFactors = FALSE))
}

normalize_median <- function(values, observed = NULL) {
  usable <- .pc_obs_for_centering(values, observed)
  medians <- vapply(seq_len(ncol(values)), function(j) { v <- values[usable[, j], j]; if (!length(v)) NA_real_ else stats::median(v) }, numeric(1))
  if (any(is.na(medians))) stop(sprintf("E_NORMALIZATION_COVERAGE: observation(s) %s have no observed values for median centering", paste(colnames(values)[is.na(medians)], collapse = ", ")), call. = FALSE)
  overall <- stats::median(medians)
  factors <- medians - overall
  out <- sweep(values, 2L, factors, "-")
  list(values = out, factors = data.frame(observation_id = colnames(values), method = "median", center = medians, factor = factors, stringsAsFactors = FALSE), reference_median = overall)
}

normalize_reference <- function(values, observed = NULL, reference_features, minimum_fraction = 1) {
  reference_features <- unlist(reference_features)
  if (!length(reference_features)) stop("E_REFERENCE_COVERAGE: normalization=reference requires declared reference_features", call. = FALSE)
  unknown <- setdiff(reference_features, rownames(values))
  if (length(unknown)) stop(sprintf("E_REFERENCE_COVERAGE: reference features not in the matrix: %s", paste(unknown, collapse = ", ")), call. = FALSE)
  usable <- .pc_obs_for_centering(values, observed)[reference_features, , drop = FALSE]
  fraction <- colSums(usable) / length(reference_features)
  short <- names(fraction)[fraction < minimum_fraction | colSums(usable) == 0]
  if (length(short)) stop(sprintf("E_REFERENCE_COVERAGE: observation(s) %s cover fewer than %.3g of the declared reference proteins; all-protein centering is not substituted", paste(short, collapse = ", "), minimum_fraction), call. = FALSE)
  ref <- values[reference_features, , drop = FALSE]
  medians <- vapply(seq_len(ncol(ref)), function(j) stats::median(ref[usable[, j], j]), numeric(1))
  overall <- stats::median(medians)
  factors <- medians - overall
  list(values = sweep(values, 2L, factors, "-"),
       factors = data.frame(observation_id = colnames(values), method = "reference", center = medians, factor = factors, reference_fraction = as.numeric(fraction), stringsAsFactors = FALSE),
       reference_median = overall)
}

# Named sensitivity only: complete-case rows, limma::normalizeBetweenArrays(method="quantile").
normalize_quantile_sensitivity <- function(values) {
  if (!requireNamespace("limma", quietly = TRUE)) stop("E_DEPENDENCY_UNAVAILABLE: limma is required for the quantile normalization sensitivity", call. = FALSE)
  complete <- stats::complete.cases(values)
  if (sum(complete) < 2L) stop("E_SENSITIVITY_UNIVERSE: fewer than two complete-case features for quantile normalization", call. = FALSE)
  out <- limma::normalizeBetweenArrays(values[complete, , drop = FALSE], method = "quantile")
  list(values = out, lost_features = rownames(values)[!complete])
}
