# Detection-only exploratory endpoint (packet R03, V029).
#
# Independent two-group Fisher exact tests on genuine observed/missing counts
# form their own DetectionResult family with a narrow detection-only BH.  This
# is not the general R05 family-adjustment API and never enters a protein
# abundance family.  Paired/repeated designs are unsupported.

detection_bh <- function(p) {
  n <- length(p); if (!n) return(numeric())
  o <- order(p, decreasing = TRUE)
  q <- pmin(1, cummin(p[o] * n / rev(seq_len(n))))
  q[order(o)]
}

detection_tests <- function(observed, observations, group_column, comparisons, blocking_mode, family_id = "detection-exploratory") {
  if (!identical(blocking_mode, "none")) stop("E_DETECTION_DESIGN_UNSUPPORTED: detection-only Fisher tests require an independent two-group design", call. = FALSE)
  if (is.null(observed)) stop("E_ORIGINAL_MASK_REQUIRED: detection analysis requires a genuine original-observed mask", call. = FALSE)
  groups <- observations[[group_column]]; rows <- list()
  for (comparison in comparisons) {
    pair <- unlist(comparison$required_groups)
    if (length(pair) != 2L) stop(sprintf("E_DETECTION_DESIGN_UNSUPPORTED: contrast %s does not compare exactly two groups", comparison$contrast_id), call. = FALSE)
    a <- groups == pair[1]; b <- groups == pair[2]
    for (i in seq_len(nrow(observed))) {
      oa <- sum(observed[i, a]); ob <- sum(observed[i, b])
      table <- matrix(c(oa, sum(a) - oa, ob, sum(b) - ob), nrow = 2L)
      rows[[length(rows) + 1L]] <- data.frame(schema_version = "1.2.0", result_type = "DetectionResult", hypothesis_type = "detection", family_id = family_id,
        engine = "fisher_exact", engine_version = as.character(getRversion()), contrast_id = comparison$contrast_id, feature_id = rownames(observed)[i],
        group_a = pair[1], group_b = pair[2], observed_a = oa, missing_a = sum(a) - oa, observed_b = ob, missing_b = sum(b) - ob,
        p_value = stats::fisher.test(table, alternative = "two.sided")$p.value, stringsAsFactors = FALSE)
    }
  }
  result <- do.call(rbind, rows)
  result$q_value <- detection_bh(result$p_value)
  result$adjustment <- "BH_detection_family_only"
  result
}
