# General family-adjustment API (packet R05, SM12, V046).
#
# A family is resolved from the frozen plan before results are seen.  BH or
# BY is applied once to the finite P values of every eligible planned
# endpoint in that family; planned rows that are excluded, nonestimable or
# numerically failed keep NA and are counted, never coerced to P=1.

adjust_pvalues <- function(p, method = c("BH", "BY")) {
  method <- match.arg(method)
  n <- length(p)
  if (!n) return(numeric())
  o <- order(p, decreasing = TRUE)
  q <- p[o] * n / rev(seq_len(n))
  if (method == "BY") q <- q * sum(1 / seq_len(n))
  q <- pmin(1, cummin(q))
  q[order(o)]
}

adjust_family <- function(rows, family) {
  members <- vapply(family$members, function(m) paste(m$model_id, m$contrast_id, sep = "\r"), "")
  in_family <- paste(rows$model_id, rows$contrast_id, sep = "\r") %in% members & rows$hypothesis_type == family$hypothesis_type
  eligible <- in_family & rows$eligibility == "tested" & is.finite(rows$p_value)
  rows$family_id[in_family] <- family$family_id
  rows$q_value[in_family] <- NA_real_
  rows$q_value[eligible] <- adjust_pvalues(rows$p_value[eligible], family$adjustment)
  n_failure <- sum(in_family & rows$eligibility == "numerical_failure")
  keys <- sort(paste(rows$model_id, rows$contrast_id, rows$feature_id, sep = "|")[in_family])
  summary <- data.frame(family_id = family$family_id, hypothesis_type = family$hypothesis_type, role = family$role,
                        definition_hash = sha256_text(paste(sort(members), collapse = "\n")), adjustment = family$adjustment,
                        dependence_assumption = family$dependence_assumption, n_planned = sum(in_family),
                        n_eligible = sum(in_family & rows$eligibility %in% c("tested", "numerical_failure")), n_finite = sum(eligible),
                        n_numerical_failure = n_failure, completeness = if (n_failure == 0L) "complete" else "incomplete",
                        q_cutoff = family$q_cutoff, rejection_count = sum(rows$q_value[eligible] <= family$q_cutoff),
                        universe_hash = sha256_text(paste(keys, collapse = "\n")), stringsAsFactors = FALSE)
  list(rows = rows, summary = summary)
}

sha256_text <- function(text) {
  path <- tempfile(); on.exit(unlink(path), add = TRUE)
  writeBin(charToRaw(enc2utf8(text)), path)
  sha256_file(path)
}
