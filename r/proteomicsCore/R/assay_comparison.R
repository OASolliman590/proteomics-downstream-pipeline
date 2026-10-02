# Method sensitivity comparison (packet R06, V058).  Matches effects by stable
# feature_id, discloses unmatched universes and never changes the primary model.

compare_with_primary <- function(rows, primary, primary_model_id) {
  base <- primary[primary$model_id == primary_model_id & primary$hypothesis_type == "protein_zero_null", , drop = FALSE]
  out <- list()
  for (key in unique(paste(rows$model_id, rows$contrast_id, sep = "\r"))) {
    r <- rows[paste(rows$model_id, rows$contrast_id, sep = "\r") == key, , drop = FALSE]
    b <- base[base$contrast_id == r$contrast_id[1], , drop = FALSE]
    alt_tested <- r$feature_id[r$eligibility == "tested"]; pri_tested <- b$feature_id[b$eligibility == "tested"]
    ids <- sort(union(alt_tested, pri_tested), method = "radix")
    out[[length(out) + 1L]] <- data.frame(model_id = r$model_id[1], engine = r$engine[1], primary_model_id = primary_model_id, contrast_id = r$contrast_id[1], feature_id = ids,
      in_alternative_universe = ids %in% alt_tested, in_primary_universe = ids %in% pri_tested,
      alternative_effect = as.numeric(r$effect[match(ids, r$feature_id)]), primary_effect = suppressWarnings(as.numeric(b$effect[match(ids, b$feature_id)])),
      alternative_p_value = as.numeric(r$p_value[match(ids, r$feature_id)]), primary_p_value = suppressWarnings(as.numeric(b$p_value[match(ids, b$feature_id)])),
      role = "sensitivity_comparison_not_primary", stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}
