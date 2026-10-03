# Leading-edge redundancy and pathway family adjustment (packet R08, SM18, V077/V079).
leading_edge_overlap <- function(fg) {
  if (is.null(fg) || nrow(fg) < 2L) return(data.frame())
  le <- lapply(fg$leading_edge, function(x) unlist(jsonlite::fromJSON(x)))
  names(le) <- paste(fg$contrast_id, fg$set_id, sep = "|")
  pairs <- utils::combn(names(le), 2L, simplify = FALSE)
  do.call(rbind, lapply(pairs, function(p) { a <- le[[p[1]]]; b <- le[[p[2]]]; inter <- length(intersect(a, b)); uni <- length(union(a, b))
    data.frame(set_a = p[1], set_b = p[2], intersection = inter, jaccard = if (uni) inter / uni else NA_real_, note = "redundancy display only; original tests unchanged", stringsAsFactors = FALSE) }))
}

pathway_family_adjust <- function(rows, families) {
  rows$family_id <- NA_character_; rows$q_value <- NA_real_
  for (f in families) {
    sel <- rows$hypothesis_type == f$hypothesis_type & rows$contrast_id %in% unlist(f$contrast_ids) & rows$resource_id %in% unlist(f$collection_ids)
    rows$family_id[sel] <- f$id
    ok <- sel & is.finite(rows$p_value)
    if (any(ok)) rows$q_value[ok] <- adjust_pvalues(rows$p_value[ok], if (identical(f$adjustment, "BY")) "BY" else "BH")
  }
  rows
}
