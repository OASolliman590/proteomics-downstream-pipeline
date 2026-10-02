# Gene sets: orthology projection, overlap eligibility and the ORA universe
# (packet R07, SM14/SM17, V064/V068/V069).

project_gene_sets <- function(sets, orthology, target_taxonomy) {
  orth <- orthology[as.character(orthology$target_taxonomy_id) == as.character(target_taxonomy), , drop = FALSE]
  rows <- list()
  for (i in seq_len(nrow(sets))) {
    hits <- orth[orth$source_gene_id == sets$gene_id[i], , drop = FALSE]
    if (!nrow(hits)) { rows[[length(rows) + 1L]] <- data.frame(set_id = sets$set_id[i], set_name = sets$set_name[i], source_gene_id = sets$gene_id[i], gene_id = NA_character_, projection_status = "no_ortholog", evidence = NA_character_, stringsAsFactors = FALSE); next }
    if (any(hits$ambiguous == "true") || length(unique(hits$target_gene_id)) > 1L) {
      rows[[length(rows) + 1L]] <- data.frame(set_id = sets$set_id[i], set_name = sets$set_name[i], source_gene_id = sets$gene_id[i], gene_id = NA_character_, projection_status = "ambiguous_ortholog_excluded", evidence = paste(hits$evidence, collapse = ";"), stringsAsFactors = FALSE)
      next
    }
    rows[[length(rows) + 1L]] <- data.frame(set_id = sets$set_id[i], set_name = sets$set_name[i], source_gene_id = sets$gene_id[i], gene_id = hits$target_gene_id[1], projection_status = "projected", evidence = hits$evidence[1], stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

# Size bounds are applied to set ∩ universe before any test result is seen.
gene_set_eligibility <- function(sets, universe, min_size, max_size) {
  ids <- unique(sets$set_id)
  summary <- do.call(rbind, lapply(ids, function(s) {
    rows <- sets[sets$set_id == s, , drop = FALSE]
    mapped <- unique(rows$gene_id[!is.na(rows$gene_id)])
    overlap <- intersect(mapped, universe)
    data.frame(set_id = s, set_name = rows$set_name[1], original_size = length(unique(rows$source_gene_id)), mapped_size = length(mapped),
               eligible_overlap_size = length(overlap), min_size = min_size, max_size = max_size,
               eligible = length(overlap) >= min_size && length(overlap) <= max_size,
               bounds_note = sprintf("test-only bounds %d-%d applied to set ∩ universe before testing", as.integer(min_size), as.integer(max_size)), stringsAsFactors = FALSE)
  }))
  membership <- data.frame(set_id = sets$set_id, set_name = sets$set_name, source_gene_id = sets$source_gene_id, gene_id = sets$gene_id,
                           projection_status = sets$projection_status, in_universe = !is.na(sets$gene_id) & sets$gene_id %in% universe, stringsAsFactors = FALSE)
  list(summary = summary, membership = membership)
}

# ORA universe: measured, mapped genes whose representative protein passed model eligibility.
ora_universe <- function(gene_mapping, eligible_feature_ids) {
  reps <- gene_mapping[gene_mapping$representative_state == "representative", , drop = FALSE]
  sort(unique(reps$gene_id[reps$feature_id %in% eligible_feature_ids]), method = "radix")
}
