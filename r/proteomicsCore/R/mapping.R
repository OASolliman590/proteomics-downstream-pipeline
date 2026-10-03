# Species-aware protein-to-gene mapping and label-independent representatives
# (packet R07, SM14/SM15, V063/V065/V066) and the resources stage handler.

map_features_to_genes <- function(features, mapping, organism) {
  mapping <- mapping[as.character(mapping$taxonomy_id) == as.character(organism), , drop = FALSE]
  rows <- list()
  for (i in seq_len(nrow(features))) {
    fid <- features$feature_id[i]
    accessions <- unlist(jsonlite::fromJSON(features$accessions[i]))
    hits <- mapping[mapping$source_id %in% accessions, , drop = FALSE]
    current <- hits[hits$status == "current", , drop = FALSE]
    genes <- unique(current$gene_id)
    if (length(genes) == 1L) { state <- "mapped"; reason <- "unique_gene" }
    else if (length(genes) > 1L) { state <- "excluded"; reason <- "multi_gene_group" }
    else if (nrow(hits)) { state <- "excluded"; reason <- "retired_identifier" }
    else { state <- "excluded"; reason <- "unmapped" }
    rows[[i]] <- data.frame(feature_id = fid, accessions = features$accessions[i], gene_id = if (length(genes) == 1L) genes else NA_character_,
                            candidate_gene_ids = as.character(jsonlite::toJSON(genes)), gene_symbol = if (length(genes) == 1L) current$gene_symbol[current$gene_id == genes][1] else NA_character_,
                            mapping_state = state, mapping_reason = reason, stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

# Lexicographic key: descending overall genuine coverage, descending overall
# observed median, ascending feature_id (C-locale).  Uses no group labels.
select_representatives <- function(mapped, values, observed) {
  mapped$coverage <- NA_real_; mapped$overall_median <- NA_real_; mapped$representative_state <- ifelse(mapped$mapping_state == "mapped", "candidate", "excluded")
  mapped$rank_in_gene <- NA_integer_
  ok <- mapped$mapping_state == "mapped"
  for (i in which(ok)) {
    f <- mapped$feature_id[i]; use <- observed[f, ] & !is.na(values[f, ])
    mapped$coverage[i] <- sum(use); mapped$overall_median[i] <- if (any(use)) stats::median(values[f, use]) else NA_real_
  }
  for (g in unique(mapped$gene_id[ok])) {
    idx <- which(ok & mapped$gene_id == g)
    o <- idx[order(-mapped$coverage[idx], -ifelse(is.na(mapped$overall_median[idx]), -Inf, mapped$overall_median[idx]), mapped$feature_id[idx], method = "radix")]
    mapped$rank_in_gene[o] <- seq_along(o)
    mapped$representative_state[o[1]] <- "representative"
    if (length(o) > 1L) mapped$representative_state[o[-1]] <- "candidate_not_selected"
  }
  mapped
}

.rs_resource_tables <- function(request, resource) {
  out <- list()
  for (stem in unlist(resource$files)) out[[stem]] <- .pc_read_tsv(.pc_find_input(request, paste0("resfile_", resource$resource_id, "__", stem)))
  out
}

mapping_stage <- function(request) .pc_run_stage(request, function(out) {
  p <- request$parameters
  outputs <- list(); warnings <- list()
  emit <- function(relative, id, type) outputs[[length(outputs) + 1L]] <<- .pc_output(out, relative, id, type)
  write_tsv <- function(df, relative, id, type) { .pc_write_tsv(df, file.path(out, relative)); emit(relative, id, type) }
  write_json <- function(value, relative, id, type) { .pc_write_json(value, file.path(out, relative)); emit(relative, id, type) }
  plan <- .pm_verify_inputs_response(request)
  for (resource in p$resources) {
    manifest <- jsonlite::fromJSON(.pc_find_input(request, paste0("resource_", resource$resource_id)), simplifyVector = FALSE)
    for (item in manifest$files) {
      stem <- sub("\\.[^.]*$", "", item$name)
      path <- .pc_find_input(request, paste0("resfile_", resource$resource_id, "__", stem))
      .pc_verify_hash(path, item$sha256, item$name)
    }
  }
  values <- .pc_matrix_from_tsv(.pc_find_input(request, "primary_matrix"), "numeric")
  observed <- .pc_matrix_from_tsv(.pc_find_input(request, "primary_observed_mask"), "logical")[rownames(values), colnames(values), drop = FALSE]
  features <- .pc_read_tsv(.pc_find_input(request, "primary_features"))
  obs <- .pc_read_tsv(.pc_find_input(request, "primary_observations")); obs <- obs[match(colnames(values), obs$observation_id), , drop = FALSE]
  resources <- list(); for (r in p$resources) resources[[r$resource_id]] <- r
  summary <- list(result_type = "ResourceValidation", resources = lapply(p$resources, function(r) list(resource_id = r$resource_id, kind = r$kind, collection_label = r$collection_label, projection = r$projection)))
  if (!is.null(p$mapping_resource_id)) {
    mapping <- .rs_resource_tables(request, resources[[p$mapping_resource_id]])[[1]]
    mapped <- map_features_to_genes(features, mapping, p$organism_taxonomy_id)
    mapped <- select_representatives(mapped, values, observed)
    write_tsv(mapped, "gene_mapping.tsv", "gene_mapping", "GeneMapping")
    gm <- build_gene_matrix(mapped, values, observed)
    write_tsv(.pc_matrix_df(gm$matrix), "gene_matrix.tsv", "gene_matrix", "GeneMatrix")
    finite <- finite_gene_matrix(gm$matrix)
    write_tsv(.pc_matrix_df(finite$matrix), "gene_matrix_finite.tsv", "gene_matrix_finite", "GeneMatrix")
    write_tsv(finite$loss, "gene_universe_loss.tsv", "gene_universe_loss", "GeneMapping")
    if (identical(p$gene_aggregation_sensitivity, "median")) {
      med <- median_aggregation_sensitivity(mapped, values, observed)
      write_tsv(.pc_matrix_df(med), "gene_matrix_median_sensitivity.tsv", "gene_matrix_median_sensitivity", "GeneMatrixSensitivity")
    }
    X <- .pc_matrix_from_tsv(.pc_find_input(request, paste0("design_", p$model$design_id)), "numeric")
    block <- if (identical(p$model$blocking_mode, "duplicate_correlation")) obs[[p$model$subject_column]] else NULL
    gene_fit <- gene_model_fit(finite$matrix, X, p$contrasts, isTRUE(p$model$trend), isTRUE(p$model$robust), block, p$model$consensus_correlation, gm$representative, p$ci_level)
    write_tsv(gene_fit$rows, "gene_zero_null.tsv", "gene_zero_null", "GeneModelFit")
    write_json(gene_fit$settings, "gene_model_settings.json", "gene_model_settings", "GeneModelFit")
    universe <- rownames(finite$matrix)
    sets <- list(); membership <- list()
    for (rid in unlist(p$gene_set_resource_ids)) {
      res <- resources[[rid]]
      table <- .rs_resource_tables(request, res)[[1]]
      if (identical(res$collection_label, "ortholog_projected")) {
        orth_id <- res$projection$orthology_resource_id
        if (is.null(orth_id) || is.null(resources[[orth_id]])) stop(sprintf("E_RESOURCE_TAXONOMY: projected collection %s needs its declared orthology resource", rid), call. = FALSE)
        table <- project_gene_sets(table, .rs_resource_tables(request, resources[[orth_id]])[[1]], p$organism_taxonomy_id)
      } else table$source_gene_id <- table$gene_id
      table$projection_status <- if (is.null(table$projection_status)) "native" else table$projection_status
      el <- gene_set_eligibility(table, universe, p$min_size, p$max_size)
      el$summary$resource_id <- rid; el$summary$collection_label <- res$collection_label
      el$membership$resource_id <- rid; el$membership$collection_label <- res$collection_label
      sets[[length(sets) + 1L]] <- el$summary; membership[[length(membership) + 1L]] <- el$membership
    }
    if (length(sets)) {
      write_tsv(do.call(rbind, sets), "gene_set_summary.tsv", "gene_set_summary", "GeneSetEligibility")
      write_tsv(do.call(rbind, membership), "gene_set_membership.tsv", "gene_set_membership", "GeneSetEligibility")
    }
    summary$mapping <- list(n_features = nrow(mapped), mapped = sum(mapped$mapping_state == "mapped"), reasons = as.list(table(mapped$mapping_reason)),
                            n_genes = nrow(gm$matrix), n_finite_genes = nrow(finite$matrix), representative_rule = "coverage desc, overall median desc, feature_id asc (label-independent)")
  }
  write_json(summary, "mapping_summary.json", "mapping_summary", "ResourceValidation")
  list(outputs = outputs, warnings = warnings, message = "resources verified and gene matrix built")
})
