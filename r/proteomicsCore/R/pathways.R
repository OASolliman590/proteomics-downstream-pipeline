# Pathway stage (packet R08: FR-071 to FR-080, SM16-SM18).  Dispatches by the
# frozen design: CAMERA only for independent designs, ROAST for paired/blocked
# designs, fgsea as an exploratory preranked gene-set null and ORA on the
# eligible measured background.  Each method/null keeps its own families.

pathway_stage <- function(request) .pc_run_stage(request, function(out) {
  p <- request$parameters
  outputs <- list(); warnings <- list()
  emit <- function(relative, id, type) outputs[[length(outputs) + 1L]] <<- .pc_output(out, relative, id, type)
  write_tsv <- function(df, relative, id, type) { .pc_write_tsv(df, file.path(out, relative)); emit(relative, id, type) }
  write_json <- function(value, relative, id, type) { .pc_write_json(value, file.path(out, relative)); emit(relative, id, type) }
  plan <- .pm_verify_inputs_response(request)
  context <- list(run_id = request$run_id, plan_hash = request$plan_hash)
  G <- .pc_matrix_from_tsv(.pc_find_input(request, "gene_matrix_finite"), "numeric")
  if (any(!is.finite(G))) stop("E_PATHWAY_MATRIX_NONFINITE: the gene matrix for matrix methods must be finite", call. = FALSE)
  membership <- .pc_read_tsv(.pc_find_input(request, "gene_set_membership"))
  summary_sets <- .pc_read_tsv(.pc_find_input(request, "gene_set_summary"))
  gene_mapping <- .pc_read_tsv(.pc_find_input(request, "gene_mapping"))
  gene_fit <- .pc_read_tsv(.pc_find_input(request, "gene_zero_null"))
  dea <- .pc_read_tsv(.pc_find_input(request, "dea_zero_null"))
  X <- .pc_matrix_from_tsv(.pc_find_input(request, paste0("design_", p$design_id)), "numeric")
  obs <- .pc_read_tsv(.pc_find_input(request, "primary_observations")); obs <- obs[match(colnames(G), obs$observation_id), , drop = FALSE]
  if (!identical(rownames(X), colnames(G))) stop("E_ID_ALIGNMENT: design rows and gene matrix columns differ", call. = FALSE)
  methods <- unlist(p$methods); blocking <- p$blocking_mode
  if ("camera" %in% methods && !identical(blocking, "none")) stop("E_CAMERA_BLOCKED_DESIGN: CAMERA is not available for paired/subject-blocked designs; use ROAST", call. = FALSE)
  if ("roast" %in% methods && identical(p$execution_profile, "production") && p$nrot < 9999) stop("E_ROTATION_RESOLUTION: production ROAST needs nrot >= 9999", call. = FALSE)
  rows <- list(); fgsea_rows <- list(); ora_rows <- list(); fg_warn <- character()
  block <- if (identical(blocking, "duplicate_correlation")) obs[[p$subject_column]] else NULL
  for (rid in unlist(p$gene_set_resource_ids)) {
    eligible <- summary_sets$set_id[summary_sets$resource_id == rid & summary_sets$eligible == "true"]
    mem <- membership[membership$resource_id == rid & membership$in_universe == "true" & membership$set_id %in% eligible, , drop = FALSE]
    sets <- split(mem$gene_id, mem$set_id)
    if (!length(sets)) next
    index <- lapply(sets, function(g) match(g, rownames(G)))
    for (contrast in p$contrasts) {
      w <- as.numeric(unlist(contrast$weights))
      base <- function(df, htype, null) { df$resource_id <- rid; df$contrast_id <- contrast$contrast_id; df$hypothesis_type <- htype; df$null_type <- null
        df$model_id <- p$gene_model_id; df$engine <- "limma"; df$role <- p$gene_model_role; df$run_id <- context$run_id; df$plan_hash <- context$plan_hash; df }
      if ("camera" %in% methods) rows[[length(rows) + 1L]] <- base(run_camera(G, index, X, w), "competitive_enrichment", "competitive_gene_set_with_estimated_inter_gene_correlation")
      if ("roast" %in% methods) {
        r <- run_roast(G, index, X, w, as.integer(p$nrot), as.integer(request$rng$seed), block, p$consensus_correlation)
        d <- base(data.frame(set_id = r$set_id, n_genes = r$n_genes, direction = r$direction, p_value = r$p_directional, native_q_value = r$native_q_directional, prop_up = r$prop_up, prop_down = r$prop_down, stringsAsFactors = FALSE), "self_contained_directional", "rotation_self_contained")
        m <- base(data.frame(set_id = r$set_id, n_genes = r$n_genes, direction = "mixed", p_value = r$p_mixed, native_q_value = r$native_q_mixed, prop_up = r$prop_up, prop_down = r$prop_down, stringsAsFactors = FALSE), "self_contained_mixed", "rotation_self_contained")
        d$nrot <- m$nrot <- p$nrot; d$resolution_label <- m$resolution_label <- if (identical(p$execution_profile, "production")) "production" else "smoke_test_only"
        rows[[length(rows) + 1L]] <- d; rows[[length(rows) + 1L]] <- m
      }
      if ("fgsea" %in% methods) {
        reps <- gene_mapping[gene_mapping$representative_state == "representative", , drop = FALSE]
        prot <- dea[dea$contrast_id == contrast$contrast_id & dea$model_id == p$primary_model_id & dea$eligibility == "tested", , drop = FALSE]
        stat <- suppressWarnings(as.numeric(prot$statistic[match(reps$feature_id, prot$feature_id)]))
        if (identical(p$rank, "signed_logp_legacy")) { pv <- suppressWarnings(as.numeric(prot$p_value[match(reps$feature_id, prot$feature_id)])); stat <- sign(stat) * -log10(pv) }
        ranks <- stats::setNames(stat, reps$gene_id); ranks <- ranks[is.finite(ranks)]
        fg <- run_fgsea(ranks, lapply(split(membership$gene_id[membership$resource_id == rid & membership$set_id %in% eligible & !is.na(membership$gene_id)],
                                            membership$set_id[membership$resource_id == rid & membership$set_id %in% eligible & !is.na(membership$gene_id)]), unique),
                        if (is.null(p$eps)) 0 else p$eps, as.integer(request$rng$seed))
        fg_warn <- c(fg_warn, fg$warnings)
        fgsea_rows[[length(fgsea_rows) + 1L]] <- base(fg$table, "preranked_gene_set", if (identical(p$rank, "signed_logp_legacy")) "preranked_gene_set_legacy_signed_logp_sensitivity" else "preranked_gene_set_exploratory")
      }
      if ("ora" %in% methods) {
        prot <- dea[dea$contrast_id == contrast$contrast_id & dea$model_id == p$primary_model_id, , drop = FALSE]
        tested <- prot$feature_id[prot$eligibility == "tested"]
        universe <- ora_universe(gene_mapping, tested)
        crit <- suppressWarnings(as.numeric(if (identical(p$ora_foreground, "family_q")) prot$q_value else prot$p_value))
        effect <- suppressWarnings(as.numeric(prot$effect))
        reps <- gene_mapping[gene_mapping$representative_state == "representative", , drop = FALSE]
        gene_of <- stats::setNames(reps$gene_id, reps$feature_id)
        allsets <- split(membership$gene_id[membership$resource_id == rid & !is.na(membership$gene_id)], membership$set_id[membership$resource_id == rid & !is.na(membership$gene_id)])
        for (direction in c("up", "down")) {
          sel <- prot$eligibility == "tested" & !is.na(crit) & crit <= p$ora_threshold & (if (direction == "up") effect > 0 else effect < 0)
          foreground <- unname(gene_of[prot$feature_id[sel]]); foreground <- foreground[!is.na(foreground)]
          o <- ora_test(universe, foreground, allsets, p$min_size, p$max_size)
          if (!is.null(o)) ora_rows[[length(ora_rows) + 1L]] <- base(o, paste0("ora_", direction), "hypergeometric_over_measured_background")
        }
      }
    }
  }
  combined <- .pc_rbind_fill(c(rows, fgsea_rows, ora_rows))
  pathway_require_finite(combined)
  combined <- pathway_family_adjust(combined, p$families)
  combined$result_type <- "GeneSetResult"; combined$schema_version <- "1.2.0"
  combined$exploratory <- combined$hypothesis_type %in% c("preranked_gene_set", "ora_up", "ora_down")
  for (h in unique(combined$hypothesis_type)) write_tsv(combined[combined$hypothesis_type == h, , drop = FALSE], paste0(h, ".tsv"), paste0("pathways_", h), "GeneSetResult")
  if (length(fgsea_rows)) {
    ov <- leading_edge_overlap(do.call(rbind, fgsea_rows))
    if (nrow(ov)) write_tsv(ov, "leading_edge_overlap.tsv", "pathways_leading_edge_overlap", "GeneSetSummary")
  }
  write_json(list(result_type = "GeneSetSummary", methods = methods, blocking_mode = blocking,
                  dispatch = "SM16: CAMERA only for independent designs; ROAST for paired/blocked; fgsea exploratory gene-set null; ORA over the eligible measured background",
                  gene_model = list(model_id = p$gene_model_id, role = p$gene_model_role), fgsea_warnings = I(unique(fg_warn)),
                  versions = list(limma = as.character(utils::packageVersion("limma")), fgsea = if (requireNamespace("fgsea", quietly = TRUE)) as.character(utils::packageVersion("fgsea")) else NULL),
                  families = lapply(p$families, function(f) f$id)), "pathway_result.json", "pathway_result", "GeneSetSummary")
  list(outputs = outputs, warnings = lapply(unique(fg_warn), function(m) .pc_warning(request, "W_FGSEA", m, "pathways/preranked_gene_set.tsv")),
       message = sprintf("pathways: %d result rows", nrow(combined)))
})

# V080: an empty eligible universe or all-nonfinite tests are failures, never a healthy zero-discovery family.
pathway_require_finite <- function(rows) {
  if (is.null(rows) || !nrow(rows)) stop("E_PATHWAY_NO_ELIGIBLE_SETS: no eligible gene set overlapped the universe", call. = FALSE)
  if (!any(is.finite(rows$p_value))) stop("E_PATHWAY_NO_FINITE_TESTS: every eligible pathway test failed numerically", call. = FALSE)
  invisible(TRUE)
}
