# Gene matrices and the gene-level model (packet R07, SM15, V066/V067).

build_gene_matrix <- function(mapped, values, observed) {
  reps <- mapped[mapped$representative_state == "representative", , drop = FALSE]
  reps <- reps[order(reps$gene_id, method = "radix"), , drop = FALSE]
  m <- values[reps$feature_id, , drop = FALSE]
  m[!observed[reps$feature_id, , drop = FALSE]] <- NA
  rownames(m) <- reps$gene_id
  list(matrix = m, representative = stats::setNames(reps$feature_id, reps$gene_id))
}

# Finite measured subset: genes observed in every required observation; no imputation.
finite_gene_matrix <- function(gene_matrix) {
  complete <- stats::complete.cases(gene_matrix)
  list(matrix = gene_matrix[complete, , drop = FALSE],
       loss = data.frame(gene_id = rownames(gene_matrix), in_finite_universe = complete,
                         reason = ifelse(complete, "complete", "incomplete_in_required_observations"), n_missing = rowSums(is.na(gene_matrix)), stringsAsFactors = FALSE))
}

median_aggregation_sensitivity <- function(mapped, values, observed) {
  ok <- mapped[mapped$mapping_state == "mapped", , drop = FALSE]
  genes <- sort(unique(ok$gene_id), method = "radix")
  out <- t(vapply(genes, function(g) {
    f <- ok$feature_id[ok$gene_id == g]
    v <- values[f, , drop = FALSE]; v[!observed[f, , drop = FALSE]] <- NA
    apply(v, 2L, function(x) if (all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE))
  }, numeric(ncol(values))))
  out <- matrix(out, nrow = length(genes), dimnames = list(genes, colnames(values)))
  out
}

# Gene-level limma refit on the finite gene matrix with the primary design (same settings).
gene_model_fit <- function(G, X, contrasts, trend, robust, block, correlation, representative, ci_level) {
  if (!requireNamespace("limma", quietly = TRUE)) stop("E_DEPENDENCY_UNAVAILABLE: limma is required for the gene model", call. = FALSE)
  if (nrow(G) < 2L) stop("E_PATHWAY_MATRIX_NONFINITE: fewer than two finite genes for the gene-level model", call. = FALSE)
  rows <- list()
  for (contrast in contrasts) {
    ex <- exact_contrast_fit(G, X, as.numeric(unlist(contrast$weights)), NULL, block, if (is.null(block)) NULL else as.numeric(correlation))
    eb <- limma::eBayes(ex$fit, trend = trend, robust = robust)
    se <- sqrt(eb$s2.post) * eb$stdev.unscaled[, 1]
    ci <- ci_moderated_t(eb$coefficients[, 1], se, eb$df.total, ci_level)
    rows[[length(rows) + 1L]] <- data.frame(contrast_id = contrast$contrast_id, gene_id = rownames(G), representative_feature_id = unname(representative[rownames(G)]),
      effect = eb$coefficients[, 1], effect_se = se, statistic = eb$t[, 1], statistic_type = "moderated_t", df_inference = eb$df.total, p_value = eb$p.value[, 1],
      ci_lower = ci[, "lower"], ci_upper = ci[, "upper"], mean_abundance = eb$Amean, stringsAsFactors = FALSE)
  }
  list(rows = do.call(rbind, rows), settings = list(model = "gene-level limma refit on the finite representative gene matrix", trend = trend, robust = robust,
       exactness = "contrast_as_coefficient_refit", block_correlation = correlation, n_genes = nrow(G), limma_version = as.character(utils::packageVersion("limma"))))
}
