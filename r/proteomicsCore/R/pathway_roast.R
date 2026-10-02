# ROAST/mroast self-contained rotation tests for paired or blocked designs (packet R08, SM16, V073/V074).
# midp = FALSE, set.statistic = "mean", approx.zscore = TRUE; directional and mixed endpoints are separate.
run_roast <- function(G, index, X, contrast, nrot, seed, block = NULL, correlation = NULL) {
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  args <- list(y = G, index = index, design = X, contrast = contrast, set.statistic = "mean", nrot = nrot, approx.zscore = TRUE, midp = FALSE, sort = "none")
  if (!is.null(block)) { args$block <- block; args$correlation <- correlation }
  res <- do.call(limma::mroast, args)
  data.frame(set_id = rownames(res), n_genes = res$NGenes, prop_down = res$PropDown, prop_up = res$PropUp, direction = res$Direction,
             p_directional = res$PValue, native_q_directional = res$FDR, p_mixed = res$PValue.Mixed, native_q_mixed = res$FDR.Mixed, stringsAsFactors = FALSE)
}
