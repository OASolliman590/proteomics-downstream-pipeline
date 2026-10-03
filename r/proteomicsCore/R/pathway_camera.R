# CAMERA competitive test for independent designs (packet R08, SM16, V072).
run_camera <- function(G, index, X, contrast) {
  res <- limma::camera(G, index = index, design = X, contrast = contrast, inter.gene.cor = NA, sort = FALSE)
  data.frame(set_id = rownames(res), n_genes = res$NGenes, correlation = res$Correlation, direction = res$Direction, p_value = res$PValue,
             native_q_value = res$FDR, stringsAsFactors = FALSE)
}
