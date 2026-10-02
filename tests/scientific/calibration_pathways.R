# Correlated pathway calibration (packet R11, V106).  Competitive (CAMERA) and
# self-contained (ROAST) nulls are evaluated separately under their own
# hypotheses; fgsea is not given a sample-label guarantee.
# Usage: Rscript --vanilla tests/scientific/calibration_pathways.R <n_datasets> <nrot> <out.json>
args <- commandArgs(TRUE); n_datasets <- as.integer(args[1]); nrot <- as.integer(args[2]); out <- args[3]
simulate_sets <- function(seed, n_genes = 600, n_sets = 12, set_size = 20, rho = 0.3, spike = FALSE, n_per_group = 5) {
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  g <- rep(c(0, 1), each = n_per_group); n <- length(g)
  Y <- matrix(stats::rnorm(n_genes * n), n_genes)
  sets <- split(seq_len(n_sets * set_size), rep(seq_len(n_sets), each = set_size))
  for (s in sets) { shared <- stats::rnorm(n); Y[s, ] <- sqrt(1 - rho) * Y[s, ] + sqrt(rho) * matrix(shared, length(s), n, byrow = TRUE) }
  if (spike) Y[sets[[1]], g == 1] <- Y[sets[[1]], g == 1] + 1
  rownames(Y) <- paste0("g", seq_len(n_genes)); names(sets) <- paste0("S", seq_along(sets))
  list(Y = Y, X = cbind(A = 1 - g, B = g), sets = sets)
}
res <- lapply(seq_len(n_datasets), function(i) {
  null <- simulate_sets(1000 + i); spike <- simulate_sets(5000 + i, spike = TRUE)
  cam <- limma::camera(null$Y, null$sets, design = null$X, contrast = c(-1, 1), inter.gene.cor = NA, sort = FALSE)
  set.seed(9000 + i, kind = "L'Ecuyer-CMRG")
  ro <- limma::mroast(null$Y, null$sets, design = null$X, contrast = c(-1, 1), nrot = nrot, midp = FALSE, set.statistic = "mean", approx.zscore = TRUE, sort = "none")
  cam_s <- limma::camera(spike$Y, spike$sets, design = spike$X, contrast = c(-1, 1), inter.gene.cor = NA, sort = FALSE)
  set.seed(9500 + i, kind = "L'Ecuyer-CMRG")
  ro_s <- limma::mroast(spike$Y, spike$sets, design = spike$X, contrast = c(-1, 1), nrot = nrot, midp = FALSE, set.statistic = "mean", approx.zscore = TRUE, sort = "none")
  c(camera_null_rate = mean(cam$PValue < 0.05), roast_null_rate = mean(ro$PValue < 0.05), roast_mixed_null_rate = mean(ro$PValue.Mixed < 0.05),
    camera_power = as.numeric(cam_s["S1", "PValue"] < 0.05), roast_power = as.numeric(ro_s["S1", "PValue"] < 0.05))
})
m <- do.call(rbind, res)
summary <- lapply(colnames(m), function(k) { v <- m[, k]; se <- stats::sd(v) / sqrt(length(v)); list(mean = mean(v), mc_se = se, upper95 = mean(v) + 1.645 * se) })
names(summary) <- colnames(m)
jsonlite::write_json(list(n_datasets = n_datasets, nrot = nrot, hypotheses = list(camera = "competitive (set vs other genes, estimated inter-gene correlation)", roast = "self-contained (set vs zero)"),
                          fgsea = "not calibrated as a sample-label test (exploratory gene-set randomization)", results = summary,
                          gate = list(null_rate_upper95_max = 0.075, camera_pass = summary$camera_null_rate$upper95 <= 0.075, roast_pass = summary$roast_null_rate$upper95 <= 0.075,
                                      roast_mixed_pass = summary$roast_mixed_null_rate$upper95 <= 0.075)),
                     out, auto_unbox = TRUE, pretty = TRUE, digits = NA)
