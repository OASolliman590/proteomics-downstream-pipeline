# Post-DE calibration (R11 extension recorded with R14f; SM35, SM36, SM40).  Runs the production R functions:
#   nested CV (pd_bm_nested) on pure-noise data: the pooled out-of-fold AUC must be centred at 0.5;
#   the whole-procedure permutation test (pd_bm_permute_labels + pd_bm_nested): P must be uniform under the null;
#   a deliberately leaky reference (selection and scaling on all data, then CV): inflated AUC on the same noise;
#   the PPI connectivity null (pd_connectivity_null) on random measured sets: uniform P; planted dense clusters: P at the floor.
# Seeds are frozen here before any result is read.
# Usage: Rscript --vanilla tests/scientific/calibration_post_de.R <n_auc> <n_perm> <B> <n_ppi> <cores> <out_dir>
args <- commandArgs(TRUE)
n_auc <- as.integer(args[1]); n_perm <- as.integer(args[2]); B <- as.integer(args[3]); n_ppi <- as.integer(args[4]); cores <- as.integer(args[5]); out <- args[6]
ns <- asNamespace("proteomicsCore"); f <- function(name) get(name, envir = ns)
SEEDS <- list(auc = 600000L, perm = 700000L, leaky = 800000L, ppi = 900000L)
settings <- function(ids) list(positive = "P", negative = "N", levels = c("P", "N"), imputation = "none", selection = "top_k_auc", k = 5L, alpha = 0.5,
                               classifier = "penalized_logistic", lambda = 0.1, cost = NULL, threshold_rule = "youden_train", scheme = "repeated_stratified_kfold",
                               outer_k = 5L, repeats = 2L, inner_k = 3L, observation_ids = ids, fixed_features = NULL)
noise <- function(seed, n = 15L, p = 200L) {
  set.seed(seed, kind = "L'Ecuyer-CMRG"); X <- matrix(stats::rnorm(2L * n * p), 2L * n, p, dimnames = list(sprintf("u%02d", seq_len(2L * n)), sprintf("F%03d", seq_len(p))))
  list(X = X, y = rep(c(1L, 0L), each = n), units = rownames(X))
}
lapply_cores <- function(x, fun) if (cores > 1L && .Platform$OS.type == "unix") parallel::mclapply(x, fun, mc.cores = cores) else lapply(x, fun)
t0 <- Sys.time()
auc <- unlist(lapply_cores(seq_len(n_auc), function(i) { d <- noise(SEEDS$auc + i); f("pd_bm_nested")(d$X, d$y, d$units, settings(d$units), f(".bm_stream")(SEEDS$auc + i, 0L), record = FALSE)$pooled }))
perm <- unlist(lapply_cores(seq_len(n_perm), function(i) { d <- noise(SEEDS$perm + i); s <- settings(d$units)
  obs <- f("pd_bm_nested")(d$X, d$y, d$units, s, f(".bm_stream")(SEEDS$perm + i, 0L), record = FALSE)$pooled
  null <- vapply(seq_len(B), function(b) { st <- f(".bm_stream")(SEEDS$perm + i, b); f(".bm_use_stream")(st); yb <- f("pd_bm_permute_labels")(d$y, d$units)
    f("pd_bm_nested")(d$X, yb, d$units, s, .Random.seed, record = FALSE)$pooled }, numeric(1))
  (sum(null >= obs - 1e-12) + 1) / (B + 1) }))
leaky <- unlist(lapply_cores(seq_len(n_auc), function(i) { d <- noise(SEEDS$leaky + i); Z <- scale(d$X)
  a <- f("pd_auc_columns")(Z, d$y); top <- names(sort(-abs(a - 0.5)))[1:5]
  set.seed(SEEDS$leaky + i, kind = "L'Ecuyer-CMRG"); fold <- sample(rep(1:5, length.out = nrow(Z))); sc <- numeric(nrow(Z))
  for (k in 1:5) { tr <- fold != k; fit <- glmnet::glmnet(Z[tr, top], d$y[tr], family = "binomial", alpha = 0, lambda = 0.1); sc[!tr] <- stats::predict(fit, Z[!tr, top, drop = FALSE], s = 0.1, type = "response") }
  f("pd_auc")(sc, d$y) }))
# PPI connectivity: one synthetic measured network; random sets (null) and planted dense clusters
set.seed(SEEDS$ppi, kind = "L'Ecuyer-CMRG")
universe <- sprintf("G%03d", 1:300)
edges <- unique(t(apply(cbind(sample(universe, 1500, TRUE), sample(universe, 1500, TRUE)), 1, sort))); edges <- edges[edges[, 1] != edges[, 2], ]
cluster <- universe[1:10]; dense <- t(utils::combn(cluster, 2))
edges <- unique(rbind(edges, dense)); edges <- data.frame(a = edges[, 1], b = edges[, 2], stringsAsFactors = FALSE)
degree <- stats::setNames(vapply(universe, function(g) sum(edges$a == g | edges$b == g), numeric(1)), universe)
ppi_null <- unlist(lapply_cores(seq_len(n_ppi), function(i) { set.seed(SEEDS$ppi + i, kind = "L'Ecuyer-CMRG"); s <- sample(universe, 10)   # random sets from the whole measured universe
  r <- f("pd_connectivity_null")(s, universe, edges, degree, 199L, SEEDS$ppi + i); (sum(r$null >= r$observed) + 1) / 200 }))
r <- f("pd_connectivity_null")(cluster, universe, edges, degree, 199L, SEEDS$ppi); cluster_p <- (sum(r$null >= r$observed) + 1) / 200
dir.create(out, recursive = TRUE, showWarnings = FALSE)
w <- function(df, name) utils::write.table(df, file.path(out, name), sep = "\t", quote = FALSE, row.names = FALSE, eol = "\n")
w(data.frame(dataset = seq_along(auc), seed = SEEDS$auc + seq_along(auc), pooled_oof_auc = auc), "nested_auc_null.tsv")
w(data.frame(dataset = seq_along(perm), seed = SEEDS$perm + seq_along(perm), permutation_p = perm, B = B), "permutation_p_null.tsv")
w(data.frame(dataset = seq_along(leaky), seed = SEEDS$leaky + seq_along(leaky), leaky_auc = leaky), "leaky_reference.tsv")
w(data.frame(set = seq_along(ppi_null), seed = SEEDS$ppi + seq_along(ppi_null), connectivity_p = ppi_null), "connectivity_p_null.tsv")
jsonlite::write_json(list(seeds = SEEDS, n_auc = n_auc, n_perm = n_perm, B = B, n_ppi = n_ppi, planted_cluster_p = cluster_p, seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")),
                          design = "15 vs 15 units, 200 pure-noise features; nested CV 2 x 5-fold outer, 3-fold inner, top-5 AUC selection, ridge glmnet lambda 0.1; leaky reference selects and scales on all data",
                          network = "300 measured genes, 1500 random edges plus a planted complete cluster of 10 genes; random sets of 10 genes drawn uniformly from the measured universe; degree-preserving null of 199 draws"),
                     file.path(out, "post_de_calibration_raw.json"), auto_unbox = TRUE, digits = NA)
