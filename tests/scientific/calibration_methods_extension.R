# R15c classifier calibration (FR-179 to FR-184; SM45; ADR 0011).  Runs the production nested procedure (pd_bm_nested,
# R14d) for every calibrated family on synthetic noise and planted-signal data, plus the leaky reference, the LOOCV
# outer scheme and the subject-blocked units.  Every count and seed comes from docs/validation/016-methods-extension-validation/
# R15c/gates-preregistration.json, which is committed before this script is run; this script does not choose counts.
#
# Usage: Rscript --vanilla tests/scientific/calibration_methods_extension.R <gates.json> <out_dir> <cores>
# Writes one TSV per scenario (rows carry the family) and calibration_raw.json (seeds, counts, R and package versions, timing).

args <- commandArgs(TRUE)
gates_path <- args[1]; out <- args[2]; cores <- as.integer(args[3])
gates <- jsonlite::fromJSON(gates_path, simplifyVector = FALSE)
ns <- asNamespace("proteomicsCore"); f <- function(name) get(name, envir = ns)
lapply_cores <- function(x, fun) {   # a failed replicate stops the run (mclapply would otherwise return try-error objects)
  res <- if (cores > 1L && .Platform$OS.type == "unix") parallel::mclapply(x, fun, mc.cores = cores) else lapply(x, fun)
  if (any(vapply(res, inherits, logical(1), "try-error")) || any(vapply(res, is.null, logical(1))))
    stop("E_CALIBRATION_REPLICATE: a calibration replicate failed; the run is not accepted", call. = FALSE)
  res
}
write_tsv <- function(df, name) utils::write.table(df, file.path(out, name), sep = "\t", quote = FALSE, row.names = FALSE, eol = "\n")
rc <- gates$replicate_counts; sd_ <- gates$seeds
families <- gates$families
classifier_of <- function(fam) fam$family

settings <- function(ids, classifier, scheme = "repeated_stratified_kfold", outer_k = 5L, repeats = 2L, fixed = NULL) list(
  positive = "P", negative = "N", levels = c("P", "N"), imputation = "none", selection = "top_k_auc", k = 5L, alpha = 0.5,
  classifier = classifier, lambda = 0.1, cost = NULL, threshold_rule = "youden_train", scheme = scheme, outer_k = outer_k,
  repeats = repeats, inner_k = 3L, observation_ids = ids, fixed_features = fixed)

noise <- function(seed, n = 15L, p = 200L, shift = 0) {
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  X <- matrix(stats::rnorm(2L * n * p), 2L * n, p, dimnames = list(sprintf("u%02d", seq_len(2L * n)), sprintf("F%03d", seq_len(p))))
  y <- rep(c(1L, 0L), each = n)
  if (shift > 0) X[y == 1L, 1:5] <- X[y == 1L, 1:5] + shift                   # planted signal: five features +1.5 SD in the positive class
  list(X = X, y = y, units = rownames(X))
}

blocked_design <- function(seed, n_subj = 10L, reps = 2L, p = 200L) {
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  subj <- sprintf("s%02d", seq_len(2L * n_subj)); units <- rep(subj, each = reps)
  S <- matrix(stats::rnorm(length(subj) * p), length(subj), p, dimnames = list(subj, NULL))
  X <- S[units, , drop = FALSE] + matrix(stats::rnorm(length(units) * p), length(units), p)
  dimnames(X) <- list(sprintf("%s_%d", units, rep(seq_len(reps), times = length(subj))), sprintf("F%03d", seq_len(p)))
  list(X = X, y = rep(rep(c(1L, 0L), each = n_subj), each = reps), units = units)
}

perm_p <- function(d, s, obs, B, seed) {
  null <- vapply(seq_len(B), function(b) {
    f(".bm_use_stream")(f(".bm_stream")(seed, b)); yb <- f("pd_bm_permute_labels")(d$y, d$units)
    f("pd_bm_nested")(d$X, yb, d$units, s, .Random.seed, record = FALSE)$pooled }, numeric(1))
  (sum(null >= obs - 1e-12) + 1) / (B + 1)
}

t0 <- Sys.time()
results <- list()
for (fam in families) {
  k <- fam$index; cl <- classifier_of(fam); off <- 1000L * k
  message(sprintf("family %s", cl))
  # 1. nested AUC under the null
  auc <- unlist(lapply_cores(seq_len(rc$nested_auc_null), function(i) { d <- noise(sd_$nested_auc_null_base + off + i)
    f("pd_bm_nested")(d$X, d$y, d$units, settings(d$units, cl), f(".bm_stream")(sd_$nested_auc_null_base + off + i, 0L), record = FALSE)$pooled }))
  results$nested_auc_null <- rbind(results$nested_auc_null, data.frame(family = cl, dataset = seq_along(auc), seed = sd_$nested_auc_null_base + off + seq_along(auc), pooled_oof_auc = auc))
  # 2. permutation P under the null (whole-procedure permutation, B from the pre-registration)
  B <- rc$permutation_B
  perm <- unlist(lapply_cores(seq_len(rc$permutation_p_null_datasets), function(i) { seed <- sd_$permutation_null_base + off + i
    d <- noise(seed); s <- settings(d$units, cl)
    obs <- f("pd_bm_nested")(d$X, d$y, d$units, s, f(".bm_stream")(seed, 0L), record = FALSE)$pooled
    perm_p(d, s, obs, B, seed) }))
  results$permutation_p_null <- rbind(results$permutation_p_null, data.frame(family = cl, dataset = seq_along(perm), seed = sd_$permutation_null_base + off + seq_along(perm), permutation_p = perm, B = B))
  # 3. leaky reference: the panel is chosen on all data before cross-validation (the leak the gate must detect)
  leaky <- unlist(lapply_cores(seq_len(rc$leaky_reference), function(i) { seed <- sd_$leaky_base + off + i; d <- noise(seed)
    a <- f("pd_auc_columns")(scale(d$X), d$y); panel <- names(sort(-abs(a - 0.5), method = "radix"))[1:5]
    f("pd_bm_nested")(d$X, d$y, d$units, settings(d$units, cl, fixed = panel), f(".bm_stream")(seed, 0L), record = FALSE)$pooled }))
  results$leaky_reference <- rbind(results$leaky_reference, data.frame(family = cl, dataset = seq_along(leaky), seed = sd_$leaky_base + off + seq_along(leaky), leaky_auc = leaky))
  # 4. planted signal: nested AUC and permutation P
  planted <- lapply_cores(seq_len(rc$planted_auc), function(i) { seed <- sd_$planted_auc_base + off + i; d <- noise(seed, shift = 1.5); s <- settings(d$units, cl)
    obs <- f("pd_bm_nested")(d$X, d$y, d$units, s, f(".bm_stream")(seed, 0L), record = FALSE)$pooled; c(obs, NA_real_) })
  results$planted_auc <- rbind(results$planted_auc, data.frame(family = cl, dataset = seq_along(planted), seed = sd_$planted_auc_base + off + seq_along(planted), pooled_oof_auc = vapply(planted, `[`, 0, 1)))
  pperm <- unlist(lapply_cores(seq_len(rc$planted_permutation_datasets), function(i) { seed <- sd_$planted_permutation_base + off + i
    d <- noise(seed, shift = 1.5); s <- settings(d$units, cl)
    obs <- f("pd_bm_nested")(d$X, d$y, d$units, s, f(".bm_stream")(seed, 0L), record = FALSE)$pooled
    perm_p(d, s, obs, B, seed) }))
  results$planted_permutation_p <- rbind(results$planted_permutation_p, data.frame(family = cl, dataset = seq_along(pperm), seed = sd_$planted_permutation_base + off + seq_along(pperm), permutation_p = pperm, B = B))
}

# 5. LOOCV outer scheme (12 vs 12 units, inner 3-fold), ridge and lasso
for (cl in gates$families[c(1, 2)] |> lapply(classifier_of) |> unlist()) {
  off <- 1000L * match(cl, sapply(families, classifier_of)) - 1000L
  null <- unlist(lapply_cores(seq_len(rc$loocv_null), function(i) { d <- noise(sd_$loocv_null_base + off + i, n = 12L)
    f("pd_bm_nested")(d$X, d$y, d$units, settings(d$units, cl, scheme = "loocv", outer_k = 5L, repeats = 1L), f(".bm_stream")(sd_$loocv_null_base + off + i, 0L), record = FALSE)$pooled }))
  results$loocv_null <- rbind(results$loocv_null, data.frame(family = cl, dataset = seq_along(null), seed = sd_$loocv_null_base + off + seq_along(null), pooled_oof_auc = null))
  planted <- unlist(lapply_cores(seq_len(rc$loocv_planted), function(i) { d <- noise(sd_$loocv_planted_base + off + i, n = 12L, shift = 1.5)
    f("pd_bm_nested")(d$X, d$y, d$units, settings(d$units, cl, scheme = "loocv", outer_k = 5L, repeats = 1L), f(".bm_stream")(sd_$loocv_planted_base + off + i, 0L), record = FALSE)$pooled }))
  results$loocv_planted <- rbind(results$loocv_planted, data.frame(family = cl, dataset = seq_along(planted), seed = sd_$loocv_planted_base + off + seq_along(planted), pooled_oof_auc = planted))
}

# 6. subject-blocked units (10 vs 10 subjects x 2 observations): grouped and ungrouped AUC, whole-subject permutation P
for (fam in families[c(1, 2, 3)]) {
  k <- fam$index; cl <- classifier_of(fam); off <- 1000L * k
  res <- lapply_cores(seq_len(rc$blocked_auc_null), function(i) { seed <- sd_$blocked_base + off + i; d <- blocked_design(seed); s <- settings(rownames(d$X), cl, repeats = 2L)
    grouped <- f("pd_bm_nested")(d$X, d$y, d$units, s, f(".bm_stream")(seed, 0L), record = FALSE)$pooled
    ungrouped <- f("pd_bm_nested")(d$X, d$y, rownames(d$X), s, f(".bm_stream")(seed, 0L), record = FALSE)$pooled
    c(grouped, ungrouped) })
  m <- do.call(rbind, res)
  results$blocked_nested_auc_null <- rbind(results$blocked_nested_auc_null, data.frame(family = cl, dataset = seq_len(nrow(m)), seed = sd_$blocked_base + off + seq_len(nrow(m)),
                                                                                       pooled_oof_auc = m[, 1], ungrouped_auc = m[, 2]))
  bperm <- unlist(lapply_cores(seq_len(rc$blocked_permutation_datasets), function(i) { seed <- sd_$blocked_base + off + i; d <- blocked_design(seed)
    s <- settings(rownames(d$X), cl, repeats = 2L); obs <- f("pd_bm_nested")(d$X, d$y, d$units, s, f(".bm_stream")(seed, 0L), record = FALSE)$pooled
    perm_p(d, s, obs, B, seed) }))
  results$blocked_permutation_p_null <- rbind(results$blocked_permutation_p_null, data.frame(family = cl, dataset = seq_along(bperm), seed = sd_$blocked_base + off + seq_along(bperm), permutation_p = bperm, B = B))
}

dir.create(out, recursive = TRUE, showWarnings = FALSE)
names_out <- c(nested_auc_null = "nested_auc_null.tsv", permutation_p_null = "permutation_p_null.tsv", leaky_reference = "leaky_reference.tsv",
               planted_auc = "planted_auc.tsv", planted_permutation_p = "planted_permutation_p.tsv", loocv_null = "loocv_null.tsv",
               loocv_planted = "loocv_planted.tsv", blocked_nested_auc_null = "blocked_nested_auc_null.tsv",
               blocked_permutation_p_null = "blocked_permutation_p_null.tsv")
for (nm in names(names_out)) if (!is.null(results[[nm]])) write_tsv(results[[nm]], names_out[[nm]])
jsonlite::write_json(list(gates = gates_path, seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")), cores = cores,
                          r_version = as.character(getRversion()), glmnet = as.character(utils::packageVersion("glmnet")),
                          e1071 = as.character(utils::packageVersion("e1071")), randomForest = as.character(utils::packageVersion("randomForest")),
                          proteomicsCore = as.character(utils::packageVersion("proteomicsCore"))),
                     file.path(out, "calibration_raw.json"), auto_unbox = TRUE, pretty = TRUE)
