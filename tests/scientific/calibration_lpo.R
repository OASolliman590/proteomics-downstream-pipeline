# R15c follow-up (D-77): calibration of the leave-pair-out AUC that the loocv scheme now reports (V181).  Runs the
# production nested procedure (pd_bm_nested with scheme loocv) on the null and planted-signal generators of the R15c
# driver at 12 vs 12 units.  Every count and seed comes from
# docs/validation/016-methods-extension-validation/R15c/followup/gates-preregistration.json, committed before this run.
#
# Usage: Rscript --vanilla tests/scientific/calibration_lpo.R <gates.json> <out_dir> <cores>
# Writes lpo_null.tsv, lpo_planted.tsv and lpo_raw.json (seeds, counts, R and package versions, timing).

args <- commandArgs(TRUE)
gates_path <- args[1]; out <- args[2]; cores <- as.integer(args[3])
gates <- jsonlite::fromJSON(gates_path, simplifyVector = FALSE)
ns <- asNamespace("proteomicsCore"); f <- function(name) get(name, envir = ns)
rc <- gates$replicate_counts; sd_ <- gates$seeds
dir.create(out, recursive = TRUE, showWarnings = FALSE)
started <- Sys.time()
lapply_cores <- function(x, fun) {
  res <- if (cores > 1L && .Platform$OS.type == "unix") parallel::mclapply(x, fun, mc.cores = cores) else lapply(x, fun)
  if (any(vapply(res, inherits, logical(1), "try-error")) || any(vapply(res, is.null, logical(1))))
    stop("E_CALIBRATION_REPLICATE: a calibration replicate failed; the run is not accepted", call. = FALSE)
  unlist(res)
}

settings <- function(ids, classifier) list(
  positive = "P", negative = "N", levels = c("P", "N"), imputation = "none", selection = "top_k_auc", k = 5L, alpha = 0.5,
  classifier = classifier, lambda = 0.1, cost = NULL, threshold_rule = "youden_train", scheme = "loocv", outer_k = 5L,
  repeats = 1L, inner_k = 3L, observation_ids = ids, fixed_features = NULL)

noise <- function(seed, n = 12L, p = 200L, shift = 0) {   # identical to the R15c driver's generators
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  X <- matrix(stats::rnorm(2L * n * p), 2L * n, p, dimnames = list(sprintf("u%02d", seq_len(2L * n)), sprintf("F%03d", seq_len(p))))
  y <- rep(c(1L, 0L), each = n)
  if (shift > 0) X[y == 1L, 1:5] <- X[y == 1L, 1:5] + shift
  list(X = X, y = y, units = rownames(X))
}

index <- c(penalized_logistic = 0L, lasso = 1L)
results <- list(lpo_null = NULL, lpo_planted = NULL)
for (cl in unlist(rc$families_lpo)) {
  off <- 1000L * index[[cl]]
  run <- function(base, shift) function(i) { seed <- base + off + i; d <- noise(seed, shift = shift)
    r <- f("pd_bm_nested")(d$X, d$y, d$units, settings(d$units, cl), f(".bm_stream")(seed, 0L), record = FALSE)
    stopifnot(identical(r$auc_estimator, "leave_pair_out")); r$pooled }
  null <- lapply_cores(seq_len(rc$lpo_null), run(sd_$lpo_null_base, 0))
  results$lpo_null <- rbind(results$lpo_null, data.frame(family = cl, dataset = seq_along(null), seed = sd_$lpo_null_base + off + seq_along(null), lpo_auc = null))
  planted <- lapply_cores(seq_len(rc$lpo_planted), run(sd_$lpo_planted_base, 1.5))
  results$lpo_planted <- rbind(results$lpo_planted, data.frame(family = cl, dataset = seq_along(planted), seed = sd_$lpo_planted_base + off + seq_along(planted), lpo_auc = planted))
}
for (name in names(results))
  utils::write.table(results[[name]], file.path(out, paste0(name, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE, eol = "\n")
raw <- list(gates = "docs/validation/016-methods-extension-validation/R15c/followup/gates-preregistration.json",
            gates_sha256 = paste(unclass(openssl::sha256(file(gates_path))), collapse = ""), seconds = round(as.numeric(difftime(Sys.time(), started, units = "secs")), 3),
            cores = cores, r_version = paste(R.version$major, R.version$minor, sep = "."),
            glmnet = as.character(utils::packageVersion("glmnet")), proteomicsCore = as.character(utils::packageVersion("proteomicsCore")))
writeBin(charToRaw(paste0(as.character(jsonlite::toJSON(raw, auto_unbox = TRUE, pretty = TRUE)), "\n")), file.path(out, "lpo_raw.json"))
cat("lpo calibration:", nrow(results$lpo_null), "null and", nrow(results$lpo_planted), "planted replicates,", raw$seconds, "s\n")
