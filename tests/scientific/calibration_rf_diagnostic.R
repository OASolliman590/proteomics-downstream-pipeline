# R15c follow-up (D-78): random-forest permutation-P diagnostic of the V180 FAIL.  Runs the production nested procedure
# (pd_bm_nested, R14d) and the production permutation rule ((k + 1) / (B + 1), k = permutations with AUC >= observed)
# on null data, and records ties and the mid-P.  Every count and seed comes from
# docs/validation/016-methods-extension-validation/R15c/followup/gates-preregistration.json, committed before this run.
#
# Usage: Rscript --vanilla tests/scientific/calibration_rf_diagnostic.R <gates.json> <out_dir> <cores>
# Writes rf_permutation_diagnostic.tsv and rf_diagnostic_raw.json (seeds, counts, R and package versions, timing).

args <- commandArgs(TRUE)
gates_path <- args[1]; out <- args[2]; cores <- as.integer(args[3])
gates <- jsonlite::fromJSON(gates_path, simplifyVector = FALSE)
ns <- asNamespace("proteomicsCore"); f <- function(name) get(name, envir = ns)
rc <- gates$replicate_counts; sd_ <- gates$seeds
dir.create(out, recursive = TRUE, showWarnings = FALSE)
started <- Sys.time()

settings <- function(ids) list(
  positive = "P", negative = "N", levels = c("P", "N"), imputation = "none", selection = "top_k_auc", k = 5L, alpha = 0.5,
  classifier = "random_forest", lambda = 0.1, cost = NULL, threshold_rule = "youden_train", scheme = "repeated_stratified_kfold", outer_k = 5L,
  repeats = 2L, inner_k = 3L, observation_ids = ids, fixed_features = NULL)

noise <- function(seed, n = 15L, p = 200L) {   # identical to the R15c driver's null generator
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  X <- matrix(stats::rnorm(2L * n * p), 2L * n, p, dimnames = list(sprintf("u%02d", seq_len(2L * n)), sprintf("F%03d", seq_len(p))))
  list(X = X, y = rep(c(1L, 0L), each = n), units = rownames(X))
}

B <- as.integer(rc$rf_diagnostic_B)
base <- sd_$rf_diagnostic_base + 1000L * 5L   # random_forest is family index 5
one <- function(i) {
  seed <- base + i; d <- noise(seed); s <- settings(d$units)
  obs <- f("pd_bm_nested")(d$X, d$y, d$units, s, f(".bm_stream")(seed, 0L), record = FALSE)$pooled
  null <- vapply(seq_len(B), function(b) {   # the production permutation job (post_de_biomarker_stage, FR-156)
    f(".bm_use_stream")(f(".bm_stream")(seed, b)); yb <- f("pd_bm_permute_labels")(d$y, d$units)
    f("pd_bm_nested")(d$X, yb, d$units, s, .Random.seed, record = FALSE)$pooled }, numeric(1))
  ties <- sum(abs(null - obs) <= 1e-12); greater <- sum(null > obs + 1e-12)
  data.frame(family = "random_forest", dataset = i, seed = seed, observed_auc = obs, n_greater = greater, n_tied = ties,
             permutation_p = (greater + ties + 1) / (B + 1), mid_p = (greater + 0.5 * (ties + 1)) / (B + 1), B = B)
}
idx <- seq_len(rc$rf_diagnostic_datasets)
res <- if (cores > 1L && .Platform$OS.type == "unix") parallel::mclapply(idx, one, mc.cores = cores) else lapply(idx, one)
if (any(vapply(res, inherits, logical(1), "try-error")) || any(vapply(res, is.null, logical(1))))
  stop("E_CALIBRATION_REPLICATE: a calibration replicate failed; the run is not accepted", call. = FALSE)
utils::write.table(do.call(rbind, res), file.path(out, "rf_permutation_diagnostic.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, eol = "\n")
raw <- list(gates = "docs/validation/016-methods-extension-validation/R15c/followup/gates-preregistration.json",
            gates_sha256 = as.character(openssl::sha256(file(gates_path))), seconds = round(as.numeric(difftime(Sys.time(), started, units = "secs")), 3),
            cores = cores, r_version = paste(R.version$major, R.version$minor, sep = "."),
            randomForest = as.character(utils::packageVersion("randomForest")), proteomicsCore = as.character(utils::packageVersion("proteomicsCore")))
writeBin(charToRaw(paste0(as.character(jsonlite::toJSON(raw, auto_unbox = TRUE, pretty = TRUE)), "\n")), file.path(out, "rf_diagnostic_raw.json"))
cat("rf diagnostic:", length(idx), "datasets, B =", B, ",", raw$seconds, "s\n")
