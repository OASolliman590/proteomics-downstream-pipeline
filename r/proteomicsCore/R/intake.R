# Canonical bundle readers and typed stage-result helpers (packet R02).
#
# Later R stages read the R02 canonical protein bundle only through these
# functions, which verify the SHA-256 manifest and the matrix/mask/metadata
# key alignment before returning any value.  The helpers also build typed
# stage results so every handler reports COMPLETED/FAILED/NOT_RUN honestly.

.pc_read_tsv <- function(path) {
  utils::read.delim(path, header = TRUE, sep = "\t", quote = "\"", comment.char = "", check.names = FALSE,
                    na.strings = character(), stringsAsFactors = FALSE, colClasses = "character",
                    encoding = "UTF-8", fill = FALSE, blank.lines.skip = TRUE, strip.white = FALSE)
}

.pc_fmt <- function(x) {
  out <- rep("NA", length(x))
  # Infinite values (for example an infinite empirical-Bayes prior df) are written literally; NaN is never written.
  out[!is.na(x) & is.infinite(x)] <- ifelse(x[!is.na(x) & is.infinite(x)] > 0, "Inf", "-Inf")
  ok <- !is.na(x) & is.finite(x)
  short <- sprintf("%.15g", x[ok])
  exact <- suppressWarnings(as.numeric(short)) == x[ok]
  short[!exact] <- sprintf("%.17g", x[ok][!exact])
  out[ok] <- short
  out
}

.pc_fmt_bool <- function(x) ifelse(is.na(x), "NA", ifelse(x, "true", "false"))

.pc_write_tsv <- function(df, path) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  for (name in names(df)) {
    column <- df[[name]]
    if (is.logical(column)) df[[name]] <- .pc_fmt_bool(column)
    else if (is.numeric(column)) df[[name]] <- .pc_fmt(as.numeric(column))
    else { column <- as.character(column); column[is.na(column)] <- "NA"; df[[name]] <- column }
  }
  special <- function(values) any(grepl("[\t\n\r\"]", values))
  quote_cols <- which(vapply(df, special, logical(1)) | special(names(df)))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  connection <- file(path, open = "wb", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)
  utils::write.table(df, connection, sep = "\t", quote = if (length(quote_cols)) quote_cols else FALSE, qmethod = "double",
                     row.names = FALSE, col.names = TRUE, na = "NA", eol = "\n", fileEncoding = "UTF-8")
  invisible(path)
}

.pc_matrix_from_tsv <- function(path, type = c("numeric", "logical")) {
  type <- match.arg(type)
  table <- .pc_read_tsv(path)
  ids <- table[[1L]]
  body <- table[-1L]
  if (type == "numeric") {
    values <- vapply(body, function(column) { missing <- column == "NA"; out <- suppressWarnings(as.numeric(column)); if (any(is.na(out) & !missing)) stop(sprintf("E_INTEGRITY: non-numeric value in %s", basename(path)), call. = FALSE); out[missing] <- NA_real_; out }, numeric(nrow(table)))
  } else {
    values <- vapply(body, function(column) { if (!all(column %in% c("true", "false"))) stop(sprintf("E_INTEGRITY: invalid mask value in %s", basename(path)), call. = FALSE); column == "true" }, logical(nrow(table)))
  }
  values <- matrix(values, nrow = nrow(table), dimnames = list(ids, names(body)))
  values
}

.pc_verify_hash <- function(path, expected, label = basename(path)) {
  if (!file.exists(path)) stop(sprintf("E_INTEGRITY: missing artifact %s", label), call. = FALSE)
  if (!identical(sha256_file(path), expected)) stop(sprintf("E_INTEGRITY: SHA-256 mismatch for %s", label), call. = FALSE)
  invisible(TRUE)
}

#' Read and verify an R02 canonical bundle directory.
read_canonical_bundle <- function(bundle_dir) {
  manifest_path <- file.path(bundle_dir, "manifest.json")
  if (!file.exists(manifest_path)) stop("E_INTEGRITY: canonical bundle manifest is missing", call. = FALSE)
  manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
  artifacts <- list()
  for (artifact in manifest$artifacts) {
    path <- file.path(bundle_dir, artifact$relative_path)
    .pc_verify_hash(path, artifact$sha256, artifact$relative_path)
    artifacts[[artifact$artifact_id]] <- path
  }
  features <- unlist(manifest$feature_ids); observations <- unlist(manifest$observation_ids)
  load_matrix <- function(id, type) {
    if (is.null(artifacts[[id]])) return(NULL)
    value <- .pc_matrix_from_tsv(artifacts[[id]], type)
    if (!identical(rownames(value), features) || !identical(colnames(value), observations)) stop(sprintf("E_INTEGRITY: %s keys do not match the manifest", id), call. = FALSE)
    value
  }
  values <- load_matrix("matrix", "numeric")
  numeric_mask <- load_matrix("numeric_mask", "logical")
  observed <- load_matrix("observed_mask", "logical")
  prior <- load_matrix("prior_imputed_mask", "logical")
  if (!identical(!is.na(values), numeric_mask)) stop("E_INTEGRITY: numeric mask disagrees with matrix availability", call. = FALSE)
  if (identical(manifest$original_observed_mask$state, "known") && (is.null(observed) || is.null(prior))) stop("E_INTEGRITY: known original-observed mask is missing", call. = FALSE)
  obs_table <- .pc_read_tsv(artifacts$observations)
  feature_table <- .pc_read_tsv(artifacts$features)
  if (!identical(obs_table$observation_id, observations) || !identical(feature_table$feature_id, features)) stop("E_INTEGRITY: metadata keys do not match the manifest", call. = FALSE)
  obs_table[obs_table == "NA"] <- NA_character_
  list(values = values, numeric_mask = numeric_mask, observed = observed, prior_imputed = prior,
       observations = obs_table, features = feature_table, manifest = manifest,
       scale = manifest$scale, mask_state = manifest$original_observed_mask$state)
}

.pc_find_input <- function(request, artifact_id, required = TRUE) {
  matches <- Filter(function(item) identical(item$artifact_id, artifact_id), request$inputs)
  if (length(matches) == 0L) { if (required) stop(sprintf("E_REFERENCE_UNKNOWN: input artifact %s", artifact_id), call. = FALSE); return(NULL) }
  if (length(matches) > 1L) stop(sprintf("E_ID_DUPLICATE: input artifact %s", artifact_id), call. = FALSE)
  item <- matches[[1L]]
  if (!is.character(item$path) || length(item$path) != 1L || !nzchar(item$path)) stop(sprintf("E_REFERENCE_UNKNOWN: input artifact %s has no usable path", artifact_id), call. = FALSE)
  path <- normalizePath(item$path, winslash = "/", mustWork = FALSE)   # forward slashes on every platform (CI run 36982402784)
  .pc_verify_hash(path, item$sha256, artifact_id)
  path
}

.pc_now <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

.pc_exit_for <- function(code) {
  if (code %in% c("E_INTEGRITY", "E_PLAN_CHANGED", "E_SOURCE_CHANGED", "E_PLAN_REQUIRED")) return(5L)
  if (code %in% c("E_CAPABILITY_NOT_AVAILABLE", "E_ENGINE_NOT_AVAILABLE", "E_DEPENDENCY_UNAVAILABLE")) return(3L)
  if (code %in% c("E_NUMERICAL", "E_ENGINE_FAILED", "E_STAGE_EXCEPTION", "E_PERMANOVA_IDENTITY")) return(4L)
  2L
}

.pc_output <- function(output_dir, relative_path, artifact_id, result_type) {
  list(artifact_id = artifact_id, relative_path = relative_path, sha256 = sha256_file(file.path(output_dir, relative_path)), result_type = result_type)
}

# Run a handler body, capturing warnings and converting typed "E_CODE: message"
# errors into a FAILED or NOT_RUN result.  sessionInfo is written in-process.
.pc_run_stage <- function(request, body) {
  started <- .pc_now()
  output_dir <- normalizePath(request$output_temp_dir, mustWork = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  warning_messages <- character()
  set.seed(as.integer(request$rng$seed), kind = "L'Ecuyer-CMRG")
  outcome <- tryCatch(withCallingHandlers(body(output_dir), warning = function(w) { warning_messages <<- c(warning_messages, conditionMessage(w)); invokeRestart("muffleWarning") }),
                      error = function(e) structure(list(message = conditionMessage(e)), class = "pc_failure"))
  session_path <- file.path(output_dir, "sessionInfo.txt")
  utils::capture.output(utils::sessionInfo(), file = session_path)
  warnings <- lapply(unique(warning_messages), function(m) list(code = "R_WARNING", stage_id = request$stage_id, severity = "warning", message = m, artifact_ref = NULL))
  base <- list(schema_version = "1.2.0", run_id = request$run_id, stage_id = request$stage_id, capability = request$capability,
               plan_hash = request$plan_hash, started_at = started, finished_at = .pc_now(), warnings = warnings, session_info_path = "sessionInfo.txt")
  session_output <- .pc_output(output_dir, "sessionInfo.txt", "session_info", "session_info")
  if (inherits(outcome, "pc_failure")) {
    code <- regmatches(outcome$message, regexpr("^E_[A-Z0-9_]+", outcome$message))
    if (length(code) == 0L) code <- "E_STAGE_EXCEPTION"
    exit <- .pc_exit_for(code)
    diagnostics <- if (!is.null(attr(outcome, "outputs"))) attr(outcome, "outputs") else list()
    return(c(base, list(state = if (exit == 3L) "NOT_RUN" else "FAILED", reason_code = code, message = outcome$message, exit_code = exit,
                        outputs = c(diagnostics, list(session_output)))))
  }
  outputs <- outcome$outputs
  extra_warnings <- outcome$warnings
  if (length(extra_warnings)) base$warnings <- c(base$warnings, extra_warnings)
  c(base, list(state = "COMPLETED", reason_code = NULL, message = outcome$message, exit_code = 0L, outputs = c(outputs, list(session_output))))
}

.pc_warning <- function(request, code, message, artifact_ref = NULL) list(code = code, stage_id = request$stage_id, severity = "warning", message = message, artifact_ref = artifact_ref)

.pc_param <- function(params, name, default = NULL) { value <- params[[name]]; if (is.null(value)) default else value }
