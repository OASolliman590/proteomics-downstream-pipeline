# Design stage handler (packet R04: FR-031 to FR-040).
#
# Receives already-compiled numeric design matrices (from the safe
# declarative grammar; nothing is evaluated), writes them as plan
# artifacts and computes rank/alias, contrast estimability and direction,
# featurewise n/df/estimability, the frozen duplicateCorrelation consensus
# for repeated designs and the unscaled contrast covariance.

.pc_design_matrix <- function(spec) {
  X <- matrix(as.numeric(unlist(spec$matrix)), nrow = length(spec$observation_ids), byrow = TRUE,
              dimnames = list(unlist(spec$observation_ids), unlist(spec$coefficients)))
  X
}

design_stage <- function(request) .pc_run_stage(request, function(out) {
  p <- request$parameters
  outputs <- list()
  emit <- function(relative, id, type) outputs[[length(outputs) + 1L]] <<- .pc_output(out, relative, id, type)
  write_tsv <- function(df, relative, id, type) { .pc_write_tsv(df, file.path(out, relative)); emit(relative, id, type) }
  write_json <- function(value, relative, id, type) { .pc_write_json(value, file.path(out, relative)); emit(relative, id, type) }
  scope <- if (is.null(p$scope)) "analysis" else p$scope
  values <- .pc_matrix_from_tsv(.pc_find_input(request, "primary_matrix"), "numeric")
  observations <- .pc_read_tsv(.pc_find_input(request, "primary_observations"))
  observed_path <- .pc_find_input(request, "primary_observed_mask", required = FALSE)
  observed <- if (is.null(observed_path)) NULL else .pc_matrix_from_tsv(observed_path, "logical")
  coverage_path <- .pc_find_input(request, "coverage", required = FALSE)
  coverage <- if (is.null(coverage_path)) data.frame(model_id = character(), contrast_id = character(), feature_id = character(), eligibility = character(), reason = character()) else .pc_read_tsv(coverage_path)
  groups_all <- stats::setNames(observations[[p$group_column]], observations$observation_id)

  design_rows <- list(); blocking_rows <- list(); designs <- list()
  for (spec in p$designs) {
    X <- .pc_design_matrix(spec)
    if (!identical(rownames(X), colnames(values))) stop(sprintf("E_ID_ALIGNMENT: design %s rows do not align with the primary matrix observations", spec$design_id), call. = FALSE)
    designs[[spec$design_id]] <- list(X = X, spec = spec)
    base <- file.path("designs", spec$design_id)
    write_tsv(data.frame(observation_id = rownames(X), X, check.names = FALSE, stringsAsFactors = FALSE), file.path(base, "design_matrix.tsv"), paste0("design_", spec$design_id), "Design")
    write_json(list(design_id = spec$design_id, coefficients = I(unlist(spec$coefficients)), term_map = spec$term_map, centers = spec$centers,
                    blocking = spec$blocking[c("mode", "subject_column")], intercept = spec$intercept,
                    encoding = "group.<level>, continuous.<column>, factor.<column>.<level>; percent-encoded UTF-8; interactions use reference-omitted factor columns joined by ':'"),
               file.path(base, "term_map.json"), paste0("term_map_", spec$design_id), "Design")
    rank <- design_rank(X, if (is.null(p$rank_tolerance)) 1e-7 else p$rank_tolerance)
    design_rows[[length(design_rows) + 1L]] <- list(design_id = spec$design_id, n_observations = nrow(X), n_coefficients = ncol(X), coefficients = I(colnames(X)),
                                                    rank = rank$rank, full_rank = rank$full_rank, aliased = I(rank$aliased), alias_relations = rank$alias_relations,
                                                    condition_number = rank$condition_number, blocking_mode = spec$blocking$mode)
    block <- NULL; correlation <- NULL
    if (identical(spec$blocking$mode, "duplicate_correlation") && rank$full_rank) {
      if (!requireNamespace("limma", quietly = TRUE)) stop("E_DEPENDENCY_UNAVAILABLE: limma is required for duplicateCorrelation", call. = FALSE)
      block <- unlist(spec$blocking$subjects)
      usable <- if (is.null(observed)) !is.na(values) else (!is.na(values) & observed)
      Y <- values; Y[!usable] <- NA
      dc <- suppressWarnings(limma::duplicateCorrelation(Y, X, block = block))
      correlation <- dc$consensus.correlation
      blocking_rows[[length(blocking_rows) + 1L]] <- list(design_id = spec$design_id, mode = "duplicate_correlation", subject_column = spec$blocking$subject_column,
        consensus_correlation = correlation, n_blocks = length(unique(block)), estimation_universe = "primary matrix, genuinely observed cells", function_called = "limma::duplicateCorrelation",
        limma_version = as.character(utils::packageVersion("limma")))
    } else if (!identical(spec$blocking$mode, "none")) {
      blocking_rows[[length(blocking_rows) + 1L]] <- list(design_id = spec$design_id, mode = spec$blocking$mode, subject_column = spec$blocking$subject_column,
        consensus_correlation = NULL, n_blocks = length(unique(unlist(spec$blocking$subjects))), estimation_universe = NULL, function_called = NULL)
    }
    designs[[spec$design_id]]$block <- block; designs[[spec$design_id]]$correlation <- correlation; designs[[spec$design_id]]$rank <- rank
  }

  contrast_rows <- list(); contrast_long <- list()
  for (contrast in p$contrasts) {
    d <- designs[[contrast$design_id]]; X <- d$X; w <- as.numeric(unlist(contrast$weights)); required <- unlist(contrast$required_groups)
    estimable <- contrast_estimable(X, w, if (is.null(p$estimability_tolerance)) 1e-8 else p$estimability_tolerance)
    implied <- implied_group_weights(X, w, d$spec$term_map, unlist(d$spec$group_levels))
    direction <- "not_applicable"
    if (implied$representable && length(required) == 2L) {
      wa <- implied$weights[[required[1]]]; wb <- implied$weights[[required[2]]]
      direction <- if (wa > 0 && wb < 0) "consistent" else if (wa < 0 && wb > 0) "mismatch" else "not_applicable"
    }
    units <- vapply(required, function(g) sum(groups_all == g), integer(1))
    subjects <- if (!is.null(d$spec$blocking$subject_column)) vapply(required, function(g) length(unique(unlist(d$spec$blocking$subjects)[groups_all == g])), integer(1)) else units
    contrast_rows[[length(contrast_rows) + 1L]] <- list(contrast_id = contrast$contrast_id, design_id = contrast$design_id, role = contrast$role, required_groups = I(required),
      estimable = estimable, implied_group_weights = if (implied$representable) as.list(implied$weights) else NULL, direction_check = direction,
      units_by_required_group = as.list(stats::setNames(as.integer(units), required)), subjects_by_required_group = as.list(stats::setNames(as.integer(subjects), required)))
    contrast_long[[length(contrast_long) + 1L]] <- data.frame(contrast_id = contrast$contrast_id, design_id = contrast$design_id, coefficient = colnames(X), weight = w, stringsAsFactors = FALSE)
  }
  if (length(contrast_long)) write_tsv(do.call(rbind, contrast_long), "contrasts.tsv", "contrast_weights", "Design")

  # unscaled contrast covariance for each full-rank design (complete design, observation weights if declared)
  for (did in names(designs)) {
    d <- designs[[did]]
    if (!d$rank$full_rank) next
    cs <- Filter(function(c) identical(c$design_id, did), p$contrasts)
    if (!length(cs)) next
    C <- sapply(cs, function(c) as.numeric(unlist(c$weights))); C <- matrix(C, ncol = length(cs)); colnames(C) <- vapply(cs, function(c) c$contrast_id, "")
    model_weights <- NULL
    for (m in p$models) if (identical(m$design_id, did) && !is.null(m$weights) && identical(m$weights$kind, "observation")) model_weights <- as.numeric(unlist(m$weights$values))
    cov <- contrast_unscaled_covariance(d$X, C, model_weights, d$block, d$correlation)
    long <- expand.grid(contrast_a = colnames(C), contrast_b = colnames(C), stringsAsFactors = FALSE)
    long$unscaled_covariance <- as.vector(cov)
    write_tsv(long, file.path("designs", did, "contrast_covariance.tsv"), paste0("contrast_covariance_", did), "DesignDiagnostics")
  }

  # featurewise estimability (analysis scope)
  if (identical(scope, "analysis")) {
    est <- list()
    for (m in p$models) {
      d <- designs[[m$design_id]]
      if (!d$rank$full_rank) next
      cs <- Filter(function(c) identical(c$design_id, m$design_id), p$contrasts)
      est[[length(est) + 1L]] <- featurewise_estimability(values, observed, d$X, cs, coverage, m$model_id, groups_all)
      if (!is.null(m$weights)) {
        wv <- m$weights$values
        wdf <- if (identical(m$weights$kind, "observation")) data.frame(observation_id = colnames(values), weight = as.numeric(unlist(wv)), stringsAsFactors = FALSE)
               else data.frame(feature_id = rownames(values), matrix(as.numeric(unlist(wv)), nrow = nrow(values), byrow = TRUE, dimnames = list(NULL, colnames(values))), check.names = FALSE)
        write_tsv(wdf, file.path("weights", paste0(m$model_id, ".tsv")), paste0("weights_", m$model_id), "WeightValidation")
      }
    }
    if (length(est)) write_tsv(do.call(rbind, est), "estimability.tsv", "estimability", "DesignDiagnostics")
  }
  availability <- lapply(unique(vapply(p$models, function(m) m$engine, "")), function(e) operational_availability(e))
  write_json(list(designs = design_rows, contrasts = contrast_rows, blocking = blocking_rows, operational_availability = availability,
                  note = "Rank/alias/estimability computed by QR/SVD on explicit numeric matrices; no formula was evaluated."),
             "design_diagnostics.json", "design_diagnostics", "DesignDiagnostics")
  list(outputs = outputs, warnings = list(), message = sprintf("design diagnostics for %d design(s) and %d contrast(s)", length(designs), length(p$contrasts)))
})
