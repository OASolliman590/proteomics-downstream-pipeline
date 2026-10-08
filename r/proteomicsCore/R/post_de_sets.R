# Post-differential analysis: shared helpers and result structure (packet R14a, SM31, SM41, FR-131 to FR-136).
#
# Every post-DE stage reads the frozen plan and completed differential tables read-only and verifies their hashes
# (.pd_verify_inputs).  Sets are built only from completed families with the finite grammar validated at plan time
# (union / intersect / difference / complement_within_tested); regions are exact, overlap P values exist only for
# contrasts that share no biological units, and concordance is descriptive (no P column).

PD_CLAIM_LABELS <- c("descriptive", "exploratory_raw_p", "in_sample", "cross_validated_nested", "fixed_panel_cv", "independently_validated", "module_level")

# ----------------------------------------------------------------------------- shared helpers (used by R14a-R14e)
.pd_verify_inputs <- function(request) {
  plan_path <- .pc_find_input(request, "plan")
  plan <- jsonlite::fromJSON(plan_path, simplifyVector = FALSE)
  if (is.null(request$plan_hash) || !identical(plan$plan_hash, request$plan_hash)) stop("E_PLAN_CHANGED: request plan_hash does not match the frozen plan", call. = FALSE)
  planned <- list(); for (a in plan$artifacts) planned[[a$artifact_id]] <- a
  for (item in request$inputs) {
    if (identical(item$artifact_id, "plan")) next
    if (startsWith(item$artifact_id, "stage__")) {   # completed-stage output: its hash was recorded by that stage (verified on read)
      if (!identical(sha256_file(item$path), item$sha256)) stop(sprintf("E_INTEGRITY: %s changed after its stage completed", item$artifact_id), call. = FALSE)
      next
    }
    if (startsWith(item$artifact_id, "input__")) next   # declared external input (hash recorded at plan time, verified on read)
    if (startsWith(item$artifact_id, "resource_") || startsWith(item$artifact_id, "resfile_")) next   # SM14 snapshot verified by verify_resources and on read (A-2026-10-01-18)
    ref <- planned[[item$artifact_id]]
    if (is.null(ref) || !identical(ref$sha256, item$sha256)) stop(sprintf("E_PLAN_CHANGED: input %s differs from the frozen plan", item$artifact_id), call. = FALSE)
    if (identical(ref$result_type, "DisplayOnlyMatrix")) stop("E_DISPLAY_MATRIX_REJECTED: display-only matrices cannot feed post-DE analysis", call. = FALSE)
  }
  plan
}

# Primary matrix restricted to genuinely observed cells (never imputed), with aligned observations and features.
.pd_primary <- function(request) {
  Y <- .pc_matrix_from_tsv(.pc_find_input(request, "primary_matrix"), "numeric")
  observed <- .pc_matrix_from_tsv(.pc_find_input(request, "primary_observed_mask"), "logical")
  Y[!observed[rownames(Y), colnames(Y)]] <- NA
  obs <- .pc_read_tsv(.pc_find_input(request, "primary_observations"))
  obs <- obs[match(colnames(Y), obs$observation_id), , drop = FALSE]
  obs[obs == "NA"] <- NA_character_
  features_path <- .pc_find_input(request, "primary_features", required = FALSE)
  features <- if (is.null(features_path)) NULL else .pc_read_tsv(features_path)
  list(Y = Y, obs = obs, features = features)
}

.pd_read_stage_table <- function(request, artifact_id) {
  path <- .pc_find_input(request, artifact_id, required = FALSE)
  if (is.null(path)) return(NULL)
  table <- .pc_read_tsv(path)
  for (column in intersect(c("effect", "p_value", "q_value", "statistic", "effect_se", "df_inference", "ci_lower", "ci_upper"), names(table)))
    table[[column]] <- suppressWarnings(as.numeric(ifelse(table[[column]] == "NA", NA, table[[column]])))
  table
}

# Biological units of each observation: the subject in subject-blocked designs, otherwise the biological unit (SM03).
.pd_units <- function(obs, subject_column = NULL) {
  if (!is.null(subject_column) && subject_column %in% names(obs)) {
    s <- obs[[subject_column]]
    if (any(is.na(s) | s == "")) stop("E_SUBJECT_MISSING: subject identifiers are missing for a subject-blocked design", call. = FALSE)
    return(s)
  }
  if ("biological_unit_id" %in% names(obs)) ifelse(is.na(obs$biological_unit_id), obs$observation_id, obs$biological_unit_id) else obs$observation_id
}

.pd_record <- function(module, state, reason_code = NULL, reason = NULL, claim_label = NULL, details = NULL)
  list(module = module, state = state, reason_code = reason_code, reason = reason, claim_label = claim_label, details = details)

.pd_input_hashes <- function(request) lapply(request$inputs, function(i) list(artifact_id = i$artifact_id, sha256 = i$sha256))

# D-59: every adaptation the planner (or the stage) made is a warning with requested vs used and the reason; the full
# records are in eligibility.json (p$eligibility$adaptations) and in the R14f eligibility report.
.pd_adaptation_warnings <- function(request, adaptations = request$parameters$eligibility$adaptations) {
  lapply(adaptations, function(a) .pc_warning(request, "W_POST_DE_ADAPTED", sprintf("%s %s: requested %s, used %s (%s)", a$analysis, a$item, a$requested, a$used, a$reason)))
}

.pd_refusal_frame <- function(rows) if (length(rows)) do.call(rbind, rows) else
  data.frame(analysis = character(), item = character(), reason_code = character(), reason = character(), stringsAsFactors = FALSE)

.pd_check_claims <- function(table, column = "claim_label") {
  if (!nrow(table)) return(invisible(TRUE))
  if (!column %in% names(table)) stop("E_BIOMARKER_CLAIM: an output table has no claim label", call. = FALSE)
  bad <- setdiff(unique(table[[column]]), PD_CLAIM_LABELS)
  if (length(bad) || any(is.na(table[[column]]) | table[[column]] == "")) stop(sprintf("E_BIOMARKER_CLAIM: unknown or missing claim label(s): %s", paste(bad, collapse = ", ")), call. = FALSE)
  invisible(TRUE)
}

.pd_safe <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)

# ----------------------------------------------------------------------------- R14a leaves and rules
.pd_leaf_rows <- function(leaf, zero, treat) {
  table <- if (identical(leaf$hypothesis_type, "protein_treat")) treat else zero
  if (is.null(table)) stop(sprintf("E_SETRULE_UNKNOWN_CONTRAST: no %s table for family %s", leaf$hypothesis_type, leaf$family_id), call. = FALSE)
  rows <- table[table$model_id == leaf$model_id & table$contrast_id == leaf$contrast_id & table$hypothesis_type == leaf$hypothesis_type, , drop = FALSE]
  if (!nrow(rows)) stop(sprintf("E_SETRULE_UNKNOWN_CONTRAST: %s/%s has no completed rows", leaf$model_id, leaf$contrast_id), call. = FALSE)
  if (!identical(unique(rows$family_id), leaf$family_id)) stop(sprintf("E_SETRULE_UNKNOWN_CONTRAST: %s/%s rows do not belong to family %s", leaf$model_id, leaf$contrast_id, leaf$family_id), call. = FALSE)
  rows
}

# A leaf: tested features of one completed endpoint and the members meeting its declared criterion and direction.
pd_leaf <- function(leaf, zero, treat) {
  rows <- .pd_leaf_rows(leaf, zero, treat)
  tested <- rows$eligibility == "tested"
  value <- if (identical(leaf$criterion, "exploratory_raw_p")) rows$p_value else rows$q_value
  meets <- tested & !is.na(value) & value < leaf$threshold   # strict: "q < threshold" (V131)
  direction <- if (is.null(leaf$direction)) "any" else leaf$direction
  if (identical(direction, "up")) meets <- meets & rows$effect > 0
  if (identical(direction, "down")) meets <- meets & rows$effect < 0
  list(tested = rows$feature_id[tested], members = rows$feature_id[meets], effect = stats::setNames(rows$effect, rows$feature_id))
}

# Evaluate a validated rule tree; returns members and the tested universe (features tested in every referenced leaf).
pd_eval_rule <- function(rule, zero, treat) {
  if (!is.null(rule$family_id)) { l <- pd_leaf(rule, zero, treat); return(list(members = l$members, universe = l$tested)) }
  args <- lapply(rule$args, pd_eval_rule, zero = zero, treat = treat)
  universe <- Reduce(intersect, lapply(args, `[[`, "universe"))
  members <- switch(rule$op,
    union = Reduce(union, lapply(args, `[[`, "members")),
    intersect = Reduce(intersect, lapply(args, `[[`, "members")),
    difference = setdiff(args[[1]]$members, args[[2]]$members),
    complement_within_tested = setdiff(universe, args[[1]]$members),
    stop(sprintf("E_SETRULE_GRAMMAR: unknown operator %s", rule$op), call. = FALSE))
  list(members = members, universe = universe)
}

.pd_leaves <- function(rule) if (!is.null(rule$family_id)) list(rule) else do.call(c, lapply(rule$args, .pd_leaves))

# Exact exclusive regions of k sets: every non-empty in/out pattern (2^k - 1), counts summing to the union.
pd_regions <- function(sets, gene_map, all_regions = TRUE) {
  k <- length(sets); ids <- names(sets)
  universe <- sort(unique(unlist(sets)))
  member <- vapply(sets, function(s) universe %in% s, logical(length(universe)))
  member <- matrix(member, nrow = length(universe), dimnames = list(universe, ids))
  patterns <- if (all_regions) as.matrix(expand.grid(rep(list(c(TRUE, FALSE)), k)))[-(2^k), , drop = FALSE] else unique(member)
  rows <- list()
  for (r in seq_len(nrow(patterns))) {
    pat <- as.logical(patterns[r, ])
    hit <- universe[apply(member, 1L, function(m) identical(as.logical(m), pat))]
    genes <- sort(unique(unlist(gene_map[hit])))
    rows[[length(rows) + 1L]] <- data.frame(region_id = paste0("R", paste(as.integer(pat), collapse = "")), defining_sets = paste(ids[pat], collapse = ";"),
      excluded_sets = paste(ids[!pat], collapse = ";"), n_sets = sum(pat), n_protein_groups = length(hit), n_unique_genes = length(genes),
      n_without_gene = sum(vapply(gene_map[hit], length, integer(1)) == 0L), members = paste(hit, collapse = ";"), genes = paste(genes, collapse = ";"),
      claim_label = "descriptive", stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

# Direction-aware overlap of two single-endpoint sets (signs from their own effects).
pd_direction_overlap <- function(a, b, effect_a, effect_b) {
  both <- intersect(a, b)
  ea <- sign(effect_a[both]); eb <- sign(effect_b[both])
  list(concordant_up = both[ea > 0 & eb > 0], concordant_down = both[ea < 0 & eb < 0], discordant = both[ea != eb],
       only_a = setdiff(a, b), only_b = setdiff(b, a))
}

pd_concordance <- function(x, y, sig) {
  common <- intersect(names(x)[is.finite(x)], names(y)[is.finite(y)])
  x <- x[common]; y <- y[common]
  sig <- intersect(sig, common)
  list(n_common = length(common), pearson = if (length(common) > 2L) stats::cor(x, y) else NA_real_,
       spearman = if (length(common) > 2L) stats::cor(x, y, method = "spearman") else NA_real_,
       origin_slope = if (length(common)) sum(x * y) / sum(x * x) else NA_real_,
       n_significant_either = length(sig), n_sign_agree = sum(sign(x[sig]) == sign(y[sig])),
       sign_agreement = if (length(sig)) mean(sign(x[sig]) == sign(y[sig])) else NA_real_)
}

# Hypergeometric upper tail of the overlap within the universe tested in both (only for unit-disjoint contrasts).
pd_overlap_p <- function(a, b, universe) {
  a <- intersect(a, universe); b <- intersect(b, universe); k <- length(intersect(a, b))
  list(universe_n = length(universe), n_a = length(a), n_b = length(b), overlap = k, expected = length(a) * length(b) / max(1, length(universe)),
       p_value = stats::phyper(k - 1, length(a), length(universe) - length(a), length(b), lower.tail = FALSE))
}

# Observations and units used by a contrast: the observations of its required groups (SM03 units).
pd_contrast_units <- function(contrast, obs, group_column, subject_column) {
  keep <- obs[[group_column]] %in% unlist(contrast$required_groups)
  list(observations = obs$observation_id[keep], units = unique(.pd_units(obs[keep, , drop = FALSE], subject_column)),
       subjects = if (!is.null(subject_column) && subject_column %in% names(obs)) unique(stats::na.omit(obs[[subject_column]][keep])) else character(),
       biological_units = unique(stats::na.omit(obs$biological_unit_id[keep])))
}

.pd_gene_map <- function(features, ids) {
  if (is.null(features) || !"gene_ids" %in% names(features)) return(stats::setNames(rep(list(character()), length(ids)), ids))
  g <- lapply(features$gene_ids, function(x) if (is.na(x) || !nzchar(x)) character() else as.character(unlist(jsonlite::fromJSON(x))))
  out <- stats::setNames(g, features$feature_id)[ids]
  out[vapply(out, is.null, logical(1))] <- list(character())
  stats::setNames(out, ids)
}

# ----------------------------------------------------------------------------- figures
.pd_venn <- function(out, stem, regions, set_ids, formats) {
  .pm_devices(out, stem, formats, 6, 5.5, function() {
    graphics::par(mar = c(1, 1, 3, 1))
    graphics::plot(0, 0, type = "n", xlim = c(-2.2, 2.2), ylim = c(-2.2, 2.2), axes = FALSE, xlab = "", ylab = "", asp = 1,
                   main = sprintf("Venn of %s (protein groups; exclusive regions)", paste(set_ids, collapse = ", ")), cex.main = 0.85)
    k <- length(set_ids)
    centres <- switch(k, cbind(0, 0), cbind(c(-0.6, 0.6), c(0, 0)), cbind(c(-0.6, 0.6, 0), c(0.4, 0.4, -0.6)))
    for (i in seq_len(k)) graphics::symbols(centres[i, 1], centres[i, 2], circles = 1.2, add = TRUE, inches = FALSE, fg = .pm_palette[i], lwd = 2)
    for (i in seq_len(k)) graphics::text(centres[i, 1] * 1.9, centres[i, 2] * 1.9 + if (k == 1) 1.5 else 0, set_ids[i], col = .pm_palette[i], font = 2)
    for (r in seq_len(nrow(regions))) {
      inside <- set_ids %in% strsplit(regions$defining_sets[r], ";")[[1]]
      xy <- colMeans(centres[inside, , drop = FALSE]) - if (sum(inside) < k && k > 1) 0.45 * colMeans(centres[!inside, , drop = FALSE]) else 0
      graphics::text(xy[1], xy[2], regions$n_protein_groups[r])
    }
  })
}

.pd_upset <- function(out, stem, regions, set_ids, formats) {
  regions <- regions[regions$n_protein_groups > 0, , drop = FALSE]
  regions <- regions[order(-regions$n_protein_groups, regions$region_id), , drop = FALSE]
  .pm_devices(out, stem, formats, max(5, 0.45 * nrow(regions) + 3), 5.5, function() {
    graphics::layout(matrix(1:2, ncol = 1), heights = c(2, 1))
    graphics::par(mar = c(0.5, 8, 3, 1))
    n <- max(1L, nrow(regions))
    graphics::barplot(if (nrow(regions)) regions$n_protein_groups else 0, names.arg = NA, col = "grey40", ylab = "protein groups", xlim = c(0, 1.2 * n), width = 1, space = 0.2,
                      main = "UpSet: exclusive region sizes (non-empty regions)", cex.main = 0.9)
    graphics::par(mar = c(1, 8, 0.5, 1))
    graphics::plot(0, 0, type = "n", xlim = c(0, 1.2 * n), ylim = c(0.5, length(set_ids) + 0.5), axes = FALSE, xlab = "", ylab = "")
    graphics::axis(2, at = seq_along(set_ids), labels = set_ids, las = 1, tick = FALSE)
    for (r in seq_len(nrow(regions))) {
      inside <- set_ids %in% strsplit(regions$defining_sets[r], ";")[[1]]
      x <- 0.2 + (r - 1) * 1.2 + 0.5
      graphics::points(rep(x, length(set_ids)), seq_along(set_ids), pch = 19, col = ifelse(inside, "black", "grey85"), cex = 1.4)
      if (sum(inside) > 1) graphics::segments(x, min(which(inside)), x, max(which(inside)), lwd = 2)
    }
  })
}

.pd_concordance_plot <- function(out, stem, source, label, formats) {
  .pm_devices(out, stem, formats, 5.5, 5.5, function() {
    graphics::par(mar = c(4.5, 4.5, 3, 1))
    graphics::plot(source$left_effect, source$right_effect, pch = 19, col = ifelse(source$significant_either, "#B2182B", "grey60"),
                   xlab = "left effect (log2)", ylab = "right effect (log2)", main = label, cex.main = 0.85)
    graphics::abline(h = 0, v = 0, col = "grey80"); graphics::abline(0, 1, lty = 3)
  })
}

# ----------------------------------------------------------------------------- stage handler
post_de_sets_stage <- function(request) .pc_run_stage(request, function(out) {
  p <- request$parameters
  outputs <- list(); warnings <- list(); refusals <- list()
  emit <- function(relative, id, type) outputs[[length(outputs) + 1L]] <<- .pc_output(out, relative, id, type)
  write_tsv <- function(df, relative, id, type) { .pc_write_tsv(df, file.path(out, relative)); emit(relative, id, type) }
  write_json <- function(value, relative, id, type) { .pc_write_json(value, file.path(out, relative)); emit(relative, id, type) }
  refuse <- function(analysis, item, code, reason) {
    refusals[[length(refusals) + 1L]] <<- data.frame(analysis = analysis, item = item, reason_code = code, reason = reason, stringsAsFactors = FALSE)
    warnings[[length(warnings) + 1L]] <<- .pc_warning(request, code, sprintf("%s %s refused: %s", analysis, item, reason), "post_de/sets/refusals.tsv")
  }
  plan <- .pd_verify_inputs(request)
  warnings <- c(warnings, .pd_adaptation_warnings(request))
  for (r in p$refused) refuse(r$analysis, r$item, r$reason_code, r$reason)   # D-59: invalid declarations refused one by one
  zero <- .pd_read_stage_table(request, "stage__dea_zero_null"); treat <- .pd_read_stage_table(request, "stage__dea_treat")
  prim <- .pd_primary(request)
  formats <- intersect(unlist(p$figure_formats), c("png", "pdf", "svg"))
  figures <- list()
  add_fig <- function(x) { for (f in x$files) emit(f$relative_path, f$artifact_id, "Figure"); figures <<- c(figures, x$records) }

  # sets
  defs <- p$definitions
  sets <- list(); universes <- list(); set_rows <- list(); leaf_info <- list()
  for (d in defs) {
    res <- pd_eval_rule(d$rule, zero, treat)
    leaves <- .pd_leaves(d$rule)
    criteria <- unique(vapply(leaves, function(l) l$criterion, ""))
    label <- if ("exploratory_raw_p" %in% criteria) "exploratory_raw_p" else "descriptive"
    sets[[d$id]] <- sort(res$members); universes[[d$id]] <- res$universe
    leaf_info[[d$id]] <- if (length(leaves) == 1L && is.null(d$rule$op)) leaves[[1]] else NULL
    genes <- unique(unlist(.pd_gene_map(prim$features, sets[[d$id]])))
    set_rows[[length(set_rows) + 1L]] <- data.frame(set_id = d$id, rule = as.character(jsonlite::toJSON(d$rule, auto_unbox = TRUE)),
      criteria = paste(criteria, collapse = ";"), thresholds = paste(unique(vapply(leaves, function(l) format(l$threshold), "")), collapse = ";"),
      directions = paste(unique(vapply(leaves, function(l) if (is.null(l$direction)) "any" else l$direction, "")), collapse = ";"),
      family_ids = paste(unique(vapply(leaves, function(l) l$family_id, "")), collapse = ";"),
      endpoints = paste(vapply(leaves, function(l) paste(l$model_id, l$contrast_id, sep = "/"), ""), collapse = ";"),
      null_type = unique(vapply(leaves, function(l) if (identical(l$hypothesis_type, "protein_treat")) "treat" else "zero_null", "")),
      input_matrix = unique(vapply(leaves, function(l) l$input_matrix, "")),
      n_members = length(sets[[d$id]]), n_unique_genes = length(genes), n_tested_universe = length(res$universe), claim_label = label, stringsAsFactors = FALSE)
  }
  set_table <- do.call(rbind, set_rows)
  write_tsv(set_table, "sets.tsv", "post_de_sets", "FeatureSetResult")

  # membership matrix: one row per feature tested in at least one referenced endpoint
  tested_any <- sort(unique(unlist(lapply(defs, function(d) unlist(lapply(.pd_leaves(d$rule), function(l) pd_leaf(l, zero, treat)$tested))))))
  gene_map_all <- .pd_gene_map(prim$features, tested_any)
  membership <- data.frame(feature_id = tested_any, gene_ids = vapply(gene_map_all, paste, "", collapse = ";"), stringsAsFactors = FALSE)
  for (id in names(sets)) membership[[id]] <- tested_any %in% sets[[id]]
  write_tsv(membership, "membership.tsv", "post_de_membership", "FeatureSetResult")

  # exact regions + UpSet (any k) + Venn (k <= 3)
  upset_ids <- if (length(p$upset)) unlist(p$upset) else names(sets)
  regions <- pd_regions(sets[upset_ids], .pd_gene_map(prim$features, unique(unlist(sets[upset_ids]))), all_regions = length(upset_ids) <= 12L)
  write_tsv(regions, "regions.tsv", "post_de_regions", "RegionPartition")
  src <- regions[, c("region_id", "defining_sets", "n_protein_groups", "n_unique_genes")]
  write_tsv(src, file.path("figure_sources", "upset_all.tsv"), "post_de_upset_source", "FigureSource")
  add_fig(.pd_upset(out, "upset_all", regions, upset_ids, formats))
  venn_ids <- unlist(p$venn)
  if (length(venn_ids)) {
    if (length(venn_ids) > 3L) refuse("venn", paste(venn_ids, collapse = ","), "E_VENN_K", sprintf("a circle Venn needs k <= 3 sets; %d were requested (the UpSet is produced)", length(venn_ids)))
    else {
      vr <- pd_regions(sets[venn_ids], .pd_gene_map(prim$features, unique(unlist(sets[venn_ids]))))
      stem <- paste0("venn_", .pd_safe(paste(venn_ids, collapse = "_")))
      write_tsv(vr[, c("region_id", "defining_sets", "n_protein_groups", "n_unique_genes")], file.path("figure_sources", paste0(stem, ".tsv")), paste0(stem, "_source"), "FigureSource")
      add_fig(.pd_venn(out, stem, vr, venn_ids, formats))
    }
  }

  # direction-aware overlap and overlap eligibility for pairs of single-endpoint sets
  single <- names(leaf_info)[!vapply(leaf_info, is.null, logical(1))]
  dir_rows <- list(); elig_rows <- list(); test_rows <- list()
  contrasts <- stats::setNames(plan$contrasts, vapply(plan$contrasts, function(c) c$contrast_id, ""))
  if (length(single) >= 2L) for (pair in utils::combn(single, 2L, simplify = FALSE)) {
    la <- leaf_info[[pair[1]]]; lb <- leaf_info[[pair[2]]]
    ea <- pd_leaf(la, zero, treat)$effect; eb <- pd_leaf(lb, zero, treat)$effect
    ov <- pd_direction_overlap(sets[[pair[1]]], sets[[pair[2]]], ea, eb)
    dir_rows[[length(dir_rows) + 1L]] <- data.frame(set_a = pair[1], set_b = pair[2], n_concordant_up = length(ov$concordant_up), n_concordant_down = length(ov$concordant_down),
      n_discordant = length(ov$discordant), n_only_a = length(ov$only_a), n_only_b = length(ov$only_b),
      concordant_up = paste(ov$concordant_up, collapse = ";"), concordant_down = paste(ov$concordant_down, collapse = ";"), discordant = paste(ov$discordant, collapse = ";"),
      only_a = paste(ov$only_a, collapse = ";"), only_b = paste(ov$only_b, collapse = ";"), claim_label = "descriptive", stringsAsFactors = FALSE)
    ua <- pd_contrast_units(contrasts[[la$contrast_id]], prim$obs, p$group_column, p$subject_column)
    ub <- pd_contrast_units(contrasts[[lb$contrast_id]], prim$obs, p$group_column, p$subject_column)
    shared <- unique(c(intersect(ua$observations, ub$observations), intersect(ua$units, ub$units), intersect(ua$subjects, ub$subjects), intersect(ua$biological_units, ub$biological_units)))
    eligible <- length(shared) == 0L && identical(la$hypothesis_type, "protein_zero_null") && identical(lb$hypothesis_type, "protein_zero_null")
    reason <- if (length(shared)) "shared_biological_units" else if (!eligible) "not_zero_null_sets" else "disjoint_biological_units"
    elig_rows[[length(elig_rows) + 1L]] <- data.frame(set_a = pair[1], set_b = pair[2], contrast_a = la$contrast_id, contrast_b = lb$contrast_id,
      n_shared_units = length(shared), overlap_p_eligible = eligible && isTRUE(p$overlap_test), reason = reason, stringsAsFactors = FALSE)
    if (eligible && isTRUE(p$overlap_test)) {
      universe <- intersect(pd_leaf(la, zero, treat)$tested, pd_leaf(lb, zero, treat)$tested)
      o <- pd_overlap_p(sets[[pair[1]]], sets[[pair[2]]], universe)
      test_rows[[length(test_rows) + 1L]] <- data.frame(set_a = pair[1], set_b = pair[2], universe = "features tested in both contrasts", universe_n = o$universe_n,
        n_a = o$n_a, n_b = o$n_b, overlap = o$overlap, expected = o$expected, p_value = o$p_value, test = "hypergeometric upper tail", family = "descriptive_overlap",
        claim_label = "descriptive", stringsAsFactors = FALSE)
    }
  }
  dir_df <- if (length(dir_rows)) do.call(rbind, dir_rows) else data.frame(set_a = character(), set_b = character(), claim_label = character())
  write_tsv(dir_df, "direction_overlap.tsv", "post_de_direction_overlap", "FeatureSetResult")
  elig_df <- if (length(elig_rows)) do.call(rbind, elig_rows) else data.frame(set_a = character(), set_b = character(), reason = character())
  write_tsv(elig_df, "overlap_eligibility.tsv", "post_de_overlap_eligibility", "PostDeEligibility")
  if (length(test_rows)) {
    tests <- do.call(rbind, test_rows); tests$q_value <- stats::p.adjust(tests$p_value, "BH")
    write_tsv(tests, "overlap_test.tsv", "post_de_overlap_test", "FeatureSetResult")
  }
  for (i in seq_len(nrow(elig_df))) if (isTRUE(p$overlap_test) && !elig_df$overlap_p_eligible[i])
    refuse("overlap_test", paste(elig_df$set_a[i], elig_df$set_b[i], sep = " x "), "W_OVERLAP_P_INELIGIBLE", elig_df$reason[i])

  # effect concordance (descriptive; no P column)
  conc_rows <- list()
  for (pair in p$concordance_pairs) {
    lx <- pd_leaf(pair$left, zero, treat); ly <- pd_leaf(pair$right, zero, treat)
    sig <- union(lx$members, ly$members)
    cc <- pd_concordance(lx$effect[lx$tested], ly$effect[ly$tested], sig)
    conc_rows[[length(conc_rows) + 1L]] <- data.frame(pair_id = pair$id, left = paste(pair$left$model_id, pair$left$contrast_id, sep = "/"),
      right = paste(pair$right$model_id, pair$right$contrast_id, sep = "/"), n_common_tested = cc$n_common, pearson = cc$pearson, spearman = cc$spearman,
      origin_slope = cc$origin_slope, significance_criterion = sprintf("%s < %s", pair$left$criterion, format(pair$left$threshold)),
      n_significant_either = cc$n_significant_either, n_sign_agree = cc$n_sign_agree, sign_agreement = cc$sign_agreement, claim_label = "descriptive", stringsAsFactors = FALSE)
    common <- intersect(lx$tested, ly$tested)
    source <- data.frame(feature_id = common, left_effect = lx$effect[common], right_effect = ly$effect[common], significant_either = common %in% sig, stringsAsFactors = FALSE)
    stem <- paste0("concordance_", .pd_safe(pair$id))
    write_tsv(source, file.path("figure_sources", paste0(stem, ".tsv")), paste0(stem, "_source"), "FigureSource")
    add_fig(.pd_concordance_plot(out, stem, source, sprintf("%s: Pearson %.2f, origin slope %.2f (descriptive)", pair$id, cc$pearson, cc$origin_slope), formats))
  }
  conc_df <- if (length(conc_rows)) do.call(rbind, conc_rows) else data.frame(pair_id = character(), claim_label = character())
  write_tsv(conc_df, "concordance.tsv", "post_de_concordance", "ConcordanceResult")

  refusals_df <- .pd_refusal_frame(refusals)
  write_tsv(refusals_df, "refusals.tsv", "post_de_sets_refusals", "PostDeEligibility")
  for (table in list(set_table, regions, dir_df, conc_df)) .pd_check_claims(table)
  write_json(list(module = "sets", state = "COMPLETED", reason_code = NULL, claim_label = if (any(set_table$claim_label == "exploratory_raw_p")) "exploratory_raw_p" else "descriptive",
                  eligibility = p$eligibility, refusals = refusals_df, input_hashes = .pd_input_hashes(request), figures = figures,
                  rule = "finite set grammar (union, intersect, difference, complement_within_tested) over completed families; q < threshold; overlap P only for unit-disjoint contrasts"),
             "eligibility.json", "post_de_sets_eligibility", "PostDeEligibility")
  list(outputs = outputs, warnings = warnings, message = sprintf("post-DE sets: %d set(s), %d region(s)", length(sets), nrow(regions)))
})
