# Post-DE co-abundance modules and protein interaction networks (packet R14e, SM39-SM40, FR-161 to FR-165).
#
# Co-abundance (opt-in, eligibility-gated): complete genuinely observed features of at least the declared number of
# biological units.  Rules (native base R, D-44): wgcna_signed = signed adjacency ((1 + r)/2)^beta with beta the
# smallest candidate whose scale-free-topology fit reaches R^2 >= 0.8, signed TOM dissimilarity, average-linkage
# clustering, dynamicTreeCut; hclust_correlation = 1 - r dissimilarity, average linkage, dynamicTreeCut.  Eigengene =
# first principal component of the module's scaled data (sign aligned with the module mean).  Stability: unit-level
# bootstrap.  Module-trait tests: limma on eigengenes in a module-level family.  Networks: hashed offline snapshots
# only; connectivity of a declared set against degree-preserving random sets from the measured, mapped universe.

PD_SCALE_FREE_TARGET <- 0.8
PD_SOFT_THRESHOLDS <- c(1:10, seq(12L, 30L, by = 2L))

# Scale-free topology fit index of a connectivity vector, as defined for WGCNA's pickSoftThreshold (Zhang and Horvath
# 2005): k is cut into `bins` equal-width bins, log10(mean k per bin) is regressed on log10(bin frequency + 1e-9)
# (empty bins enter at their midpoint with frequency 0), and the index is -sign(slope) * R^2.
pd_scale_free_fit <- function(k, bins = 10L) {
  if (length(unique(k)) < 2L) return(list(r2 = NA_real_, slope = NA_real_))
  discretized <- cut(k, bins)
  dk <- tapply(k, discretized, mean)
  p_dk <- as.vector(tapply(k, discretized, length) / length(k))
  mids <- graphics::hist(k, breaks = seq(min(k), max(k), length.out = bins + 1L), plot = FALSE, right = TRUE)$mids
  dk <- ifelse(is.na(dk) | dk == 0, mids, dk); p_dk <- ifelse(is.na(p_dk), 0, p_dk)
  x <- as.vector(log10(dk)); y <- log10(p_dk + 1e-9)
  fit <- stats::lm(y ~ x)
  slope <- unname(stats::coef(fit)[2]); r2 <- summary(fit)$r.squared
  list(r2 = -sign(slope) * r2, slope = slope)
}

pd_signed_adjacency <- function(C, beta) { A <- ((1 + C) / 2)^beta; diag(A) <- 0; A }

pd_signed_tom <- function(A) {
  k <- rowSums(A); L <- A %*% A
  tom <- (L + A) / (outer(k, k, pmin) + 1 - A)
  diag(tom) <- 1
  tom
}

# Smallest candidate reaching the scale-free fit target; when none does, the published signed-network default by
# sample size (n < 20: 18, < 30: 16, < 40: 14, otherwise 12), recorded with a warning (D-49).
pd_choose_beta <- function(C, n_units, candidates = PD_SOFT_THRESHOLDS, target = PD_SCALE_FREE_TARGET) {
  fits <- do.call(rbind, lapply(candidates, function(b) { f <- pd_scale_free_fit(rowSums(pd_signed_adjacency(C, b))); data.frame(beta = b, scale_free_r2 = f$r2, slope = f$slope) }))
  ok <- which(is.finite(fits$scale_free_r2) & fits$scale_free_r2 >= target)
  chosen <- if (length(ok)) fits$beta[ok[1]] else if (n_units < 20) 18L else if (n_units < 30) 16L else if (n_units < 40) 14L else 12L
  list(beta = chosen, table = fits, reached = length(ok) > 0L)
}

# Module labels (0 = unassigned) for a units x features matrix under a rule.
pd_modules <- function(X, rule, beta, min_size, deep_split, cut_height) {
  C <- stats::cor(X)
  D <- if (identical(rule, "wgcna_signed")) 1 - pd_signed_tom(pd_signed_adjacency(C, beta)) else 1 - C
  D[D < 0] <- 0
  tree <- stats::hclust(stats::as.dist(D), method = "average")
  labels <- dynamicTreeCut::cutreeDynamic(tree, distM = D, deepSplit = deep_split, minClusterSize = min_size, pamStage = FALSE, cutHeight = cut_height, verbose = 0)   # no PAM stage: unassigned features stay unassigned
  stats::setNames(as.integer(labels), colnames(X))
}

pd_eigengenes <- function(X, labels) {
  mods <- sort(setdiff(unique(labels), 0L))
  out <- list(); variance <- c()
  for (m in mods) {
    Z <- scale(X[, labels == m, drop = FALSE])
    pc <- stats::prcomp(Z, center = FALSE, scale. = FALSE)
    e <- pc$x[, 1]
    if (stats::cor(e, rowMeans(Z)) < 0) e <- -e
    out[[paste0("M", m)]] <- e; variance[paste0("M", m)] <- pc$sdev[1]^2 / sum(pc$sdev^2)
  }
  list(eigengenes = do.call(cbind, out), variance = variance)
}

pd_adjusted_rand <- function(a, b) {
  tab <- table(a, b); n <- sum(tab); c2 <- function(x) sum(x * (x - 1) / 2)
  s <- c2(tab); sa <- c2(rowSums(tab)); sb <- c2(colSums(tab)); e <- sa * sb / (n * (n - 1) / 2)
  (s - e) / (0.5 * (sa + sb) - e)
}

# Undirected edges: each unordered gene pair once (highest score kept), self-loops and edges below min_score dropped.
pd_undirected_edges <- function(edges, min_score) {
  e <- edges[!is.na(edges$score) & edges$score >= min_score & edges$gene_a != edges$gene_b, , drop = FALSE]
  if (!nrow(e)) return(data.frame(gene_a = character(), gene_b = character(), score = numeric(), stringsAsFactors = FALSE))
  a <- pmin(e$gene_a, e$gene_b); b <- pmax(e$gene_a, e$gene_b)
  e <- data.frame(gene_a = a, gene_b = b, score = e$score, stringsAsFactors = FALSE)
  e <- e[order(e$gene_a, e$gene_b, -e$score), , drop = FALSE]
  e[!duplicated(e[, c("gene_a", "gene_b")]), , drop = FALSE]
}

# Degree-preserving null: each set gene is replaced by an unused universe gene of the same degree bin.
pd_degree_bins <- function(degree, min_bin = 5L) {
  d <- sort(unique(degree)); bin <- stats::setNames(integer(length(degree)), names(degree)); current <- 1L; members <- 0L
  for (v in d) { bin[degree == v] <- current; members <- members + sum(degree == v); if (members >= min_bin) { current <- current + 1L; members <- 0L } }
  if (members > 0L && current > 1L) bin[bin == current] <- current - 1L
  bin
}

pd_connectivity_null <- function(set_genes, universe, edges_in, degree, draws, seed) {
  bins <- pd_degree_bins(degree[universe])
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  count_edges <- function(g) sum(edges_in$a %in% g & edges_in$b %in% g)
  null <- vapply(seq_len(draws), function(i) {
    chosen <- character()
    for (g in set_genes) {
      width <- 0L
      repeat { pool <- setdiff(universe[abs(bins[universe] - bins[[g]]) <= width], chosen); if (length(pool) || width > max(bins)) break; width <- width + 1L }   # adjacent bins only when a bin is exhausted
      chosen <- c(chosen, if (length(pool) == 1L) pool else sample(pool, 1L))
    }
    if (!all(chosen %in% universe)) stop("E_NETWORK_UNIVERSE: a null draw left the measured mapped universe", call. = FALSE)
    count_edges(chosen) }, numeric(1))
  list(observed = count_edges(set_genes), null = null)
}

# Unit-level bootstrap stability: biological units (never features) are resampled with replacement from one
# L'Ecuyer stream; every draw re-runs module construction on all features with fixed parameters and records, per
# reference module, the best Jaccard overlap with any module of the draw.  A draw whose construction fails is NA
# and counted in `failed`.
pd_module_stability <- function(X, labels, units, rule, beta, min_size, deep_split, cut_height, draws, seed) {
  mods <- paste0("M", sort(setdiff(unique(labels), 0L)))
  jac <- matrix(NA_real_, draws, length(mods), dimnames = list(NULL, mods))
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  u <- unique(units); picks <- vector("list", draws); failed <- 0L
  for (b in seq_len(draws)) {
    pick <- sample(u, length(u), replace = TRUE); picks[[b]] <- pick
    rows <- unlist(lapply(pick, function(x) which(units == x)))
    lb <- tryCatch(pd_modules(X[rows, , drop = FALSE], rule, beta, min_size, deep_split, cut_height), error = function(e) NULL)
    if (is.null(lb)) { failed <- failed + 1L; next }
    for (m in mods) {
      a <- names(labels)[paste0("M", labels) == m]; best <- 0
      for (bm in setdiff(unique(lb), 0L)) { g <- names(lb)[lb == bm]; best <- max(best, length(intersect(a, g)) / length(union(a, g))) }
      jac[b, m] <- best
    }
  }
  list(jaccard = jac, picks = picks, failed = failed)
}

post_de_networks_stage <- function(request) .pc_run_stage(request, function(out) {
  p <- request$parameters
  outputs <- list(); warnings <- list(); refusals <- list(); figures <- list()
  emit <- function(relative, id, type) outputs[[length(outputs) + 1L]] <<- .pc_output(out, relative, id, type)
  write_tsv <- function(df, relative, id, type) { .pc_write_tsv(df, file.path(out, relative)); emit(relative, id, type) }
  write_json <- function(value, relative, id, type) { .pc_write_json(value, file.path(out, relative)); emit(relative, id, type) }
  add_fig <- function(x) { for (f in x$files) emit(f$relative_path, f$artifact_id, "Figure"); figures <<- c(figures, x$records) }
  plan <- .pd_verify_inputs(request)
  prim <- .pd_primary(request); obs <- prim$obs; fm <- unlist(p$figure_formats)
  for (s in p$eligibility$subanalyses) refusals[[length(refusals) + 1L]] <- data.frame(analysis = s$analysis, item = s$item, reason_code = s$reason_code, reason = s$reason, stringsAsFactors = FALSE)
  co <- p$coabundance; summary <- list()
  if (!is.null(co) && identical(co$state, "ELIGIBLE")) {
    if (!requireNamespace("dynamicTreeCut", quietly = TRUE)) stop("E_ENGINE_NOT_AVAILABLE: R package dynamicTreeCut is required for co-abundance module cutting", call. = FALSE)
    feats <- unlist(co$features)
    X <- t(prim$Y[feats, , drop = FALSE])
    if (anyNA(X)) stop("E_INTEGRITY: co-abundance features must be observed in every unit (no imputation)", call. = FALSE)
    beta <- NA_integer_; fit_table <- NULL
    if (identical(co$rule, "wgcna_signed")) {
      if (identical(co$soft_threshold, "auto")) { ch <- pd_choose_beta(stats::cor(X), nrow(X)); beta <- ch$beta; fit_table <- ch$table
        if (!ch$reached) warnings[[length(warnings) + 1L]] <- .pc_warning(request, "W_SCALE_FREE_FIT", sprintf("no soft threshold reached a scale-free fit R2 >= %.1f; the signed-network default beta = %d for %d units is used", PD_SCALE_FREE_TARGET, beta, nrow(X)), "post_de/networks/soft_threshold.tsv") }
      else beta <- as.integer(co$soft_threshold)
    }
    labels <- pd_modules(X, co$rule, beta, as.integer(co$min_module_size), as.integer(co$deep_split), as.numeric(co$cut_height))
    eg <- pd_eigengenes(X, labels)
    modules <- data.frame(feature_id = names(labels), module = ifelse(labels == 0L, "unassigned", paste0("M", labels)), rule = co$rule, soft_threshold = beta, cut_height = as.numeric(co$cut_height),
                          min_module_size = as.integer(co$min_module_size), claim_label = "module_level", stringsAsFactors = FALSE)
    write_tsv(modules, "modules.tsv", "post_de_modules", "CoabundanceModules")
    if (!is.null(fit_table)) write_tsv(cbind(fit_table, claim_label = "descriptive"), "soft_threshold.tsv", "post_de_soft_threshold", "CoabundanceModules")
    E <- eg$eigengenes
    if (is.null(E) || !ncol(E)) {
      warnings[[length(warnings) + 1L]] <- .pc_warning(request, "W_NO_MODULES", "no co-abundance module reached the minimum module size", "post_de/networks/modules.tsv")
    } else {
      write_tsv(data.frame(observation_id = rownames(X), E, check.names = FALSE, claim_label = "module_level"), "eigengenes.tsv", "post_de_eigengenes", "CoabundanceModules")
      # unit-level bootstrap stability (units resampled with replacement; beta and parameters fixed)
      units <- .pd_units(obs[match(rownames(X), obs$observation_id), , drop = FALSE], p$subject_column)
      boot <- pd_module_stability(X, labels, units, co$rule, beta, as.integer(co$min_module_size), as.integer(co$deep_split), as.numeric(co$cut_height), as.integer(co$bootstrap), as.integer(co$seed))
      jac <- boot$jaccard; mods <- colnames(E)
      stability <- data.frame(module = mods, n_features = as.integer(table(labels[labels != 0L])[sub("^M", "", mods)]), variance_explained = unname(eg$variance[mods]),
                              bootstrap = as.integer(co$bootstrap), seed = as.integer(co$seed), resampling = "biological units with replacement",
                              failed_draws = boot$failed, mean_jaccard = unname(colMeans(jac[, mods, drop = FALSE], na.rm = TRUE)), sd_jaccard = unname(apply(jac[, mods, drop = FALSE], 2L, stats::sd, na.rm = TRUE)),
                              claim_label = "descriptive", stringsAsFactors = FALSE)
      write_tsv(stability, "module_stability.tsv", "post_de_module_stability", "CoabundanceModules")
      # module-trait association: limma on eigengenes (module-level family per trait)
      trait_rows <- list()
      Em <- t(E); colnames(Em) <- rownames(X)
      observedE <- matrix(TRUE, nrow(Em), ncol(Em), dimnames = dimnames(Em))
      settings <- list(group_column = p$group_column, subject_column = p$subject_column, blocking_mode = p$blocking_mode, consensus_correlation = p$consensus_correlation,
                       trend = FALSE, robust = isTRUE(p$robust), ci_level = p$ci_level)   # eigengenes are centred: their means are identically 0, so a mean-variance trend is undefined (D-49)
      for (tr in co$trait_designs) {
        Xd <- .pd_design_from(tr$design); ids <- rownames(Xd)
        fam <- list(family_id = tr$family_id, hypothesis_type = "protein_zero_null", role = "secondary", adjustment = "BH", q_cutoff = 0.05, dependence_assumption = "BH",
                    members = lapply(tr$contrasts, function(c) list(model_id = tr$model_id, contrast_id = c$contrast_id)))
        res <- pd_fit_model(tr$model_id, Em[, ids, drop = FALSE], observedE[, ids, drop = FALSE], obs[match(ids, obs$observation_id), , drop = FALSE], Xd, tr$contrasts, settings, list(fam), list(run_id = request$run_id, plan_hash = request$plan_hash), p$coverage)
        r <- res$rows
        trait_rows[[length(trait_rows) + 1L]] <- data.frame(trait = tr$trait, model_id = tr$model_id, contrast_id = r$contrast_id, module = r$feature_id, eligibility = r$eligibility,
          effect = r$effect, statistic = r$statistic, p_value = r$p_value, q_value = r$q_value, family_id = r$family_id, inference = "module-level (eigengene), not individual proteins", moderation = sprintf("limma eBayes trend=FALSE (eigengene means are 0), robust=%s", if (isTRUE(p$robust)) "TRUE" else "FALSE"),
          claim_label = "module_level", stringsAsFactors = FALSE)
      }
      if (length(trait_rows)) write_tsv(do.call(rbind, trait_rows), "module_trait.tsv", "post_de_module_trait", "ModuleTraitResult")
      src <- data.frame(observation_id = rownames(X), group = obs[[p$group_column]][match(rownames(X), obs$observation_id)], E, check.names = FALSE, claim_label = "module_level")
      write_tsv(src, file.path("figure_sources", "eigengenes.tsv"), "eigengenes_source", "FigureSource")
      add_fig(.pm_devices(out, "eigengenes", fm, 6.5, 4.5, function() { graphics::par(mar = c(4.5, 4.5, 3, 1))
        graphics::boxplot(stats::as.formula(paste(sprintf("`%s`", colnames(E)[1]), "~ group")), data = src, col = "grey80", main = sprintf("%s eigengene by group (module-level)", colnames(E)[1]), cex.main = 0.85) }))
    }
    summary$coabundance <- list(n_units = nrow(X), n_features = ncol(X), rule = co$rule, soft_threshold = beta, n_modules = if (is.null(E)) 0L else ncol(E))
  } else if (!is.null(co) && !is.null(co$state)) refusals[[length(refusals) + 1L]] <- data.frame(analysis = "coabundance", item = "modules", reason_code = co$reason_code, reason = co$reason, stringsAsFactors = FALSE)

  # interaction networks from the hashed offline snapshot
  pp <- p$ppi
  if (!is.null(pp) && identical(pp$state, "ELIGIBLE")) {
    edges <- .pc_read_tsv(.pc_find_input(request, paste0("resfile_", pp$snapshot_id, "__", pp$edge_file)))
    edges$score <- suppressWarnings(as.numeric(edges$score))
    if (anyNA(edges$score)) stop("E_RESOURCE_SCHEMA: the interaction snapshot has non-numeric scores", call. = FALSE)
    edges <- pd_undirected_edges(edges, as.numeric(pp$min_score))
    mapping <- .pc_read_tsv(.pc_find_input(request, paste0("resfile_", pp$mapping_resource_id, "__", pp$mapping_file)))
    mapped <- map_features_to_genes(prim$features, mapping, p$organism_taxonomy_id)
    measured <- unique(stats::na.omit(mapped$gene_id[mapped$mapping_state == "mapped"]))
    universe <- sort(intersect(measured, unique(c(edges$gene_a, edges$gene_b))))
    edges_in <- data.frame(a = edges$gene_a, b = edges$gene_b, stringsAsFactors = FALSE)
    edges_in <- edges_in[edges_in$a %in% universe & edges_in$b %in% universe, , drop = FALSE]
    degree <- stats::setNames(vapply(universe, function(g) sum(edges_in$a == g | edges_in$b == g), numeric(1)), universe)
    conn <- list(); subnet <- list(); hubs <- list()
    for (set in pp$sets) {
      genes <- unique(stats::na.omit(mapped$gene_id[mapped$feature_id %in% unlist(set$members) & mapped$mapping_state == "mapped"]))
      genes <- intersect(genes, universe)
      if (length(genes) < 2L) { refusals[[length(refusals) + 1L]] <- data.frame(analysis = "ppi", item = set$id, reason_code = "E_NETWORK_SET_SIZE", reason = "fewer than two set genes in the measured mapped network universe", stringsAsFactors = FALSE); next }
      res <- pd_connectivity_null(genes, universe, edges_in, degree, as.integer(pp$null_draws), as.integer(pp$seed))
      k <- sum(res$null >= res$observed)
      conn[[length(conn) + 1L]] <- data.frame(set_id = set$id, n_set_features = length(unlist(set$members)), n_set_genes_in_universe = length(genes), universe = "measured mapped genes in the snapshot",
        universe_size = length(universe), observed_edges = res$observed, null_mean_edges = mean(res$null), null_draws = length(res$null), seed = as.integer(pp$seed), k = k,
        p_value = (k + 1) / (length(res$null) + 1), null = "degree-preserving random sets from the measured mapped universe", snapshot_release = pp$release, snapshot_sha256 = pp$snapshot_sha256,
        claim_label = "descriptive", stringsAsFactors = FALSE)
      inside <- edges[edges$gene_a %in% genes & edges$gene_b %in% genes, , drop = FALSE]
      if (nrow(inside)) subnet[[length(subnet) + 1L]] <- data.frame(set_id = set$id, gene_a = inside$gene_a, gene_b = inside$gene_b, score = inside$score, claim_label = "descriptive", stringsAsFactors = FALSE)
      hubs[[length(hubs) + 1L]] <- data.frame(set_id = set$id, gene_id = genes, degree_in_universe = unname(degree[genes]),
        degree_in_set = vapply(genes, function(g) sum((inside$gene_a == g) | (inside$gene_b == g)), numeric(1)), claim_label = "descriptive", stringsAsFactors = FALSE)
    }
    cdf <- .pd_rows(conn, c("set_id", "claim_label"))
    if (nrow(cdf) > 1L) cdf$q_value <- stats::p.adjust(cdf$p_value, "BH")
    write_tsv(cdf, "connectivity.tsv", "post_de_connectivity", "NetworkResult")
    sdf <- .pd_rows(subnet, c("set_id", "gene_a", "gene_b", "claim_label")); write_tsv(sdf, "ppi_subnetwork.tsv", "post_de_ppi_subnetwork", "NetworkResult")
    hdf <- .pd_rows(hubs, c("set_id", "gene_id", "claim_label")); write_tsv(hdf, "hub_degree.tsv", "post_de_hub_degree", "NetworkResult")
    if (nrow(sdf)) {
      first <- sdf[sdf$set_id == sdf$set_id[1], , drop = FALSE]; nodes <- sort(unique(c(first$gene_a, first$gene_b))); ang <- 2 * pi * seq_along(nodes) / length(nodes)
      lay <- data.frame(set_id = first$set_id[1], gene_id = nodes, x = cos(ang), y = sin(ang), claim_label = "descriptive")
      write_tsv(lay, file.path("figure_sources", "ppi_subnetwork.tsv"), "ppi_subnetwork_source", "FigureSource")
      add_fig(.pm_devices(out, "ppi_subnetwork", fm, 5.5, 5.5, function() { graphics::par(mar = c(1, 1, 3, 1))
        graphics::plot(lay$x, lay$y, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1, xlim = c(-1.3, 1.3), ylim = c(-1.3, 1.3), main = sprintf("Induced interactions of %s (descriptive)", lay$set_id[1]), cex.main = 0.85)
        for (i in seq_len(nrow(first))) graphics::segments(lay$x[lay$gene_id == first$gene_a[i]], lay$y[lay$gene_id == first$gene_a[i]], lay$x[lay$gene_id == first$gene_b[i]], lay$y[lay$gene_id == first$gene_b[i]], col = "grey50")
        graphics::points(lay$x, lay$y, pch = 21, bg = "#0072B2", cex = 1.5); graphics::text(lay$x * 1.15, lay$y * 1.15, lay$gene_id, cex = 0.6) }))
    }
    summary$ppi <- list(snapshot_id = pp$snapshot_id, source = pp$source, release = pp$release, species = pp$species, score_type = pp$score_type, snapshot_sha256 = pp$snapshot_sha256,
                        universe_size = length(universe), n_edges_in_universe = nrow(edges_in), min_score = pp$min_score)
  } else if (!is.null(pp) && !is.null(pp$state)) refusals[[length(refusals) + 1L]] <- data.frame(analysis = "ppi", item = "connectivity", reason_code = pp$reason_code, reason = pp$reason, stringsAsFactors = FALSE)
  refusals_df <- .pd_refusal_frame(refusals)
  write_tsv(cbind(refusals_df, claim_label = rep("descriptive", nrow(refusals_df))), "refusals.tsv", "post_de_networks_refusals", "PostDeEligibility")
  for (r in seq_len(nrow(refusals_df))) warnings[[length(warnings) + 1L]] <- .pc_warning(request, refusals_df$reason_code[r], sprintf("networks %s %s: %s", refusals_df$analysis[r], refusals_df$item[r], refusals_df$reason[r]), "post_de/networks/refusals.tsv")
  write_json(list(module = "networks", state = "COMPLETED", reason_code = NULL, claim_label = "module_level", summary = summary, eligibility = p$eligibility, refusals = refusals_df,
                  input_hashes = .pd_input_hashes(request), figures = figures,
                  rule = "co-abundance only when declared with >= min units of complete observed features; modules from the declared rule; eigengene module-trait tests in a module-level family; PPI connectivity against degree-preserving random sets from the measured mapped universe of a hashed offline snapshot",
                  limitations = c("Hub and degree results are descriptive; network structure is not evidence of mechanism or causation.",
                                  "Module-trait results are module-level (eigengene) inference, not inference on individual proteins.")),
             "eligibility.json", "post_de_networks_eligibility", "PostDeEligibility")
  list(outputs = outputs, warnings = warnings, message = "post-DE networks completed")
})
