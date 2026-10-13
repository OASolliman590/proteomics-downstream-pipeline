# Offline GMT pathway libraries (packet R15b: FR-174 to FR-178, SM44, SM17; ADR 0011).
#
# Reads user-supplied local GMT files only; the caller verifies every file SHA-256 before the stage runs and
# the stage verifies it again through the standard input check.  ORA uses the measured, mapped and eligible
# universe (SM17) and the existing exact hypergeometric adapter (ora_test, pathway_ora.R).  Pre-ranked
# enrichment ranks the Welch t of R15a, tie-free under the SM15 stable-ID rule; it enumerates every subset when
# that is feasible and otherwise falls back to the existing exploratory fgsea adapter (run_fgsea).  Each
# library x method x contrast is its own multiplicity family (FR-178).  Library content is never copied into an
# output: outputs hold the manifest records and the result tables only (FR-175).

# GMT: name <TAB> description <TAB> gene ... (one set per line).
gmt_parse <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- lines[nzchar(trimws(lines))]
  if (!length(lines)) stop("E_LIBRARY_FORMAT: the GMT file has no sets", call. = FALSE)
  sets <- lapply(lines, function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 3L) stop("E_LIBRARY_FORMAT: a GMT row needs a name, a description and at least one gene", call. = FALSE)
    genes <- unique(fields[-(1:2)]); genes <- genes[nzchar(genes)]
    if (!length(genes)) stop("E_LIBRARY_FORMAT: a GMT set has no genes", call. = FALSE)
    genes
  })
  names_ <- vapply(strsplit(lines, "\t", fixed = TRUE), function(f) f[1], "")
  if (anyDuplicated(names_)) stop("E_LIBRARY_FORMAT: duplicated set names", call. = FALSE)
  stats::setNames(sets, names_)
}

# Licence rule (FR-175): KEGG and MSigDB libraries are never redistributable; a library that claims redistribution
# must be the synthetic, marked fixture.  Returns the typed refusal reason or NULL.
library_redistribution_refusal <- function(licence_id, redistributable, synthetic) {
  restricted <- grepl("^(KEGG|MSIGDB)", toupper(licence_id))
  if (restricted && isTRUE(redistributable)) return("restricted licence declared redistributable")
  if (isTRUE(redistributable) && !isTRUE(synthetic)) return("non-synthetic library declared redistributable")
  NULL
}

check_artifact_redistribution <- function(artifact_sha256, libraries) {
  for (library in libraries) {
    hit <- library$sha256 %in% artifact_sha256
    refusal <- library_redistribution_refusal(library$licence_id, library$redistributable, library$synthetic)
    if (!is.null(refusal) && (hit || isTRUE(library$redistributable))) stop(sprintf("E_LIBRARY_REDISTRIBUTION: %s (library %s); nothing is written", refusal, library$id), call. = FALSE)
    if (hit && grepl("^(KEGG|MSIGDB)", toupper(library$licence_id))) stop(sprintf("E_LIBRARY_REDISTRIBUTION: a restricted library (%s) is part of the artifact", library$id), call. = FALSE)
  }
  invisible(TRUE)
}

# SM17: the ORA background is the measured, mapped and eligible universe; any other background is refused.
ora_background_check <- function(background, measured_universe) {
  if (!identical(background$kind, "measured")) stop(sprintf("E_ORA_UNIVERSE: background %s is not the measured, mapped and eligible universe", background$kind), call. = FALSE)
  if (!identical(as.integer(background$size), as.integer(length(unique(measured_universe))))) stop(sprintf("E_ORA_UNIVERSE: declared universe size %s differs from the measured universe of %d genes", background$size, length(unique(measured_universe))), call. = FALSE)
  invisible(TRUE)
}

# Family discipline (FR-178): a family key that spans two libraries is a pooled family and is refused.
enrichment_family_adjust <- function(rows, family_key) {
  key <- family_key(rows)
  for (k in unique(key)) {
    if (length(unique(rows$library_id[key == k])) > 1L) stop(sprintf("E_ENRICHMENT_FAMILY_POOLED: family %s spans more than one library", k), call. = FALSE)
  }
  rows$family_id <- key
  rows$q_value <- NA_real_
  for (k in unique(key)) {
    sel <- key == k & is.finite(rows$p_value)
    if (any(sel)) rows$q_value[sel] <- adjust_pvalues(rows$p_value[sel], "BH")
  }
  rows
}

gene_set_family_key <- function(rows) paste(rows$library_id, rows$method, rows$contrast_id, sep = "|")

gene_set_libraries_stage <- function(request) .pc_run_stage(request, function(out) {
  p <- request$parameters
  outputs <- list(); warnings <- list()
  emit <- function(relative, id, type) outputs[[length(outputs) + 1L]] <<- .pc_output(out, relative, id, type)
  write_tsv <- function(df, relative, id, type) { .pc_write_tsv(df, file.path(out, relative)); emit(relative, id, type) }
  write_json <- function(value, relative, id, type) { .pc_write_json(value, file.path(out, relative)); emit(relative, id, type) }
  welch <- .pc_read_tsv(.pc_find_input(request, "welch_results"))
  welch <- welch[welch$eligibility == "tested", , drop = FALSE]
  welch$effect <- as.numeric(welch$effect); welch$statistic <- as.numeric(welch$statistic); welch$q_value <- as.numeric(welch$q_value)
  ora_threshold <- if (is.null(p$ora_threshold)) 0.05 else as.numeric(p$ora_threshold)
  tie_break <- if (is.null(p$tie_break)) "none" else p$tie_break
  libraries <- p$libraries
  for (library in libraries) {
    sets <- gmt_parse(.pc_find_input(request, paste0("library_", library$id)))
    rows <- list()
    for (contrast_id in unique(welch$contrast_id)) {
      w <- welch[welch$contrast_id == contrast_id, , drop = FALSE]
      universe <- sort(unique(w$gene_id[nzchar(w$gene_id)]), method = "radix")
      ora_background_check(p$background, universe)
      eligible_sets <- lapply(sets, function(genes) intersect(genes, universe))
      # ORA: measured background, foreground = Welch family q within threshold and the stated direction (SM17)
      for (direction in c("up", "down")) {
        sel <- w$q_value <= ora_threshold & (if (direction == "up") w$effect > 0 else w$effect < 0)
        foreground <- unique(w$gene_id[sel & nzchar(w$gene_id)])
        o <- ora_test(universe, foreground, sets, as.integer(p$min_size), as.integer(p$max_size))
        if (!is.null(o)) rows[[length(rows) + 1L]] <- data.frame(library_id = library$id, method = paste0("ora_", direction), contrast_id = contrast_id,
                                                               set_id = o$set_id, N = o$N, K = o$K, n = o$n, k = o$k, statistic = NA_real_, es = NA_real_,
                                                               p_value = o$p_value, exploratory = FALSE, estimator = "hypergeometric_measured_background",
                                                               stringsAsFactors = FALSE)
      }
      # Pre-ranked enrichment on the Welch t, gene ranks tie-free under SM15 (exploratory)
      stat <- stats::setNames(w$statistic, w$gene_id)
      ranked <- prerank_welch_ranks(stat[!is.na(names(stat)) & nzchar(names(stat))], tie_break)
      for (set_id in names(eligible_sets)) {
        members <- eligible_sets[[set_id]]
        k <- length(members)
        if (k < as.integer(p$min_size) || k > as.integer(p$max_size) || k >= length(ranked)) next
        idx <- which(names(ranked) %in% members)
        res <- prerank_set_result(ranked, idx, set_id, as.integer(request$rng$seed))
        rows[[length(rows) + 1L]] <- data.frame(library_id = library$id, method = "prerank_welch", contrast_id = contrast_id, set_id = set_id, N = length(ranked),
                                               K = k, n = NA_integer_, k = length(idx), statistic = NA_real_, es = res$es, p_value = res$p_value,
                                               exploratory = TRUE, estimator = res$method, stringsAsFactors = FALSE)
      }
    }
    combined <- do.call(rbind, rows)
    if (is.null(combined)) stop("E_PATHWAY_NO_ELIGIBLE_SETS: no eligible set for library ", library$id, call. = FALSE)
    combined$family_id <- NA_character_; combined$q_value <- NA_real_
    combined <- enrichment_family_adjust(combined, gene_set_family_key)
    combined$claim_label <- ifelse(combined$exploratory, "exploratory", "measured_background_ora")
    combined$result_type <- "GeneSetResult"; combined$schema_version <- "1.2.0"
    ora_part <- combined[startsWith(combined$method, "ora_"), , drop = FALSE]
    prerank_part <- combined[startsWith(combined$method, "prerank"), , drop = FALSE]
    if (nrow(ora_part)) write_tsv(ora_part, file.path("pathways", paste0("ora_", library$id, ".tsv")), paste0("pathways_ora_", library$id), "GeneSetResult")
    if (nrow(prerank_part)) write_tsv(prerank_part, file.path("pathways", paste0("prerank_", library$id, ".tsv")), paste0("pathways_prerank_", library$id), "GeneSetResult")
  }
  manifest <- lapply(libraries, function(library) list(id = library$id, sha256 = library$sha256, source = library$source, version = library$version,
                                                       species = library$species, identifier_type = library$identifier_type, licence_id = library$licence_id,
                                                       redistributable = isTRUE(library$redistributable), synthetic = isTRUE(library$synthetic)))
  write_json(list(result_type = "LibraryManifest", libraries = manifest, library_content_written = FALSE,
                  network = "none", versions = list(R = as.character(getRversion()))), file.path("pathways", "libraries.json"), "pathways_libraries", "LibraryManifest")
  list(outputs = outputs, warnings = warnings, message = sprintf("gene-set libraries: %d librar(ies) read from hashed local files", length(libraries)))
})
