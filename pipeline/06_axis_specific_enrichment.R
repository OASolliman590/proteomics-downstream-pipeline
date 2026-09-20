#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(fgsea)
  library(msigdbr)
})

set.seed(48)
ROOT <- Sys.getenv("DAPA_ANALYSIS_ROOT", unset = "")
if (!nzchar(ROOT)) stop("DAPA_ANALYSIS_ROOT must be set")
ROOT <- normalizePath(ROOT, mustWork = TRUE)
OUT <- file.path(ROOT, "results", "tables", "axis_specific")

load_sets <- function(collection, subcollection = NULL) {
  x <- suppressMessages(msigdbr(species = "Rattus norvegicus", collection = collection,
                                subcollection = subcollection))
  lapply(split(x$gene_symbol, x$gs_name), unique)
}

collections <- list(
  KEGG_LEGACY = load_sets("C2", "CP:KEGG_LEGACY"),
  REACTOME = load_sets("C2", "CP:REACTOME"),
  HALLMARK = load_sets("H")
)

dea <- read_tsv(file.path(OUT, "axis_dea_all.tsv"), show_col_types = FALSE)
primary <- c("PreDM_disease_vs_CTRL", "DM_disease_vs_CTRL",
             "PreDM_Dapa_vs_PreDM", "DM_Dapa_vs_DM")

collapse_ranks <- function(x) {
  ranked <- x |>
    filter(!is.na(gene_symbol), gene_symbol != "", is.finite(P.Value), is.finite(logFC)) |>
    mutate(rank_value = sign(logFC) * -log10(pmax(P.Value, 1e-300))) |>
    group_by(gene_symbol) |>
    slice_max(abs(rank_value), n = 1, with_ties = FALSE) |>
    ungroup()
  sort(setNames(ranked$rank_value, ranked$gene_symbol), decreasing = TRUE)
}

rows <- list()
for (cn in primary) {
  message("[06] ", cn)
  ranks <- collapse_ranks(dea |> filter(contrast == cn))
  for (coll in names(collections)) {
    sets <- lapply(collections[[coll]], intersect, y = names(ranks))
    sets <- sets[lengths(sets) >= 5 & lengths(sets) <= 500]
    if (!length(sets)) next
    fg <- suppressWarnings(fgseaMultilevel(pathways = sets, stats = ranks,
                                           minSize = 5, maxSize = 500, eps = 0)) |>
      as_tibble() |>
      mutate(contrast = cn, collection = coll,
             leadingEdge = vapply(leadingEdge, paste, collapse = ";", character(1))) |>
      select(contrast, collection, pathway, size, ES, NES, pval, padj, leadingEdge)
    rows[[paste(cn, coll)]] <- fg
  }
}

gsea <- bind_rows(rows)
write_tsv(gsea, file.path(OUT, "axis_gsea_all.tsv"))
summary <- gsea |>
  group_by(contrast, collection) |>
  summarise(n_tested = n(), n_FDR_0_05 = sum(padj < 0.05, na.rm = TRUE),
            n_FDR_0_10 = sum(padj < 0.10, na.rm = TRUE),
            min_FDR = min(padj, na.rm = TRUE), .groups = "drop")
write_tsv(summary, file.path(OUT, "axis_gsea_summary.tsv"))
message("[06] axis-specific GSEA complete")
