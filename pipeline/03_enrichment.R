#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(fgsea)
  library(msigdbr)
  library(limma)
})
set.seed(48)

args <- commandArgs(FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
ROOT <- Sys.getenv("DAPA_ANALYSIS_ROOT", unset="")
if (!nzchar(ROOT)) ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[1])), ".."))
ROOT <- normalizePath(ROOT, mustWork=TRUE)
TAB <- file.path(ROOT, "results", "tables")
FIG <- file.path(ROOT, "results", "figures", "enrichment")
REPORT <- file.path(ROOT, "report")
dir.create(FIG, recursive=TRUE, showWarnings=FALSE)

primary <- c("D1_PreDM_vs_CTRL", "D2_DM_vs_CTRL", "T1_PreDM_Dapa_vs_PreDM", "T2_DM_Dapa_vs_DM")

load_sets <- function(collection, subcollection=NULL) {
  tryCatch({
    message(sprintf("[03] loading %s/%s", collection, ifelse(is.null(subcollection), "all", subcollection)))
    x <- suppressMessages(msigdbr(species="Rattus norvegicus", collection=collection, subcollection=subcollection))
    out <- lapply(split(x$gene_symbol, x$gs_name), unique)
    message(sprintf("[03] loaded %d pathways", length(out)))
    out
  }, error=function(e) {
    message(sprintf("[03] collection unavailable %s/%s: %s", collection, subcollection, conditionMessage(e)))
    list()
  })
}
collections <- list(
  KEGG_LEGACY=load_sets("C2", "CP:KEGG_LEGACY"),
  REACTOME=load_sets("C2", "CP:REACTOME"),
  HALLMARK=load_sets("H")
)
collections <- collections[lengths(collections) > 0]
if (length(collections) == 0) stop("No local msigdbr collections were available")

collapse_ranks <- function(dea) {
  ranked <- dea |>
    filter(!is.na(gene_symbol), gene_symbol != "", is.finite(P.Value), is.finite(logFC)) |>
    mutate(rank_value=sign(logFC) * -log10(pmax(P.Value, 1e-300))) |>
    group_by(gene_symbol) |>
    slice_max(abs(rank_value), n=1, with_ties=FALSE) |>
    ungroup()
  sort(setNames(ranked$rank_value, ranked$gene_symbol), decreasing=TRUE)
}

gsea_rows <- list()
ora_rows <- list()
for (cn in primary) {
  message(sprintf("[03] contrast %s", cn))
  dea <- read_tsv(file.path(TAB, "dea", paste0("dea_", cn, ".tsv")), show_col_types=FALSE)
  ranks <- collapse_ranks(dea)
  background <- names(ranks)
  candidate_genes <- unique(dea$gene_symbol[dea$effect_candidate & !is.na(dea$gene_symbol) & dea$gene_symbol != ""])
  for (coll in names(collections)) {
    sets <- collections[[coll]]
    sets <- lapply(sets, intersect, y=background)
    sets <- sets[lengths(sets) >= 5 & lengths(sets) <= 500]
    if (length(sets) == 0) next
    fg <- suppressWarnings(fgseaMultilevel(pathways=sets, stats=ranks, minSize=5, maxSize=500, eps=0)) |>
      as_tibble() |>
      mutate(contrast=cn, collection=coll, leadingEdge=vapply(leadingEdge, paste, collapse=";", character(1))) |>
      select(contrast, collection, pathway, size, ES, NES, pval, padj, leadingEdge)
    gsea_rows[[paste(cn,coll)]] <- fg

    if (length(candidate_genes) >= 3) {
      N <- length(background); n <- length(intersect(candidate_genes, background))
      ora <- bind_rows(lapply(names(sets), function(pathway) {
        genes <- sets[[pathway]]
        K <- length(genes); overlap <- intersect(candidate_genes, genes); k <- length(overlap)
        tibble(pathway=pathway, set_size=K, candidate_size=n, overlap_size=k,
               pvalue=phyper(k-1, K, N-K, n, lower.tail=FALSE),
               overlap_genes=paste(overlap, collapse=";"))
      })) |>
        mutate(padj=p.adjust(pvalue, method="BH"), contrast=cn, collection=coll) |>
        filter(overlap_size > 0) |>
        select(contrast, collection, everything())
      ora_rows[[paste(cn,coll)]] <- ora
    }
  }
}

gsea_all <- bind_rows(gsea_rows)
ora_all <- bind_rows(ora_rows)
write_tsv(gsea_all, file.path(TAB, "gsea_all.tsv"))
write_tsv(ora_all, file.path(TAB, "candidate_ora_all.tsv"))

# Focused leave-one-sample stability for pathways significant in either direct
# dapagliflozin contrast. Re-test only the full-run significant pathways; the
# resulting nominal P values are sensitivity diagnostics, not a second FDR screen.
sig_treatment <- gsea_all |>
  filter(contrast %in% c("T1_PreDM_Dapa_vs_PreDM", "T2_DM_Dapa_vs_DM"), padj < 0.05)
if (nrow(sig_treatment) > 0) {
  meta <- read_csv(file.path(ROOT, "data", "processed", "sample_meta.csv"), show_col_types=FALSE)
  meta$group <- factor(meta$group, levels=c("CTRL","PreDM","DM","PreDM_Dapa","DM_Dapa"))
  mat_tbl <- read_tsv(file.path(ROOT, "data", "processed", "abundance_analysis_log2.tsv"), show_col_types=FALSE)
  mat <- as.matrix(mat_tbl[, meta$sample_id]); storage.mode(mat) <- "numeric"; rownames(mat) <- mat_tbl$Accession
  annot <- read_tsv(file.path(ROOT, "data", "processed", "protein_annot.tsv"), show_col_types=FALSE)
  symbol <- setNames(annot$gene_symbol, annot$Accession)
  cfg <- list(
    T1_PreDM_Dapa_vs_PreDM=list(expr="PreDM_Dapa-PreDM", groups=c("PreDM","PreDM_Dapa")),
    T2_DM_Dapa_vs_DM=list(expr="DM_Dapa-DM", groups=c("DM","DM_Dapa"))
  )
  loo_pathway_rows <- list()
  for (cn in names(cfg)) {
    sig <- sig_treatment |> filter(contrast == cn)
    if (nrow(sig) == 0) next
    selected_sets <- list()
    for (i in seq_len(nrow(sig))) {
      key <- paste(sig$collection[i], sig$pathway[i], sep="::")
      selected_sets[[key]] <- collections[[sig$collection[i]]][[sig$pathway[i]]]
    }
    omissions <- meta$sample_id[meta$group %in% cfg[[cn]]$groups]
    for (sample_out in omissions) {
      keep <- meta$sample_id != sample_out
      design <- model.matrix(~0+group, data=meta[keep,])
      colnames(design) <- sub("^group", "", colnames(design))
      cm <- makeContrasts(contrasts=cfg[[cn]]$expr, levels=design)
      fit <- eBayes(contrasts.fit(lmFit(mat[,keep,drop=FALSE], design), cm), trend=TRUE, robust=TRUE)
      tt <- topTable(fit, coef=1, number=Inf, sort.by="none") |> as_tibble(rownames="Accession") |>
        mutate(gene_symbol=symbol[Accession])
      ranks <- collapse_ranks(tt)
      sets_use <- lapply(selected_sets, intersect, y=names(ranks))
      fg <- suppressWarnings(fgseaMultilevel(pathways=sets_use, stats=ranks, minSize=5, maxSize=500, eps=0)) |> as_tibble()
      loo_pathway_rows[[paste(cn,sample_out)]] <- fg |>
        separate(pathway, into=c("collection","pathway"), sep="::", extra="merge") |>
        mutate(contrast=cn, omitted_sample=sample_out) |>
        select(contrast, omitted_sample, collection, pathway, NES, pval)
    }
  }
  loo_pathway <- bind_rows(loo_pathway_rows) |>
    left_join(sig_treatment |> select(contrast, collection, pathway, full_NES=NES, full_FDR=padj), by=c("contrast","collection","pathway")) |>
    mutate(same_direction=sign(NES)==sign(full_NES))
  write_tsv(loo_pathway, file.path(TAB, "gsea_leave_one_sample_results.tsv"))
  loo_pathway_summary <- loo_pathway |>
    group_by(contrast, collection, pathway, full_NES, full_FDR) |>
    summarise(n_omissions=n(), fraction_same_direction=mean(same_direction,na.rm=TRUE), fraction_nominal_p_0_05=mean(pval<0.05,na.rm=TRUE), min_NES=min(NES,na.rm=TRUE), max_NES=max(NES,na.rm=TRUE), .groups="drop")
  write_tsv(loo_pathway_summary, file.path(TAB, "gsea_leave_one_sample_stability.tsv"))
}

gsea_summary <- gsea_all |>
  group_by(contrast, collection) |>
  summarise(n_tested=n(), n_fdr_0_05=sum(padj<0.05, na.rm=TRUE), n_fdr_0_10=sum(padj<0.10, na.rm=TRUE), min_fdr=min(padj, na.rm=TRUE), .groups="drop")
write_tsv(gsea_summary, file.path(TAB, "gsea_summary.tsv"))

top_plot <- gsea_all |>
  group_by(contrast) |>
  arrange(padj, pval, .by_group=TRUE) |>
  slice_head(n=12) |>
  ungroup() |>
  mutate(pathway_short=gsub("^(KEGG_|REACTOME_|GOBP_)", "", pathway), pathway_short=gsub("_", " ", pathway_short),
         contrast_label=recode(contrast,
           D1_PreDM_vs_CTRL="PreDM",
           D2_DM_vs_CTRL="DM",
           T1_PreDM_Dapa_vs_PreDM="PreDM+D",
           T2_DM_Dapa_vs_DM="DM+D"))
if (nrow(top_plot) > 0) {
  p <- ggplot(top_plot, aes(NES, reorder(pathway_short, NES), colour=padj, size=-log10(pmax(pval,1e-300)))) +
    geom_point() + facet_wrap(~contrast_label, scales="free_y", ncol=2) +
    scale_colour_viridis_c(direction=-1) + theme_bw(base_size=9) +
    labs(title="Top rank-based pathway results", x="Normalized enrichment score", y=NULL)
  ggsave(file.path(FIG, "top_gsea_pathways.pdf"), p, width=12, height=10)
  ggsave(file.path(FIG, "top_gsea_pathways.png"), p, width=12, height=10, dpi=180)
}

notes <- c(
  "# Pathway Analysis Notes", "",
  "- Rank-based fgsea was run independently for PreDM disease, DM disease, PreDM dapagliflozin, and DM dapagliflozin contrasts.",
  "- Ranks use signed -log10(raw P), with duplicate gene symbols collapsed to the largest absolute rank.",
  "- Candidate ORA uses raw-P/effect-size protein candidates and is exploratory, not confirmatory.",
  "- FDR-significant direct-treatment pathways were subjected to focused leave-one-sample stability analysis.",
  "- Pathway evidence does not convert raw-P protein candidates into formal DEPs.",
  "- Tissue is not supplied, so pathway interpretation must remain tissue-agnostic.", "",
  paste(capture.output(print(gsea_summary)), collapse="\n")
)
writeLines(notes, file.path(REPORT, "enrichment_notes.md"))
message("[03] rank-based and candidate pathway analyses complete")
