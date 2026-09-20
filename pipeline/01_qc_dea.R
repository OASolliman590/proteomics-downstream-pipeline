#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(limma)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(impute)
  library(matrixStats)
})
set.seed(48)

args <- commandArgs(FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
ROOT <- Sys.getenv("DAPA_ANALYSIS_ROOT", unset="")
if (!nzchar(ROOT)) ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[1])), ".."))
ROOT <- normalizePath(ROOT, mustWork=TRUE)
PROC <- file.path(ROOT, "data", "processed")
TAB <- file.path(ROOT, "results", "tables")
FIG <- file.path(ROOT, "results", "figures")
REPORT <- file.path(ROOT, "report")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(REPORT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(TAB, "dea"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(FIG, "qc"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(FIG, "dea"), recursive = TRUE, showWarnings = FALSE)

meta <- read_csv(file.path(PROC, "sample_meta.csv"), show_col_types = FALSE)
meta$group <- factor(meta$group, levels = c("CTRL", "PreDM", "DM", "PreDM_Dapa", "DM_Dapa"))
raw_tbl <- read_tsv(file.path(PROC, "abundance_input_log2.tsv"), show_col_types = FALSE)
annot <- read_tsv(file.path(PROC, "protein_annot.tsv"), show_col_types = FALSE)
stopifnot(all(meta$sample_id %in% names(raw_tbl)), !anyDuplicated(meta$sample_id))
mat <- as.matrix(raw_tbl[, meta$sample_id])
storage.mode(mat) <- "numeric"
rownames(mat) <- raw_tbl$Accession

# This input is already log2 processed. Preserve it as the primary scale.
coverage <- sapply(levels(meta$group), function(g) {
  rowSums(!is.na(mat[, meta$sample_id[meta$group == g], drop = FALSE]))
})
colnames(coverage) <- levels(meta$group)
keep <- apply(coverage, 1, max) >= 3
mat_f <- mat[keep, , drop = FALSE]
annot_f <- annot[match(rownames(mat_f), annot$Accession), ]

write_tsv(
  as_tibble(mat_f, rownames = "Accession"),
  file.path(PROC, "abundance_analysis_log2.tsv")
)
write_tsv(
  bind_cols(tibble(Accession = rownames(mat)), as_tibble(coverage), tibble(keep_analysis = keep)),
  file.path(TAB, "protein_coverage.tsv")
)

# QC summaries. PCA-only median imputation does not alter the primary DEA matrix.
sample_qc <- tibble(
  sample_id = colnames(mat_f),
  group = meta$group[match(colnames(mat_f), meta$sample_id)],
  n_missing = colSums(is.na(mat_f)),
  pct_missing = colMeans(is.na(mat_f)) * 100,
  median_log2 = colMedians(mat_f, na.rm = TRUE),
  iqr_log2 = colIQRs(mat_f, na.rm = TRUE)
)

mat_pca <- mat_f
row_medians <- rowMedians(mat_pca, na.rm = TRUE)
for (j in seq_len(ncol(mat_pca))) {
  missing <- is.na(mat_pca[, j])
  mat_pca[missing, j] <- row_medians[missing]
}
variable <- rowSds(mat_pca) > 0
pca <- prcomp(t(mat_pca[variable, , drop = FALSE]), center = TRUE, scale. = FALSE)
pca_tbl <- as_tibble(pca$x[, 1:min(5, ncol(pca$x)), drop = FALSE], rownames = "sample_id") |>
  left_join(meta |> select(sample_id, group), by = "sample_id")
variance_explained <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
cor_mat <- cor(mat_f, use = "pairwise.complete.obs", method = "pearson")

within_group_cor <- sapply(seq_len(nrow(meta)), function(i) {
  peers <- which(meta$group == meta$group[i] & seq_len(nrow(meta)) != i)
  median(cor_mat[i, peers], na.rm = TRUE)
})
pca_xy <- pca_tbl |> select(sample_id, PC1, PC2, group)
centroids <- pca_xy |> group_by(group) |> summarise(across(c(PC1, PC2), mean), .groups = "drop")
pca_dist <- pca_xy |> left_join(centroids, by = "group", suffix = c("", "_centroid")) |>
  mutate(pca_distance = sqrt((PC1-PC1_centroid)^2 + (PC2-PC2_centroid)^2)) |>
  select(sample_id, pca_distance)
sample_qc <- sample_qc |>
  mutate(median_within_group_correlation = within_group_cor) |>
  left_join(pca_dist, by = "sample_id") |>
  mutate(
    missingness_robust_z = (pct_missing - median(pct_missing)) / pmax(mad(pct_missing), 1e-8),
    pca_distance_robust_z = (pca_distance - median(pca_distance)) / pmax(mad(pca_distance), 1e-8),
    technical_watch = abs(missingness_robust_z) >= 3 | abs(pca_distance_robust_z) >= 3
  )
write_tsv(sample_qc, file.path(TAB, "sample_qc.tsv"))
write_tsv(pca_tbl, file.path(TAB, "pca_scores.tsv"))
write_tsv(as_tibble(cor_mat, rownames = "sample_id"), file.path(TAB, "sample_correlations.tsv"))

group_colors <- c(CTRL="#4D4D4D", PreDM="#E69F00", DM="#D55E00", PreDM_Dapa="#56B4E9", DM_Dapa="#0072B2")
p_pca <- ggplot(pca_tbl, aes(PC1, PC2, colour = group, label = sample_id)) +
  geom_point(size = 3) + geom_text_repel(size = 2.5, max.overlaps = Inf) +
  scale_colour_manual(values = group_colors) + theme_bw(base_size = 11) +
  labs(title = "PCA of processed log2 proteome", x = sprintf("PC1 (%.1f%%)", variance_explained[1]), y = sprintf("PC2 (%.1f%%)", variance_explained[2]))
ggsave(file.path(FIG, "qc", "pca.pdf"), p_pca, width = 8, height = 6)
ggsave(file.path(FIG, "qc", "pca.png"), p_pca, width = 8, height = 6, dpi = 180)

p_miss <- ggplot(sample_qc, aes(reorder(sample_id, pct_missing), pct_missing, fill = group)) +
  geom_col() + coord_flip() + scale_fill_manual(values = group_colors) + theme_bw(base_size = 10) +
  labs(title = "Sample-level missingness", x = NULL, y = "Missing proteins (%)")
ggsave(file.path(FIG, "qc", "sample_missingness.pdf"), p_miss, width = 8, height = 6)

pdf(file.path(FIG, "qc", "sample_correlation_heatmap.pdf"), width = 9, height = 8)
pheatmap(cor_mat, annotation_col = data.frame(group=meta$group, row.names=meta$sample_id), border_color = NA)
dev.off()

# Joint five-group limma model. No batch term is invented because none was supplied.
design <- model.matrix(~ 0 + group, data = meta)
colnames(design) <- sub("^group", "", colnames(design))
contrast_mat <- makeContrasts(
  D1_PreDM_vs_CTRL = PreDM - CTRL,
  D2_DM_vs_CTRL = DM - CTRL,
  T1_PreDM_Dapa_vs_PreDM = PreDM_Dapa - PreDM,
  T2_DM_Dapa_vs_DM = DM_Dapa - DM,
  R1_PreDM_Dapa_vs_CTRL = PreDM_Dapa - CTRL,
  R2_DM_Dapa_vs_CTRL = DM_Dapa - CTRL,
  S1_DM_vs_PreDM = DM - PreDM,
  S2_DM_Dapa_vs_PreDM_Dapa = DM_Dapa - PreDM_Dapa,
  I1_stage_by_dapa = (DM_Dapa - DM) - (PreDM_Dapa - PreDM),
  levels = design
)

fit_matrix <- function(x) {
  fit <- lmFit(x, design)
  eBayes(contrasts.fit(fit, contrast_mat), trend = TRUE, robust = TRUE)
}
fit2 <- fit_matrix(mat_f)
saveRDS(fit2, file.path(TAB, "joint_limma_fit.rds"))
write_tsv(as_tibble(design, rownames = "sample_id"), file.path(TAB, "design_matrix.tsv"))
write_tsv(as_tibble(contrast_mat, rownames = "coefficient"), file.path(TAB, "contrast_matrix.tsv"))

contrast_names <- colnames(contrast_mat)
dea_tables <- list()
summary_rows <- list()
for (i in seq_along(contrast_names)) {
  cn <- contrast_names[i]
  tt <- topTable(fit2, coef = i, number = Inf, sort.by = "none") |>
    as_tibble(rownames = "Accession") |>
    left_join(annot_f |> select(Accession, gene_symbol, description), by = "Accession")
  for (g in colnames(coverage)) tt[[paste0("n_", g)]] <- coverage[tt$Accession, g]
  tt <- tt |>
    mutate(
      contrast = cn,
      formal_DEP = adj.P.Val < 0.05 & abs(logFC) >= 0.25,
      supportive_FDR = adj.P.Val < 0.10 & abs(logFC) >= 0.25,
      effect_candidate = P.Value < 0.05 & abs(logFC) >= 0.58,
      evidence_tier = case_when(
        formal_DEP ~ "formal_FDR_DEP",
        supportive_FDR ~ "supportive_FDR_0.10",
        effect_candidate ~ "rawP_effect_candidate",
        TRUE ~ "not_selected"
      )
    ) |>
    relocate(Accession, gene_symbol, description, contrast)
  dea_tables[[cn]] <- tt
  write_tsv(tt, file.path(TAB, "dea", paste0("dea_", cn, ".tsv")))
  summary_rows[[cn]] <- tibble(
    contrast = cn,
    n_tested = nrow(tt),
    n_formal_up = sum(tt$formal_DEP & tt$logFC > 0, na.rm = TRUE),
    n_formal_down = sum(tt$formal_DEP & tt$logFC < 0, na.rm = TRUE),
    n_supportive_fdr = sum(tt$supportive_FDR, na.rm = TRUE),
    n_rawP_effect_candidates = sum(tt$effect_candidate, na.rm = TRUE),
    min_p_value = min(tt$P.Value, na.rm = TRUE),
    min_adj_p_value = min(tt$adj.P.Val, na.rm = TRUE)
  )
  p_volcano <- ggplot(tt, aes(logFC, -log10(pmax(P.Value, 1e-300)), colour = evidence_tier)) +
    geom_point(alpha = 0.65, size = 1) +
    geom_vline(xintercept = c(-0.58, 0.58), linetype = 2, colour = "grey45") +
    geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey45") +
    scale_colour_manual(values = c(formal_FDR_DEP="#B2182B", supportive_FDR_0.10="#EF8A62", rawP_effect_candidate="#2166AC", not_selected="#BDBDBD")) +
    theme_bw(base_size = 10) + labs(title = cn, x = "log2 fold change", y = "-log10(raw P)")
  ggsave(file.path(FIG, "dea", paste0("volcano_", cn, ".pdf")), p_volcano, width = 7, height = 5.5)
}
dea_summary <- bind_rows(summary_rows)
write_tsv(dea_summary, file.path(TAB, "dea_summary.tsv"))

master <- reduce(lapply(names(dea_tables), function(cn) {
  dea_tables[[cn]] |> select(Accession, logFC, P.Value, adj.P.Val) |>
    rename(!!paste0("logFC_", cn) := logFC, !!paste0("P_", cn) := P.Value, !!paste0("adjP_", cn) := adj.P.Val)
}), full_join, by = "Accession") |>
  left_join(annot_f |> select(Accession, gene_symbol, description), by = "Accession") |>
  relocate(Accession, gene_symbol, description)
write_tsv(master, file.path(TAB, "dea_master.tsv"))

# Top-candidate heatmap, explicitly exploratory when formal DEPs are absent.
primary <- c("D1_PreDM_vs_CTRL", "D2_DM_vs_CTRL", "T1_PreDM_Dapa_vs_PreDM", "T2_DM_Dapa_vs_DM")
top_acc <- unique(unlist(lapply(primary, function(cn) dea_tables[[cn]] |> arrange(P.Value) |> slice_head(n=15) |> pull(Accession))))
top_acc <- intersect(top_acc, rownames(mat_pca))
heat <- t(scale(t(mat_pca[top_acc, , drop = FALSE])))
rownames(heat) <- make.unique(ifelse(is.na(annot_f$gene_symbol[match(top_acc, annot_f$Accession)]) | annot_f$gene_symbol[match(top_acc, annot_f$Accession)] == "", top_acc, annot_f$gene_symbol[match(top_acc, annot_f$Accession)]))
pdf(file.path(FIG, "dea", "top_candidate_heatmap.pdf"), width = 10, height = 11)
pheatmap(heat, annotation_col = data.frame(group=meta$group, row.names=meta$sample_id), show_colnames=TRUE, fontsize_row=6, border_color=NA)
dev.off()

# Sensitivity analyses: sample median centering, KNN, and MinDet imputation.
grand_median <- median(mat_f, na.rm = TRUE)
median_centered <- sweep(mat_f, 2, colMedians(mat_f, na.rm = TRUE) - grand_median, "-")
knn_mat <- tryCatch(impute.knn(mat_f, k = 5, rng.seed = 48)$data, error = function(e) NULL)
mindet_mat <- mat_f
for (j in seq_len(ncol(mindet_mat))) {
  observed <- mindet_mat[, j][is.finite(mindet_mat[, j])]
  mu <- median(observed) - 1.8 * sd(observed)
  sigma <- 0.3 * sd(observed)
  miss <- is.na(mindet_mat[, j])
  if (any(miss)) mindet_mat[miss, j] <- rnorm(sum(miss), mu, sigma)
}
sensitivity_mats <- list(source_log2_primary = mat_f, median_centered = median_centered, MinDet_1.8_0.3 = mindet_mat)
if (!is.null(knn_mat)) sensitivity_mats$KNN_k5 <- knn_mat
sensitivity_rows <- list()
for (method in names(sensitivity_mats)) {
  sf <- fit_matrix(sensitivity_mats[[method]])
  for (i in seq_along(contrast_names)) {
    cn <- contrast_names[i]
    base_lfc <- fit2$coefficients[, i]
    sens_lfc <- sf$coefficients[, i]
    sensitivity_rows[[paste(method, cn)]] <- tibble(
      method = method,
      contrast = cn,
      pearson_logFC_vs_primary = cor(base_lfc, sens_lfc, use = "pairwise.complete.obs"),
      spearman_logFC_vs_primary = cor(base_lfc, sens_lfc, use = "pairwise.complete.obs", method = "spearman"),
      n_fdr_0_05_abs_lfc_0_25 = sum(sf$p.value[,i] |> p.adjust(method="BH") < 0.05 & abs(sens_lfc) >= 0.25, na.rm=TRUE),
      n_rawP_0_05_abs_lfc_0_58 = sum(sf$p.value[,i] < 0.05 & abs(sens_lfc) >= 0.58, na.rm=TRUE)
    )
  }
}
write_tsv(bind_rows(sensitivity_rows), file.path(TAB, "method_sensitivity_summary.tsv"))

# Leave-one-sample-out stability for the four primary contrasts.
loo_summary <- list()
candidate_stability <- list()
for (sample_out in meta$sample_id) {
  keep_samples <- meta$sample_id != sample_out
  design_loo <- model.matrix(~ 0 + group, data = droplevels(meta[keep_samples, ]))
  colnames(design_loo) <- sub("^group", "", colnames(design_loo))
  cmat_loo <- makeContrasts(
    D1_PreDM_vs_CTRL = PreDM - CTRL,
    D2_DM_vs_CTRL = DM - CTRL,
    T1_PreDM_Dapa_vs_PreDM = PreDM_Dapa - PreDM,
    T2_DM_Dapa_vs_DM = DM_Dapa - DM,
    levels = design_loo
  )
  lf <- eBayes(contrasts.fit(lmFit(mat_f[, keep_samples, drop=FALSE], design_loo), cmat_loo), trend=TRUE, robust=TRUE)
  for (cn in primary) {
    idx_full <- match(cn, contrast_names)
    idx_loo <- match(cn, colnames(cmat_loo))
    full_lfc <- fit2$coefficients[, idx_full]
    loo_lfc <- lf$coefficients[, idx_loo]
    loo_p <- lf$p.value[, idx_loo]
    loo_summary[[paste(sample_out, cn)]] <- tibble(
      omitted_sample = sample_out,
      omitted_group = as.character(meta$group[meta$sample_id == sample_out]),
      contrast = cn,
      pearson_logFC = cor(full_lfc, loo_lfc, use="pairwise.complete.obs"),
      sign_agreement = mean(sign(full_lfc) == sign(loo_lfc), na.rm=TRUE),
      n_rawP_effect_candidates = sum(loo_p < 0.05 & abs(loo_lfc) >= 0.58, na.rm=TRUE)
    )
    full_candidates <- dea_tables[[cn]]$Accession[dea_tables[[cn]]$effect_candidate]
    if (length(full_candidates) > 0) {
      candidate_stability[[paste(sample_out, cn)]] <- tibble(
        omitted_sample = sample_out,
        contrast = cn,
        Accession = full_candidates,
        same_direction = sign(loo_lfc[full_candidates]) == sign(full_lfc[full_candidates]),
        retains_candidate = loo_p[full_candidates] < 0.05 & abs(loo_lfc[full_candidates]) >= 0.58
      )
    }
  }
}
write_tsv(bind_rows(loo_summary), file.path(TAB, "leave_one_sample_summary.tsv"))
if (length(candidate_stability) > 0) {
  stability <- bind_rows(candidate_stability) |>
    group_by(contrast, Accession) |>
    summarise(fraction_same_direction = mean(same_direction, na.rm=TRUE), fraction_retains_candidate = mean(retains_candidate, na.rm=TRUE), .groups="drop") |>
    left_join(annot_f |> select(Accession, gene_symbol), by="Accession")
  write_tsv(stability, file.path(TAB, "candidate_leave_one_sample_stability.tsv"))
}

qc_notes <- c(
  "# QC and Differential Analysis Notes", "",
  sprintf("- Input proteins: %d; retained after >=3/4 coverage in at least one group: %d.", nrow(mat), nrow(mat_f)),
  "- Input was already processed on a log2 scale; no additional log transformation was applied.",
  "- Primary limma model: ~0 + group. Batch/wave was not supplied and was not invented.",
  "- Primary model retains missing values and lets limma use available observations.",
  "- Median-centering, KNN, and MinDet are sensitivity analyses only.",
  sprintf("- Samples on the technical watchlist: %d. No sample was automatically excluded.", sum(sample_qc$technical_watch)),
  "- Formal DEP: BH-FDR < 0.05 and |log2FC| >= 0.25.",
  "- Raw-P/effect-size candidates are exploratory and are not called DEPs."
)
writeLines(qc_notes, file.path(REPORT, "qc_dea_notes.md"))
message("[01] QC and joint limma analysis complete")
