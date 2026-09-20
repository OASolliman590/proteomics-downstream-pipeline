#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(limma)
  library(ggplot2)
})

set.seed(48)
ROOT <- Sys.getenv("DAPA_ANALYSIS_ROOT", unset = "")
if (!nzchar(ROOT)) stop("DAPA_ANALYSIS_ROOT must be set")
ROOT <- normalizePath(ROOT, mustWork = TRUE)
PROC <- file.path(ROOT, "data", "processed")
OUT <- file.path(ROOT, "results", "tables", "axis_specific")
FIG <- file.path(ROOT, "results", "figures", "axis_specific")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

meta <- read_csv(file.path(PROC, "sample_meta.csv"), show_col_types = FALSE)
abund <- read_tsv(file.path(PROC, "abundance_analysis_log2.tsv"), show_col_types = FALSE)
annot <- read_tsv(file.path(PROC, "protein_annot.tsv"), show_col_types = FALSE)
mat_all <- as.matrix(abund[, meta$sample_id])
storage.mode(mat_all) <- "numeric"
rownames(mat_all) <- abund$Accession

axis_cfg <- list(
  PreDM = list(groups = c("CTRL", "PreDM", "PreDM_Dapa"), disease = "PreDM-CTRL",
               treatment = "PreDM_Dapa-PreDM", residual = "PreDM_Dapa-CTRL"),
  DM = list(groups = c("CTRL", "DM", "DM_Dapa"), disease = "DM-CTRL",
            treatment = "DM_Dapa-DM", residual = "DM_Dapa-CTRL")
)

safe_top <- function(fit, coef_name, label, coverage) {
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none") |>
    as_tibble() |>
    mutate(Accession = rownames(fit$coefficients), .before = 1) |>
    left_join(annot, by = "Accession") |>
    mutate(
      contrast = label,
      formal_DEP = adj.P.Val < 0.05 & abs(logFC) >= 0.25,
      supportive_FDR = adj.P.Val < 0.10 & abs(logFC) >= 0.25,
      nominal_effect_candidate = P.Value < 0.05 & abs(logFC) >= 0.58
    ) |>
    relocate(Accession, gene_symbol, description, contrast)
  for (nm in names(coverage)) tt[[paste0("n_", nm)]] <- coverage[[nm]]
  tt
}

all_de <- list()
all_rev <- list()
all_omnibus <- list()
summary_rows <- list()
reversal_summary_rows <- list()
axis_score_rows <- list()
axis_test_rows <- list()

exact_permutation_p <- function(x, y) {
  observed <- mean(y) - mean(x)
  pooled <- c(x, y)
  assignments <- combn(seq_along(pooled), length(x))
  null <- apply(assignments, 2, function(idx) mean(pooled[-idx]) - mean(pooled[idx]))
  (sum(abs(null) >= abs(observed)) + 1) / (length(null) + 1)
}

for (axis_name in names(axis_cfg)) {
  cfg <- axis_cfg[[axis_name]]
  m <- meta |> filter(group %in% cfg$groups)
  m$group <- factor(m$group, levels = cfg$groups)
  x <- mat_all[, m$sample_id, drop = FALSE]
  keep <- apply(x, 1, function(z) {
    all(vapply(cfg$groups, function(g) sum(is.finite(z[m$group == g])) >= 3, logical(1)))
  })
  x <- x[keep, , drop = FALSE]
  coverage <- as.data.frame(sapply(cfg$groups, function(g) {
    rowSums(is.finite(x[, m$group == g, drop = FALSE]))
  }))
  names(coverage) <- cfg$groups

  design <- model.matrix(~0 + group, data = m)
  colnames(design) <- levels(m$group)
  cm <- makeContrasts(contrasts = c(cfg$disease, cfg$treatment, cfg$residual), levels = design)
  colnames(cm) <- c("disease", "treatment", "residual")
  fit <- lmFit(x, design)
  fit2 <- eBayes(contrasts.fit(fit, cm), trend = TRUE, robust = TRUE)
  disease <- safe_top(fit2, "disease", paste0(axis_name, "_disease_vs_CTRL"), coverage)
  treatment <- safe_top(fit2, "treatment", paste0(axis_name, "_Dapa_vs_", axis_name), coverage)
  residual <- safe_top(fit2, "residual", paste0(axis_name, "_Dapa_vs_CTRL"), coverage)
  axis_de <- bind_rows(disease, treatment, residual)
  write_tsv(axis_de, file.path(OUT, paste0("axis_dea_", axis_name, ".tsv")))
  all_de[[axis_name]] <- axis_de

  omnibus <- topTable(fit2, coef = c("disease", "treatment"), number = Inf, sort.by = "none") |>
    as_tibble() |>
    mutate(Accession = rownames(fit2$coefficients), .before = 1) |>
    left_join(annot, by = "Accession") |>
    mutate(axis = axis_name, omnibus_FDR_0_05 = adj.P.Val < 0.05) |>
    relocate(Accession, gene_symbol, description, axis)
  write_tsv(omnibus, file.path(OUT, paste0("axis_omnibus_", axis_name, ".tsv")))
  all_omnibus[[axis_name]] <- omnibus

  rev <- disease |>
    select(Accession, gene_symbol, description, disease_logFC = logFC,
           disease_P = P.Value, disease_FDR = adj.P.Val) |>
    inner_join(treatment |>
                 select(Accession, treatment_logFC = logFC,
                        treatment_P = P.Value, treatment_FDR = adj.P.Val),
               by = "Accession") |>
    mutate(
      axis = axis_name,
      formal_disease_DEP = disease_FDR < 0.05 & abs(disease_logFC) >= 0.25,
      candidate_disease_protein = disease_P < 0.05 & abs(disease_logFC) >= 0.58,
      reversal_index = if_else(abs(disease_logFC) >= 0.25,
                               -treatment_logFC / disease_logFC, NA_real_),
      reversal_class = case_when(
        is.na(reversal_index) ~ "Not estimable",
        reversal_index >= 0.80 ~ "Full reversal",
        reversal_index >= 0.30 ~ "Partial reversal",
        reversal_index >= 0 ~ "Minimal/no reversal",
        TRUE ~ "Exacerbation"
      ),
      formal_statistical_rescue = formal_disease_DEP & treatment_FDR < 0.05 &
        abs(treatment_logFC) >= 0.25 & sign(treatment_logFC) == -sign(disease_logFC)
    )
  write_tsv(rev, file.path(OUT, paste0("axis_reversal_", axis_name, ".tsv")))
  all_rev[[axis_name]] <- rev
  cand <- rev |> filter(candidate_disease_protein, is.finite(reversal_index))
  reversal_summary_rows[[axis_name]] <- tibble(
    axis = axis_name,
    n_formal_disease_DEPs = sum(rev$formal_disease_DEP, na.rm = TRUE),
    n_candidate_disease_proteins = nrow(cand),
    n_full_reversal = sum(cand$reversal_class == "Full reversal"),
    n_partial_reversal = sum(cand$reversal_class == "Partial reversal"),
    n_minimal_no_reversal = sum(cand$reversal_class == "Minimal/no reversal"),
    n_exacerbation = sum(cand$reversal_class == "Exacerbation"),
    n_formal_statistical_rescues = sum(rev$formal_statistical_rescue, na.rm = TRUE),
    median_RI_candidate_set = median(cand$reversal_index, na.rm = TRUE)
  )

  if (nrow(cand) >= 5) {
    xs <- x[cand$Accession, , drop = FALSE]
    row_med <- apply(xs, 1, median, na.rm = TRUE)
    for (j in seq_len(ncol(xs))) {
      miss <- !is.finite(xs[, j])
      xs[miss, j] <- row_med[miss]
    }
    z <- t(scale(t(xs)))
    z[!is.finite(z)] <- 0
    disease_dir <- sign(cand$disease_logFC)
    names(disease_dir) <- cand$Accession
    score <- colMeans(sweep(z, 1, disease_dir[rownames(z)], "*"), na.rm = TRUE)
    score_tbl <- tibble(sample_id = names(score), disease_axis_score = as.numeric(score)) |>
      left_join(m |> select(sample_id, group), by = "sample_id") |>
      mutate(axis = axis_name, n_axis_proteins = nrow(cand))
    axis_score_rows[[axis_name]] <- score_tbl
    untreated_group <- cfg$groups[2]
    treated_group <- cfg$groups[3]
    untreated <- score_tbl$disease_axis_score[score_tbl$group == untreated_group]
    treated <- score_tbl$disease_axis_score[score_tbl$group == treated_group]
    wt <- t.test(treated, untreated)
    axis_test_rows[[axis_name]] <- tibble(
      axis = axis_name, untreated_group = untreated_group, treated_group = treated_group,
      n_axis_proteins = nrow(cand), mean_untreated = mean(untreated),
      mean_treated = mean(treated), mean_difference = mean(treated) - mean(untreated),
      welch_p = wt$p.value, exact_permutation_p = exact_permutation_p(untreated, treated),
      selection_rule = "axis_specific_rawP<0.05_abs_log2FC>=0.58_with_3of4_coverage"
    )
  }

  for (cn in unique(axis_de$contrast)) {
    z <- axis_de |> filter(contrast == cn)
    summary_rows[[paste(axis_name, cn)]] <- tibble(
      axis = axis_name,
      contrast = cn,
      n_tested = nrow(z),
      n_rawP_0_05 = sum(z$P.Value < 0.05, na.rm = TRUE),
      n_nominal_effect = sum(z$nominal_effect_candidate, na.rm = TRUE),
      n_FDR_0_10 = sum(z$supportive_FDR, na.rm = TRUE),
      n_FDR_0_05 = sum(z$formal_DEP, na.rm = TRUE),
      min_rawP = min(z$P.Value, na.rm = TRUE),
      min_FDR = min(z$adj.P.Val, na.rm = TRUE)
    )
  }

  plot_tbl <- bind_rows(disease |> mutate(panel = "Disease vs CTRL"),
                        treatment |> mutate(panel = "Dapagliflozin vs disease")) |>
    mutate(display = case_when(formal_DEP ~ "FDR < 0.05",
                               supportive_FDR ~ "FDR < 0.10",
                               nominal_effect_candidate ~ "Nominal P + effect",
                               TRUE ~ "Not selected"))
  p <- ggplot(plot_tbl, aes(logFC, -log10(pmax(P.Value, 1e-300)), colour = display)) +
    geom_point(size = 1.15, alpha = 0.68) +
    geom_vline(xintercept = c(-0.58, 0.58), linetype = 2, colour = "grey55") +
    geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey55") +
    facet_wrap(~panel) +
    scale_colour_manual(values = c("FDR < 0.05" = "#D55E00", "FDR < 0.10" = "#E69F00",
                                   "Nominal P + effect" = "#0072B2", "Not selected" = "#BDBDBD")) +
    theme_bw(base_size = 10) + theme(legend.position = "bottom") +
    labs(title = paste(axis_name, "axis-specific limma"), x = "log2 fold change",
         y = expression(-log[10](italic(P))), colour = NULL)
  ggsave(file.path(FIG, paste0("axis_volcano_", axis_name, ".pdf")), p, width = 8.2, height = 4.8)
  ggsave(file.path(FIG, paste0("axis_volcano_", axis_name, ".png")), p, width = 8.2, height = 4.8, dpi = 400)
}

write_tsv(bind_rows(all_de), file.path(OUT, "axis_dea_all.tsv"))
write_tsv(bind_rows(all_rev), file.path(OUT, "axis_reversal_all.tsv"))
write_tsv(bind_rows(all_omnibus), file.path(OUT, "axis_omnibus_all.tsv"))
write_tsv(bind_rows(summary_rows), file.path(OUT, "axis_dea_summary.tsv"))
write_tsv(bind_rows(reversal_summary_rows), file.path(OUT, "axis_reversal_summary.tsv"))
write_tsv(bind_rows(axis_score_rows), file.path(OUT, "axis_disease_axis_sample_scores.tsv"))
write_tsv(bind_rows(axis_test_rows), file.path(OUT, "axis_disease_axis_treatment_tests.tsv"))

message("[05] independent PreDM and DM limma models complete")
