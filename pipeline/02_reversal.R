#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
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
FIG <- file.path(ROOT, "results", "figures", "reversal")
REPORT <- file.path(ROOT, "report")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

meta <- read_csv(file.path(PROC, "sample_meta.csv"), show_col_types = FALSE)
meta$group <- factor(meta$group, levels=c("CTRL","PreDM","DM","PreDM_Dapa","DM_Dapa"))
mat_tbl <- read_tsv(file.path(PROC, "abundance_analysis_log2.tsv"), show_col_types = FALSE)
mat <- as.matrix(mat_tbl[, meta$sample_id]); storage.mode(mat) <- "numeric"; rownames(mat) <- mat_tbl$Accession
annot <- read_tsv(file.path(PROC, "protein_annot.tsv"), show_col_types = FALSE)

read_dea <- function(cn) read_tsv(file.path(TAB, "dea", paste0("dea_", cn, ".tsv")), show_col_types = FALSE)
axes <- list(
  PreDM = list(disease="D1_PreDM_vs_CTRL", treatment="T1_PreDM_Dapa_vs_PreDM", untreated="PreDM", treated="PreDM_Dapa"),
  DM = list(disease="D2_DM_vs_CTRL", treatment="T2_DM_Dapa_vs_DM", untreated="DM", treated="DM_Dapa")
)

exact_permutation_p <- function(x, y) {
  observed <- mean(y) - mean(x)
  pooled <- c(x, y)
  n_x <- length(x)
  assignments <- combn(seq_along(pooled), n_x)
  null <- apply(assignments, 2, function(idx) mean(pooled[-idx]) - mean(pooled[idx]))
  (sum(abs(null) >= abs(observed)) + 1) / (length(null) + 1)
}

rescue_tables <- list()
rescue_summary <- list()
axis_scores_all <- list()
axis_tests <- list()

for (axis_name in names(axes)) {
  axis_cfg <- axes[[axis_name]]
  disease <- read_dea(axis_cfg$disease) |>
    select(Accession, gene_symbol, description, disease_logFC=logFC, disease_P=P.Value, disease_FDR=adj.P.Val)
  treatment <- read_dea(axis_cfg$treatment) |>
    select(Accession, treatment_logFC=logFC, treatment_P=P.Value, treatment_FDR=adj.P.Val)
  joined <- inner_join(disease, treatment, by="Accession") |>
    mutate(
      axis = axis_name,
      formal_disease_DEP = disease_FDR < 0.05 & abs(disease_logFC) >= 0.25,
      candidate_disease_protein = disease_P < 0.05 & abs(disease_logFC) >= 0.58,
      reversal_direction = sign(treatment_logFC) == -sign(disease_logFC),
      reversal_index = if_else(abs(disease_logFC) >= 0.25, -treatment_logFC / disease_logFC, NA_real_),
      effect_size_reversal = candidate_disease_protein & reversal_direction & reversal_index >= 0.30,
      formal_statistical_rescue = formal_disease_DEP & treatment_FDR < 0.05 & abs(treatment_logFC) >= 0.25 & reversal_direction,
      evidence_tier = case_when(
        formal_statistical_rescue ~ "formal_statistical_rescue",
        effect_size_reversal ~ "candidate_effect_size_reversal",
        candidate_disease_protein ~ "candidate_disease_not_reversed",
        TRUE ~ "not_in_candidate_disease_set"
      )
    )
  rescue_tables[[axis_name]] <- joined
  write_tsv(joined, file.path(TAB, paste0("reversal_", axis_name, ".tsv")))
  rescue_summary[[axis_name]] <- tibble(
    axis = axis_name,
    n_formal_disease_DEPs = sum(joined$formal_disease_DEP, na.rm=TRUE),
    n_candidate_disease_proteins = sum(joined$candidate_disease_protein, na.rm=TRUE),
    n_candidate_effect_size_reversals = sum(joined$effect_size_reversal, na.rm=TRUE),
    pct_candidate_effect_size_reversal = ifelse(sum(joined$candidate_disease_protein)>0, 100*sum(joined$effect_size_reversal)/sum(joined$candidate_disease_protein), NA_real_),
    n_formal_statistical_rescues = sum(joined$formal_statistical_rescue, na.rm=TRUE),
    median_RI_candidate_set = median(joined$reversal_index[joined$candidate_disease_protein], na.rm=TRUE)
  )

  p_scatter <- ggplot(joined, aes(disease_logFC, treatment_logFC, colour=evidence_tier)) +
    geom_hline(yintercept=0, colour="grey70") + geom_vline(xintercept=0, colour="grey70") +
    geom_abline(slope=-0.3, intercept=0, linetype=2, colour="#2166AC") +
    geom_point(alpha=0.65, size=1.2) +
    scale_colour_manual(values=c(formal_statistical_rescue="#B2182B", candidate_effect_size_reversal="#2166AC", candidate_disease_not_reversed="#E69F00", not_in_candidate_disease_set="#BDBDBD")) +
    theme_bw(base_size=10) +
    labs(title=paste(axis_name, "disease and dapagliflozin effects"), x=paste(axis_name, "vs CTRL log2FC"), y=paste(axis_name, "+ dapagliflozin vs", axis_name, "log2FC"))
  ggsave(file.path(FIG, paste0("reversal_scatter_", axis_name, ".pdf")), p_scatter, width=7, height=5.5)

  # Rat-level continuous disease-axis scores. Candidate selection is explicitly exploratory.
  selected <- joined |> filter(candidate_disease_protein) |> arrange(disease_P)
  selection_rule <- "rawP<0.05_and_abs_log2FC>=0.58"
  if (nrow(selected) < 5) {
    selected <- joined |> arrange(disease_P) |> slice_head(n=min(50, n()))
    selection_rule <- "fallback_top50_by_disease_rawP"
  }
  selected <- selected |> filter(Accession %in% rownames(mat))
  x <- mat[selected$Accession, , drop=FALSE]
  row_med <- rowMedians(x, na.rm=TRUE)
  for (j in seq_len(ncol(x))) x[is.na(x[,j]),j] <- row_med[is.na(x[,j])]
  z <- t(scale(t(x))); z[!is.finite(z)] <- 0
  direction <- sign(selected$disease_logFC); names(direction) <- selected$Accession
  score <- colMeans(sweep(z, 1, direction[rownames(z)], "*"), na.rm=TRUE)
  score_tbl <- tibble(sample_id=names(score), disease_axis_score=as.numeric(score)) |>
    left_join(meta |> select(sample_id, group), by="sample_id") |>
    mutate(axis=axis_name, selection_rule=selection_rule, n_axis_proteins=nrow(selected))
  axis_scores_all[[axis_name]] <- score_tbl
  untreated_scores <- score_tbl$disease_axis_score[score_tbl$group == axis_cfg$untreated]
  treated_scores <- score_tbl$disease_axis_score[score_tbl$group == axis_cfg$treated]
  welch <- t.test(treated_scores, untreated_scores)
  axis_tests[[axis_name]] <- tibble(
    axis=axis_name,
    untreated_group=axis_cfg$untreated,
    treated_group=axis_cfg$treated,
    n_axis_proteins=nrow(selected),
    selection_rule=selection_rule,
    mean_untreated=mean(untreated_scores),
    mean_treated=mean(treated_scores),
    mean_difference=mean(treated_scores)-mean(untreated_scores),
    welch_p=welch$p.value,
    exact_permutation_p=exact_permutation_p(untreated_scores, treated_scores)
  )
}

summary_tbl <- bind_rows(rescue_summary)
scores_tbl <- bind_rows(axis_scores_all)
tests_tbl <- bind_rows(axis_tests)
write_tsv(summary_tbl, file.path(TAB, "reversal_summary.tsv"))
write_tsv(scores_tbl, file.path(TAB, "disease_axis_sample_scores.tsv"))
write_tsv(tests_tbl, file.path(TAB, "disease_axis_treatment_tests.tsv"))

p_scores <- ggplot(scores_tbl, aes(group, disease_axis_score, colour=group)) +
  geom_boxplot(outlier.shape=NA, alpha=0.15) + geom_jitter(width=0.08, size=2.5) +
  facet_wrap(~axis, scales="free_x") + theme_bw(base_size=10) +
  theme(axis.text.x=element_text(angle=35, hjust=1), legend.position="none") +
  labs(title="Exploratory disease-axis scores", x=NULL, y="Direction-weighted standardized abundance")
ggsave(file.path(FIG, "disease_axis_scores.pdf"), p_scores, width=9, height=5.5)

notes <- c(
  "# Dapagliflozin Reversal Analysis", "",
  "Two independent disease axes were evaluated: PreDM and DM.", "",
  "- Formal disease DEP: FDR < 0.05 and |log2FC| >= 0.25.",
  "- Candidate disease protein: raw P < 0.05 and |log2FC| >= 0.58.",
  "- Candidate effect-size reversal: opposing dapagliflozin effect with RI >= 0.30.",
  "- Formal statistical rescue additionally requires FDR-supported disease and treatment effects.",
  "- Candidate sets and disease-axis scores are exploratory because selection uses raw P when formal disease sets are empty.",
  "- Reversal indices and axis scores may be influenced by regression to the mean; direct treatment-vs-disease contrasts remain the primary treatment tests.", "",
  paste(capture.output(print(summary_tbl)), collapse="\n"), "",
  paste(capture.output(print(tests_tbl)), collapse="\n")
)
writeLines(notes, file.path(REPORT, "reversal_notes.md"))
message("[02] stage-specific dapagliflozin reversal analysis complete")
