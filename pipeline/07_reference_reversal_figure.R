#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

set.seed(48)
ROOT <- Sys.getenv("DAPA_ANALYSIS_ROOT", unset = "")
if (!nzchar(ROOT)) stop("DAPA_ANALYSIS_ROOT must be set")
ROOT <- normalizePath(ROOT, mustWork = TRUE)
TAB <- file.path(ROOT, "results", "tables", "axis_specific")
FIG <- file.path(ROOT, "results", "figures", "axis_specific")
SRC <- file.path(TAB, "figure_source_data")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(SRC, recursive = TRUE, showWarnings = FALSE)

AXES <- c("PreDM", "DM")
axis_cols <- c(PreDM = "#0072B2", DM = "#CC79A7")
class_levels <- c("Exacerbation", "Minimal/no reversal", "Partial reversal", "Full reversal", "Not estimable")
class_cols <- c("Exacerbation" = "#3B0F70", "Minimal/no reversal" = "#D95F59",
                "Partial reversal" = "#F28E2B", "Full reversal" = "#F0C808",
                "Not estimable" = "#D9D9D9")

rev <- bind_rows(lapply(AXES, function(a) {
  read_tsv(file.path(TAB, paste0("axis_reversal_", a, ".tsv")), show_col_types = FALSE)
}))
candidate <- rev |>
  filter(candidate_disease_protein, is.finite(reversal_index)) |>
  mutate(axis = factor(axis, levels = AXES))
write_tsv(candidate, file.path(SRC, "panel_A_candidate_reversal_indices.tsv"))

selected <- bind_rows(lapply(AXES, function(a) {
  candidate |>
    filter(axis == a) |>
    arrange(disease_P, desc(abs(disease_logFC))) |>
    slice_head(n = 25) |>
    transmute(Accession, selected_from_axis = a, selection_disease_P = disease_P)
})) |>
  distinct(Accession, .keep_all = TRUE)
selected_acc <- selected$Accession

gene_labels <- rev |>
  filter(Accession %in% selected_acc) |>
  distinct(Accession, .keep_all = TRUE) |>
  transmute(Accession, label = ifelse(is.na(gene_symbol) | gene_symbol == "", Accession, gene_symbol))
gene_label_map <- setNames(gene_labels$label, gene_labels$Accession)

heat <- expand_grid(axis = AXES, Accession = selected_acc) |>
  left_join(rev |>
              select(axis, Accession, gene_symbol, description, disease_logFC, disease_P,
                     disease_FDR, treatment_logFC, treatment_P, treatment_FDR,
                     reversal_index, reversal_class, candidate_disease_protein),
            by = c("axis", "Accession")) |>
  left_join(selected, by = "Accession") |>
  mutate(
    axis = factor(axis, levels = base::rev(AXES)),
    Accession = factor(Accession, levels = selected_acc),
    reversal_class = factor(ifelse(is.na(reversal_class), "Not estimable", reversal_class),
                            levels = class_levels)
  )
write_tsv(heat |> mutate(axis = as.character(axis), Accession = as.character(Accession)),
          file.path(SRC, "panel_B_selected_protein_reversal.tsv"))

gsea <- read_tsv(file.path(TAB, "axis_gsea_all.tsv"), show_col_types = FALSE) |>
  mutate(key = paste(collection, pathway, sep = "::"))
direct_treatment <- c("PreDM_Dapa_vs_PreDM", "DM_Dapa_vs_DM")
chosen <- union(gsea$key[gsea$padj < 0.05],
                gsea$key[gsea$contrast %in% direct_treatment & gsea$padj < 0.10])
path_rank <- gsea |>
  group_by(key) |>
  summarise(min_fdr = min(padj, na.rm = TRUE), max_abs_nes = max(abs(NES), na.rm = TRUE), .groups = "drop") |>
  arrange(min_fdr, desc(max_abs_nes))
if (length(chosen) < 20) chosen <- union(chosen, head(path_rank$key[!path_rank$key %in% chosen], 20 - length(chosen)))
path_keys <- path_rank |> filter(key %in% chosen) |> slice_head(n = 24) |> pull(key)

pretty_path <- function(collection, pathway) {
  x <- sub("^(REACTOME_|KEGG_|HALLMARK_)", "", pathway)
  x <- tools::toTitleCase(tolower(gsub("_", " ", x)))
  prefix <- recode(collection, KEGG_LEGACY = "KEGG", REACTOME = "Reactome", HALLMARK = "Hallmark")
  paste0(prefix, ": ", x)
}
contrast_levels <- c("PreDM_disease_vs_CTRL", "DM_disease_vs_CTRL",
                     "PreDM_Dapa_vs_PreDM", "DM_Dapa_vs_DM")
contrast_labels <- c("PreDM vs\nCTRL", "DM vs\nCTRL", "PreDM+Dapa\nvs PreDM", "DM+Dapa\nvs DM")
path <- gsea |>
  filter(key %in% path_keys, contrast %in% contrast_levels) |>
  mutate(
    contrast = factor(contrast, levels = contrast_levels),
    pathway_label = pretty_path(collection, pathway),
    pathway_label = factor(pathway_label,
                           levels = base::rev(pretty_path(
                             sub("::.*", "", path_keys), sub("^[^:]+::", "", path_keys)
                           ))),
    mark = case_when(padj < 0.05 ~ "*", padj < 0.10 ~ "·", TRUE ~ "")
  )
write_tsv(path |> mutate(contrast = as.character(contrast), pathway_label = as.character(pathway_label)),
          file.path(SRC, "panel_C_selected_pathway_GSEA.tsv"))

base_theme <- theme_classic(base_size = 9) +
  theme(text = element_text(family = "Helvetica"), plot.title = element_text(face = "bold", size = 10))

# A: protein-level RI distribution. No inferential star is shown because proteins are not animals.
sum_a <- candidate |>
  group_by(axis) |>
  summarise(median = median(reversal_index), q1 = quantile(reversal_index, 0.25),
            q3 = quantile(reversal_index, 0.75), n = n(), .groups = "drop")
p_a <- ggplot(candidate, aes(axis, reversal_index, colour = axis)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_hline(yintercept = 1, linetype = 3, colour = "grey50") +
  geom_jitter(width = 0.13, height = 0, size = 1.9, alpha = 0.82) +
  geom_linerange(data = sum_a, aes(x = axis, ymin = q1, ymax = q3), inherit.aes = FALSE,
                 linewidth = 0.65, colour = "black") +
  geom_crossbar(data = sum_a, aes(x = axis, y = median, ymin = median, ymax = median),
                inherit.aes = FALSE, width = 0.40,
                linewidth = 0.75, colour = "black") +
  geom_text(data = sum_a, aes(x = axis, y = 2.10, label = paste0("n=", n)), inherit.aes = FALSE,
            size = 2.6) +
  scale_colour_manual(values = axis_cols, guide = "none") +
  scale_x_discrete(labels = c(PreDM = "PreDM+Dapa", DM = "DM+Dapa")) +
  coord_cartesian(ylim = c(-0.8, 2.2)) +
  labs(title = "A  Protein-level reversal distributions", x = NULL, y = "Reversal index",
       caption = "Median and IQR; exploratory nominal disease sets") +
  base_theme + theme(axis.text.x = element_text(angle = 28, hjust = 1),
                     plot.caption = element_text(size = 6.5, hjust = 0))

# B: categorical strips and continuous RI heatmap.
p_b_class <- ggplot(heat, aes(Accession, axis, fill = reversal_class)) +
  geom_tile(colour = "white", linewidth = 0.15) +
  scale_fill_manual(values = class_cols, drop = FALSE, name = "Reversal class") +
  scale_x_discrete(labels = gene_label_map, expand = c(0, 0)) +
  scale_y_discrete(labels = c(DM = "DM+Dapa", PreDM = "PreDM+Dapa"), expand = c(0, 0)) +
  labs(title = "B  Protein reversal classification", x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 10),
        legend.position = "right", legend.key.height = grid::unit(3.2, "mm"))

p_b_ri <- ggplot(heat, aes(Accession, axis, fill = reversal_index)) +
  geom_tile(colour = "white", linewidth = 0.15) +
  scale_fill_viridis_c(option = "D", limits = c(-0.75, 1.50), oob = squish,
                       na.value = "#D9D9D9", name = "Reversal index") +
  scale_x_discrete(labels = gene_label_map, expand = c(0, 0)) +
  scale_y_discrete(labels = c(DM = "DM+Dapa", PreDM = "PreDM+Dapa"), expand = c(0, 0)) +
  labs(x = "Protein gene name", y = NULL) +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        axis.ticks = element_blank(), panel.grid = element_blank(), legend.position = "right")
p_b <- p_b_class / p_b_ri + plot_layout(heights = c(0.30, 0.70))

# C: pathway NES across both disease axes and direct treatment effects.
p_c <- ggplot(path, aes(contrast, pathway_label, fill = NES)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = mark, colour = abs(NES) > 1.4), size = 3.2, fontface = "bold") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  scale_fill_gradient2(low = "#5E3C99", mid = "#F7F7F7", high = "#E66101",
                       midpoint = 0, limits = c(-2.5, 2.5), oob = squish, name = "NES") +
  scale_x_discrete(labels = contrast_labels, expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0), labels = function(x) stringr::str_wrap(x, width = 58)) +
  labs(title = "C  Rank-based pathway enrichment", x = NULL, y = NULL,
       caption = "* FDR < 0.05; · FDR < 0.10") +
  theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", size = 10),
        axis.text.y = element_text(size = 6.3), axis.text.x = element_text(size = 7),
        plot.caption = element_text(size = 7, hjust = 1), legend.position = "right")

top <- p_a | p_b
composite <- top / p_c + plot_layout(heights = c(0.54, 0.46)) +
  plot_annotation(title = "Stage-specific proteomic response to dapagliflozin",
                  theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))

base <- file.path(FIG, "stage_specific_reversal_pathway_composite")
ggsave(paste0(base, ".pdf"), composite, width = 16, height = 10.5, device = cairo_pdf)
ggsave(paste0(base, ".png"), composite, width = 16, height = 10.5, dpi = 400)
ggsave(paste0(base, ".tiff"), composite, width = 16, height = 10.5, dpi = 400,
       compression = "lzw")

methods <- c(
  "# Composite figure methods", "",
  "Panel A shows per-protein reversal indices for exploratory nominal disease sets (raw P < 0.05 and |log2FC| >= 0.58). The horizontal bar is the median and the vertical segment is the interquartile range. Protein observations are not independent biological replicates, so no protein-level significance stars are shown.", "",
  "Panel B selects the 25 strongest nominal disease proteins from each axis, ordered by disease P value. RI = -(dapagliflozin vs disease log2FC)/(disease vs CTRL log2FC). Full reversal is RI >= 0.80, partial reversal is 0.30 <= RI < 0.80, minimal/no reversal is 0 <= RI < 0.30, and exacerbation is RI < 0. Values are clipped only for color display; source tables retain uncapped values.", "",
  "Panel C shows fgsea normalized enrichment scores from the two independently fitted three-group limma models. Rows include every pathway with FDR < 0.05 in any displayed contrast and pathways with FDR < 0.10 in a direct dapagliflozin contrast, up to 24 pathways. Asterisks denote FDR < 0.05 and centered dots denote FDR < 0.10."
)
writeLines(methods, file.path(SRC, "composite_figure_methods.md"))
message("[07] composite figure written: ", paste0(base, ".png"))
