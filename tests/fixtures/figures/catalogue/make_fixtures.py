"""Synthetic figure-source tables for the R16e catalogue (V206-V214). Run from the repository root:

    .venv/bin/python tests/fixtures/figures/catalogue/make_fixtures.py

Every value comes from fixed trigonometric formulas with planted structure (shifts, missing blocks, up and down effects,
DEP clusters, overlapping sets). Nothing is study data. The PCA scores are the output of stats::prcomp on the synthetic
matrix (the display matrix of SM06: missing cells filled with the feature median, then centred); they are the figure's
source table, and the tests recompute them independently. Files are written with LF line endings.
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "src"))

from proteomics_pipeline.runtime import run_r_code  # noqa: E402


def _tsv(name: str, header: list[str], rows: list[list[object]]) -> Path:
    lines = ["\t".join(header)]
    for row in rows:
        cells = []
        for v in row:
            if isinstance(v, bool):
                cells.append("TRUE" if v else "FALSE")
            elif isinstance(v, float):
                cells.append("NA" if math.isnan(v) else repr(v))
            else:
                cells.append(str(v))
        lines.append("\t".join(cells))
    path = HERE / name
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return path


NAN = float("nan")


def _matrix_value(f: int, s: int, shift: float = 0.0) -> float:
    return 20.0 + 0.5 * math.sin(0.37 * f) + 0.4 * math.cos(0.91 * s + 0.13 * f) + shift


def qc_intensity() -> None:
    rows = []
    samples = [f"S{s:02d}" for s in range(1, 13)]
    for si, sample in enumerate(samples):
        shift = 1.5 if sample in ("S03", "S09") else 0.0
        for f in range(1, 61):
            missing = f % 7 == 0 and (si + 1) % 5 == 0
            value = NAN if missing else _matrix_value(f, si + 1, shift)
            rows.append([sample, f"F{f:03d}", value, not missing])
    _tsv("qc_intensity_distributions_source.tsv", ["sample", "feature", "value", "observed"], rows)
    # the same cells with the observed mask removed and the missing cells imputed: refused by V206 negative
    imputed = [[r[0], r[1], (r[2] if not math.isnan(r[2]) else 20.0)] for r in rows]
    _tsv("qc_intensity_distributions_imputed.tsv", ["sample", "feature", "value"], imputed)


def qc_missingness() -> None:
    rows = []
    for f in range(1, 31):
        for s in range(1, 13):
            missing = f <= 8 and s <= 4
            value = NAN if missing else _matrix_value(f, s)
            rows.append([f"F{f:03d}", f"S{s:02d}", not missing, value])
    _tsv("qc_missingness_heatmap_source.tsv", ["feature", "sample", "observed", "value"], rows)


def qc_correlation() -> None:
    rows = []
    for s in range(1, 11):
        block = 0 if s <= 5 else 1
        for f in range(1, 51):
            missing = (f + s) % 9 == 0
            base = 20.0 + 0.4 * math.sin(0.5 * f) + (1.2 * math.cos(0.7 * f) if block else 0.0)
            value = NAN if missing else base + 0.3 * math.sin(1.1 * s + 0.2 * f)
            rows.append([f"S{s:02d}", f"F{f:03d}", not missing, value])
    _tsv("qc_sample_correlation_source.tsv", ["sample", "feature", "observed", "value"], rows)
    _tsv("qc_sample_correlation_imputed.tsv", ["sample", "feature", "value"],
         [[r[0], r[1], r[3] if not math.isnan(r[3]) else 20.0] for r in rows])


def qc_pca() -> None:
    """Matrix (features x samples) for the PCA figure, and the scores from prcomp on the SM06 display matrix."""
    samples = [f"S{s:02d}" for s in range(1, 17)]
    matrix_rows = []
    for f in range(1, 41):
        row = []
        for s in range(1, 17):
            shift = 0.0 if s <= 8 else 1.1 * math.sin(0.3 * f)
            missing = (f * 7 + s) % 11 == 0
            row.append(NAN if missing else _matrix_value(f, s, shift))
        matrix_rows.append([f"F{f:03d}"] + row)
    _tsv("qc_pca_matrix.tsv", ["feature"] + samples, matrix_rows)
    r_code = (
        "args <- commandArgs(trailingOnly = TRUE)\n"
        "m <- read.delim(args[1], check.names = FALSE, row.names = 1)\n"
        "m <- as.matrix(m)\n"
        "filled <- t(apply(m, 1, function(x) { x[is.na(x)] <- median(x, na.rm = TRUE); x }))\n"
        "p <- stats::prcomp(t(filled), center = TRUE, scale. = FALSE)\n"
        "scores <- p$x[, 1:2]\n"
        "ve <- 100 * p$sdev^2 / sum(p$sdev^2)\n"
        "write.table(data.frame(sample = rownames(scores), PC1 = scores[, 1], PC2 = scores[, 2], PC1_var_pct = ve[1], PC2_var_pct = ve[2]), args[2], sep = '\\t', quote = FALSE, row.names = FALSE)\n"
    )
    scores_path = HERE / "qc_pca_scores_prcomp.tsv"
    result = run_r_code(r_code, [HERE / "qc_pca_matrix.tsv", scores_path])
    if result.returncode != 0:
        raise SystemExit(f"prcomp failed: {result.stderr}")
    text = scores_path.read_text(encoding="utf-8").splitlines()
    scores = {line.split("\t")[0]: line.split("\t")[1:] for line in text[1:]}
    variance = (float(scores[samples[0]][2]), float(scores[samples[0]][3]))   # the same prcomp, variance explained in per cent
    group_main = {s: ("G1" if i < 8 else "G2") for i, s in enumerate(samples)}
    group_variant = {s: ("G1" if i < 8 else ("G2" if i < 14 else "G3")) for i, s in enumerate(samples)}
    for name, groups in (("qc_pca_ellipses_source.tsv", group_main), ("qc_pca_ellipses_variant_source.tsv", group_variant)):
        rows = [[s, groups[s], float(scores[s][0]), float(scores[s][1]), variance[0], variance[1]] for s in samples]
        _tsv(name, ["sample", "group", "PC1", "PC2", "PC1_var_pct", "PC2_var_pct"], rows)
    scores_path.unlink()


def volcano() -> None:
    rows = []
    for f in range(1, 201):
        effect = 2.0 * math.sin(0.61 * f) * 0.5
        if f <= 30:
            effect, p = 2.0 + 0.1 * math.sin(f), 10 ** -(3.5 + 0.02 * f)
        elif f <= 60:
            effect, p = -2.0 + 0.1 * math.cos(f), 10 ** -(3.2 + 0.02 * f)
        else:
            p = 10 ** -(0.2 + 2.0 * abs(math.sin(0.23 * f)))
        q = min(1.0, p * 200 / (1 + f % 37))
        rows.append([f"F{f:03d}", effect, p, q, 6, 6])
    _tsv("diff_volcano_labelled_source.tsv", ["feature", "effect", "P", "q", "n_Control", "n_Acute"], rows)
    labels = [f"F{f:03d}" for f in (1, 2, 3, 4, 5, 31, 32, 33, 34, 35)]
    _tsv("diff_volcano_labels.tsv", ["feature"], [[x] for x in labels])
    _tsv("diff_volcano_labels_absent.tsv", ["feature"], [[x] for x in labels[:9] + ["F999"]])


def ma_plot() -> None:
    rows = []
    for contrast, numerator, denominator in (("Acute - Control", "Acute", "Control"), ("Chronic - Control", "Chronic", "Control")):
        for f in range(1, 201):
            mean_den = 18.0 + 6.0 * abs(math.sin(0.11 * f))
            known = 1.5 * math.sin(0.4 * f) if contrast.startswith("Acute") else -1.0 * math.cos(0.3 * f)
            mean_num = mean_den + known
            p = 10 ** -(1.0 + abs(known))
            rows.append([contrast, numerator, denominator, f"F{f:03d}", mean_num, mean_den, p, min(1.0, p * 200), 6, 6])
    _tsv("diff_ma_plot_source.tsv", ["contrast", "numerator", "denominator", "feature", "mean_numerator", "mean_denominator", "P", "q", "n_numerator", "n_denominator"], rows)


GROUPS3 = ["Control", "Acute", "Chronic"]


def dotplot_brackets() -> None:
    values, family = [], []
    for k, protein in enumerate(("P1", "P2", "P3")):
        means = {}
        for g, group in enumerate(GROUPS3):
            vals = [18.0 + 1.2 * k + 1.6 * g + 0.35 * math.sin(1.3 * r + 0.7 * g + k) for r in range(5)]
            means[group] = vals
            for r, v in enumerate(vals):
                values.append([protein, group, r + 1, v])
        for a, b in (("Control", "Acute"), ("Control", "Chronic"), ("Acute", "Chronic")):
            diff = sum(means[b]) / 5 - sum(means[a]) / 5
            p = 10 ** -(1.0 + 0.8 * abs(diff))
            family.append([protein, f"{b} - {a}", a, b, diff, p, min(1.0, p * 9), 5, 5, 5])
    _tsv("diff_dotplot_values.tsv", ["protein", "group", "replicate", "value"], values)
    _tsv("diff_dotplot_family.tsv", ["protein", "comparison", "group_a", "group_b", "effect", "P", "q", "n_Control", "n_Acute", "n_Chronic"], family)


def dep_heatmap() -> None:
    samples = [f"S{s:02d}" for s in range(1, 13)]
    rows = []
    for f in range(1, 41):
        cluster = 0 if f <= 20 else 1
        for si, s in enumerate(samples):
            shift = (1.8 if (cluster == 0) == (si < 6) else -1.8)
            rows.append([f"D{f:03d}", s, 20.0 + shift + 0.2 * math.sin(0.9 * f + si), True, "Control" if si < 6 else "Acute"])
    _tsv("diff_dep_heatmap_source.tsv", ["feature", "sample", "value", "declared", "group"], rows)
    non_dep = [[f"N{f:03d}", s, 20.0 + 0.3 * math.sin(0.5 * f + si), False, "Control" if si < 6 else "Acute"] for f in range(1, 41) for si, s in enumerate(samples)]
    _tsv("diff_dep_heatmap_with_nondep.tsv", ["feature", "sample", "value", "declared", "group"], rows + non_dep)


def upset_venn() -> None:
    k3 = []
    for f in range(1, 121):
        k3.append([f"F{f:03d}", int(f % 2 == 0 or f % 5 == 0), int(f % 3 == 0 or f % 7 == 0), int(f % 4 == 0 or f % 11 == 0)])
    _tsv("diff_upset_venn_k3.tsv", ["feature", "A", "B", "C"], k3)
    k5 = []
    for f in range(1, 151):
        k5.append([f"F{f:03d}", int(f % 2 == 0), int(f % 3 == 0), int(f % 5 == 0), int(f % 7 == 0), int(f % 11 < 4)])
    _tsv("diff_upset_venn_k5.tsv", ["feature", "A", "B", "C", "D", "E"], k5)


# ---------------------------------------------------------------- stage 2 (R16e, V215-V222)
#
# Biomarker, pathway, network and PERMANOVA figures. The analysis values of each source (ROC curve and bootstrap band,
# Welch intervals, permutation P, hypergeometric P and BH q, module eigengenes, principal coordinates, PERMDISP P) are
# computed here with the R packages that the pipeline uses (pROC, stats, vegan), and written as the figure's source
# table. The R renderers draw only from those tables; the tests recompute the same values independently.

CV_LABEL = "cross_validated_nested"     # the frozen claim label of cross-validated results (post_de_biomarker.R)
ROC_SEED = 20261008
ROC_BOOTSTRAPS = 2000
PERMDISP_SEED = 20261003
PERMDISP_B = 199


def _r_json(code: str, args: list) -> dict:
    result = run_r_code(code, args, cwd=HERE)
    if result.returncode != 0:
        raise SystemExit(f"R fixture step failed: {result.stderr[-600:]}")
    return json.loads(result.stdout.strip().splitlines()[-1])


def bio_roc() -> None:
    """V215: 15 controls and 15 cases with planted separation; out-of-fold scores (the analysis input), then the ROC curve,
    the stratified bootstrap band (pROC ci.se, seed recorded) and the AUC (the figure source)."""
    scores = []
    for i in range(15):
        scores.append([f"C{i + 1:02d}", "Control", round(0.8 * math.sin(1.7 * i + 0.4), 6), "out_of_fold"])
        scores.append([f"K{i + 1:02d}", "Case", round(1.0 + 0.9 * math.sin(1.3 * i + 1.1), 6), "out_of_fold"])
    _tsv("bio_roc_scores.tsv", ["unit_id", "class", "score", "score_source"], scores)
    code = (
        "args <- commandArgs(trailingOnly = TRUE)\n"
        "suppressPackageStartupMessages(library(pROC))\n"
        "d <- read.delim(args[1], stringsAsFactors = FALSE, fileEncoding = 'UTF-8')\n"
        "r <- roc(response = d$class, predictor = d$score, levels = c('Control', 'Case'), direction = '<', quiet = TRUE)\n"
        "ci <- as.numeric(ci.auc(r, method = 'delong'))\n"
        "sp_grid <- seq(0, 1, 0.02)\n"
        f"set.seed({ROC_SEED})\n"
        f"se <- ci.se(r, specificities = sp_grid, boot.n = {ROC_BOOTSTRAPS})\n"
        "curve <- data.frame(fpr = 1 - r$specificities, sensitivity = r$sensitivities)\n"
        "curve <- curve[order(curve$fpr, curve$sensitivity), ]\n"
        "out <- list(auc = as.numeric(pROC::auc(r)), auc_lower = ci[1], auc_upper = ci[3],\n"
        "  curve_fpr = curve$fpr, curve_sensitivity = curve$sensitivity,\n"
        "  band_fpr = 1 - sp_grid, band_lower = se[, 1], band_upper = se[, 3])\n"
        "cat(jsonlite::toJSON(out, auto_unbox = TRUE, digits = NA))\n"
    )
    result = _r_json(code, [HERE / "bio_roc_scores.tsv"])
    n_control, n_case = 15, 15
    header = ["component", "fpr", "sensitivity", "band_lower", "band_upper", "auc", "auc_lower", "auc_upper", "n_control", "n_case",
              "claim_label", "score_source", "seed", "B"]
    common = [result["auc"], result["auc_lower"], result["auc_upper"], n_control, n_case, CV_LABEL, "out_of_fold", ROC_SEED, ROC_BOOTSTRAPS]
    rows = [["curve", f, s, NAN, NAN] + common for f, s in zip(result["curve_fpr"], result["curve_sensitivity"])]
    rows += [["band", f, NAN, lo, hi] + common for f, lo, hi in zip(result["band_fpr"], result["band_lower"], result["band_upper"])]
    _tsv("bio_roc_ci_band_source.tsv", header, rows)
    # the same figure from a training-fold (in-sample) score set: refused by V215 (E_BIOMARKER_CLAIM)
    training = [list(r) for r in rows]
    for r in training:
        r[header.index("claim_label")] = "in_sample"
        r[header.index("score_source")] = "training_fold"
    _tsv("bio_roc_ci_band_training_source.tsv", header, training)


PANEL_FEATURES = ["PANEL1", "PANEL2", "PANEL3", "PANEL4", "PANEL5"]
PANEL_EFFECTS = [1.2, 0.8, 0.0, -0.6, 0.4]


def bio_panel_forest() -> None:
    """V216: five panel features, ten units per group; Welch intervals and P from stats::t.test, BH q across the panel."""
    values = []
    for f, feature in enumerate(PANEL_FEATURES):
        for g, group in enumerate(("Control", "Case")):
            for u in range(10):
                value = 20.0 + 0.3 * f + (PANEL_EFFECTS[f] if group == "Case" else 0.0) + 0.4 * math.sin(0.7 * u + f + 1.3 * g)
                values.append([feature, f"{group[0]}{u + 1:02d}", group, value])
    _tsv("bio_panel_values.tsv", ["feature", "unit_id", "group", "value"], values)
    code = (
        "args <- commandArgs(trailingOnly = TRUE)\n"
        "d <- read.delim(args[1], stringsAsFactors = FALSE, fileEncoding = 'UTF-8')\n"
        "feats <- unique(d$feature)\n"
        "res <- lapply(feats, function(f) {\n"
        "  x <- d[d$feature == f, ]\n"
        "  t <- stats::t.test(x$value[x$group == 'Case'], x$value[x$group == 'Control'])\n"
        "  data.frame(feature = f, estimate = unname(t$estimate[1] - t$estimate[2]), ci_lower = t$conf.int[1], ci_upper = t$conf.int[2],\n"
        "    P = t$p.value, n_Control = sum(x$group == 'Control'), n_Case = sum(x$group == 'Case'))\n"
        "})\n"
        "res <- do.call(rbind, res)\n"
        "res$q <- stats::p.adjust(res$P, method = 'BH')\n"
        "cat(jsonlite::toJSON(as.list(res), auto_unbox = FALSE, digits = NA))\n"
    )
    result = _r_json(code, [HERE / "bio_panel_values.tsv"])
    header = ["feature", "estimate", "ci_lower", "ci_upper", "P", "q", "n_Control", "n_Case", "claim_label"]
    rows = [[result["feature"][i], result["estimate"][i], result["ci_lower"][i], result["ci_upper"][i], result["P"][i], result["q"][i],
             result["n_Control"][i], result["n_Case"][i], CV_LABEL] for i in range(len(PANEL_FEATURES))]
    _tsv("bio_panel_forest_source.tsv", header, rows)
    _tsv("bio_panel_forest_unlabelled_source.tsv", header, [r[:-1] + [""] for r in rows])


CV_REPEATS = [0.84, 0.79, 0.86, 0.81, 0.83]


def bio_cv_auc() -> None:
    """V217: five nested-CV repeat AUCs and a permutation null of B = 19 AUCs (planted near 0.5). The displayed P is
    (k + 1) / (B + 1) with k the null values at or above the observed mean of the repeats (0 here: P = 0.05)."""
    null = [0.5 + 0.09 * math.sin(2.3 * j + 0.7) for j in range(19)]
    observed = sum(CV_REPEATS) / len(CV_REPEATS)
    k = sum(1 for v in null if v >= observed)
    p_display = (k + 1) / (19 + 1)
    header = ["series", "value", "B", "seed", "P_display", "claim_label"]
    rows = [["cv_repeat", v, 19, 20261009, p_display, CV_LABEL] for v in CV_REPEATS]
    rows += [["null", v, 19, 20261009, p_display, CV_LABEL] for v in null]
    _tsv("bio_cv_auc_vs_permutation_source.tsv", header, rows)
    _tsv("bio_cv_auc_vs_permutation_wrong_p_source.tsv", header, [r[:4] + [0.1] + r[5:] for r in rows])


def main() -> None:
    qc_intensity()
    qc_missingness()
    qc_correlation()
    qc_pca()
    volcano()
    ma_plot()
    dotplot_brackets()
    dep_heatmap()
    upset_venn()
    bio_roc()
    bio_panel_forest()
    bio_cv_auc()


if __name__ == "__main__":
    main()
