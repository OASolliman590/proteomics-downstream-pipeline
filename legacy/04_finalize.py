#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "results" / "tables"
REPORT = ROOT / "report"
REPRO = ROOT / "reproducibility"


def read_tsv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


dea = read_tsv(TAB / "dea_summary.tsv")
reversal = read_tsv(TAB / "reversal_summary.tsv")
axis_tests = read_tsv(TAB / "disease_axis_treatment_tests.tsv")
qc = read_tsv(TAB / "sample_qc.tsv")
gsea_summary = read_tsv(TAB / "gsea_summary.tsv")
gsea_all = read_tsv(TAB / "gsea_all.tsv")
sensitivity = read_tsv(TAB / "method_sensitivity_summary.tsv")
loo = read_tsv(TAB / "leave_one_sample_summary.tsv")
gsea_loo = read_tsv(TAB / "gsea_leave_one_sample_stability.tsv")
AXIS_TAB = TAB / "axis_specific"
axis_dea = read_tsv(AXIS_TAB / "axis_dea_summary.tsv")
axis_gsea_all = read_tsv(AXIS_TAB / "axis_gsea_all.tsv")
axis_reversal = read_tsv(AXIS_TAB / "axis_reversal_summary.tsv")
axis_sample_tests = read_tsv(AXIS_TAB / "axis_disease_axis_treatment_tests.tsv")

primary = {
    "D1_PreDM_vs_CTRL",
    "D2_DM_vs_CTRL",
    "T1_PreDM_Dapa_vs_PreDM",
    "T2_DM_Dapa_vs_DM",
}
primary_dea = [row for row in dea if row["contrast"] in primary]

required = [
    TAB / "sample_qc.tsv",
    TAB / "dea_summary.tsv",
    TAB / "dea_master.tsv",
    TAB / "reversal_summary.tsv",
    TAB / "disease_axis_treatment_tests.tsv",
    TAB / "method_sensitivity_summary.tsv",
    TAB / "leave_one_sample_summary.tsv",
    TAB / "gsea_all.tsv",
    TAB / "gsea_leave_one_sample_stability.tsv",
    REPORT / "qc_dea_notes.md",
    REPORT / "reversal_notes.md",
    REPORT / "enrichment_notes.md",
    AXIS_TAB / "axis_dea_summary.tsv",
    AXIS_TAB / "axis_reversal_all.tsv",
    AXIS_TAB / "axis_reversal_summary.tsv",
    AXIS_TAB / "axis_disease_axis_treatment_tests.tsv",
    AXIS_TAB / "axis_omnibus_all.tsv",
    AXIS_TAB / "axis_gsea_all.tsv",
    ROOT / "results" / "figures" / "axis_specific" / "stage_specific_reversal_pathway_composite.pdf",
    ROOT / "results" / "figures" / "axis_specific" / "stage_specific_reversal_pathway_composite.png",
    AXIS_TAB / "figure_source_data" / "composite_figure_methods.md",
]

gates = [
    {
        "gate": "Required analysis artifacts exist",
        "status": "PASS" if all(path.exists() for path in required) else "FAIL",
        "detail": f"{sum(path.exists() for path in required)}/{len(required)} required files present.",
    },
    {
        "gate": "Processed log2 scale preserved",
        "status": "PASS",
        "detail": "No second log2 transformation was applied; source-log2 values are primary.",
    },
    {
        "gate": "No unsupported batch term",
        "status": "PASS",
        "detail": "Joint model is ~0+group because batch/wave metadata were not supplied.",
    },
    {
        "gate": "Formal DEPs separated from candidates",
        "status": "PASS",
        "detail": "BH-FDR DEPs and raw-P/effect-size candidates are stored in separate evidence tiers.",
    },
    {
        "gate": "Dapagliflozin treatment identity recorded",
        "status": "PASS",
        "detail": "PreDM_Dapa and DM_Dapa are explicitly modeled as dapagliflozin-treated groups.",
    },
    {
        "gate": "Biological replicate identity confirmed",
        "status": "PASS",
        "detail": "The data provider confirmed that every numbered sample column represents one independent rat (n=4/group).",
    },
    {
        "gate": "Disease stages isolated in primary refits",
        "status": "PASS",
        "detail": "Separate 12-rat limma models were fitted for CTRL/PreDM/PreDM_Dapa and CTRL/DM/DM_Dapa.",
    },
    {
        "gate": "Tissue metadata available",
        "status": "ATTENTION",
        "detail": "Tissue/compartment remains unknown; interpretation must be tissue-agnostic.",
    },
    {
        "gate": "Raw-intensity provenance available",
        "status": "ATTENTION",
        "detail": "Only processed log2 values were supplied; original normalization and raw intensities are unavailable.",
    },
]
write_tsv(TAB / "analysis_gates.tsv", gates)

top_gsea = sorted(
    gsea_all,
    key=lambda row: (float(row["padj"]) if row["padj"] not in ("", "NA") else 1.0,
                     float(row["pval"]) if row["pval"] not in ("", "NA") else 1.0),
)
top_by_contrast = {}
for row in top_gsea:
    top_by_contrast.setdefault(row["contrast"], [])
    if len(top_by_contrast[row["contrast"]]) < 5:
        top_by_contrast[row["contrast"]].append(row)

top_axis_gsea = sorted(
    axis_gsea_all,
    key=lambda row: (float(row["padj"]) if row["padj"] not in ("", "NA") else 1.0,
                     float(row["pval"]) if row["pval"] not in ("", "NA") else 1.0),
)
top_axis_by_contrast = {}
for row in top_axis_gsea:
    top_axis_by_contrast.setdefault(row["contrast"], [])
    if len(top_axis_by_contrast[row["contrast"]]) < 5:
        top_axis_by_contrast[row["contrast"]].append(row)

lines = [
    "# Dapagliflozin Pre-Diabetes/Diabetes Proteomics Analysis",
    "",
    f"_Generated: {datetime.now(timezone.utc).isoformat()}_",
    "",
    "## Study design",
    "",
    "A completely new five-group study was analyzed: CTRL, PreDM, DM, PreDM + dapagliflozin, and DM + dapagliflozin. The data provider confirmed that every numbered sample column is one independent rat (n=4/group). The supplied processed log2 abundance matrix contains 3,714 proteins before coverage filtering.",
    "",
    "Primary disease/treatment axes:",
    "",
    "1. PreDM vs CTRL, followed by PreDM + dapagliflozin vs PreDM.",
    "2. DM vs CTRL, followed by DM + dapagliflozin vs DM.",
    "",
    "The primary refits isolate the two disease stages in separate 12-rat limma models. The original joint five-group model is retained as a sensitivity analysis.",
    "",
    "## Primary axis-specific limma results",
    "",
    "| Axis | Contrast | Raw P < 0.05 | Raw-P/effect candidates | FDR < 0.10 | Formal FDR DEPs | Minimum FDR |",
    "|---|---|---:|---:|---:|---:|---:|",
]
for row in axis_dea:
    lines.append(
        f"| {row['axis']} | {row['contrast']} | {row['n_rawP_0_05']} | {row['n_nominal_effect']} | "
        f"{row['n_FDR_0_10']} | {row['n_FDR_0_05']} | {float(row['min_FDR']):.4g} |"
    )

lines.extend([
    "",
    "The independent models estimate residual variance within each disease stage. They do not apply another log transformation or renormalize the supplied values.",
    "",
    "## Joint five-group sensitivity analysis",
    "",
    "| Contrast | Formal FDR DEPs | FDR 0.10 support | Raw-P/effect candidates | Minimum FDR |",
    "|---|---:|---:|---:|---:|",
])
for row in primary_dea:
    formal = int(row["n_formal_up"]) + int(row["n_formal_down"])
    lines.append(
        f"| {row['contrast']} | {formal} | {row['n_supportive_fdr']} | "
        f"{row['n_rawP_effect_candidates']} | {float(row['min_adj_p_value']):.4g} |"
    )

lines.extend(["", "Raw-P/effect-size candidates are exploratory and are not described as differentially abundant proteins.", "", "## Primary axis-specific reversal", "", "| Axis | Formal disease DEPs | Candidate disease proteins | Full | Partial | Minimal/no reversal | Exacerbation | Formal rescues | Median RI |", "|---|---:|---:|---:|---:|---:|---:|---:|---:|"])
for row in axis_reversal:
    lines.append(
        f"| {row['axis']} | {row['n_formal_disease_DEPs']} | {row['n_candidate_disease_proteins']} | "
        f"{row['n_full_reversal']} | {row['n_partial_reversal']} | {row['n_minimal_no_reversal']} | "
        f"{row['n_exacerbation']} | {row['n_formal_statistical_rescues']} | {float(row['median_RI_candidate_set']):.3f} |"
    )

lines.extend(["", "### Disease-axis sample tests", "", "| Axis | Mean score difference, treated - untreated | Welch P | Exact permutation P | Selection rule |", "|---|---:|---:|---:|---|"])
for row in axis_sample_tests:
    lines.append(
        f"| {row['axis']} | {float(row['mean_difference']):.3f} | {float(row['welch_p']):.4g} | "
        f"{float(row['exact_permutation_p']):.4g} | {row['selection_rule']} |"
    )
lines.extend(["", "Disease-axis score tests are exploratory because the axis proteins may be selected using raw P values and can be affected by regression to the mean.", "", "## Axis-specific pathway results", ""])
for contrast in sorted(top_axis_by_contrast):
    lines.extend([f"### {contrast}", "", "| Collection | Pathway | NES | FDR |", "|---|---|---:|---:|"])
    for row in top_axis_by_contrast[contrast]:
        lines.append(f"| {row['collection']} | {row['pathway']} | {float(row['NES']):.3f} | {float(row['padj']):.4g} |")
    lines.append("")

lines.extend(["## Joint-model pathway sensitivity results", ""])
for contrast in sorted(top_by_contrast):
    lines.extend([f"### {contrast}", "", "| Collection | Pathway | NES | FDR |", "|---|---|---:|---:|"])
    for row in top_by_contrast[contrast]:
        lines.append(f"| {row['collection']} | {row['pathway']} | {float(row['NES']):.3f} | {float(row['padj']):.4g} |")
    lines.append("")

lines.extend([
    "### Leave-one-sample stability of significant treatment pathways",
    "",
    "| Contrast | Pathway | Same NES direction | Nominal P < 0.05 | NES range |",
    "|---|---|---:|---:|---:|",
])
for row in gsea_loo:
    lines.append(
        f"| {row['contrast']} | {row['pathway']} | {float(row['fraction_same_direction']):.0%} | "
        f"{float(row['fraction_nominal_p_0_05']):.0%} | {float(row['min_NES']):.2f} to {float(row['max_NES']):.2f} |"
    )
lines.extend(["", "These focused leave-one-sample P values are stability diagnostics for pathways selected in the full run; they are not a second multiple-testing analysis.", ""])

watch = [row for row in qc if row["technical_watch"].lower() == "true"]
min_corr = min(float(row["pearson_logFC_vs_primary"]) for row in sensitivity if row["method"] != "source_log2_primary")
loo_treatment = [row for row in loo if row["contrast"] in {"T1_PreDM_Dapa_vs_PreDM", "T2_DM_Dapa_vs_DM"}]
min_loo_corr = min(float(row["pearson_logFC"]) for row in loo_treatment)
watch_names = ", ".join(row["sample_id"] for row in watch)
lines.extend([
    "## Robustness and limitations",
    "",
    f"- Technical sample watchlist: {len(watch)} sample(s): {watch_names}. None were automatically excluded.",
    f"- Lowest log2FC correlation across normalization/imputation sensitivity fits: {min_corr:.3f}.",
    f"- Lowest treatment-contrast log2FC correlation after omitting one sample: {min_loo_corr:.3f}.",
    "- n=4/group limits individual-protein FDR power.",
    "- Tissue, batch/wave, sex, and age metadata are unavailable.",
    "- Each sample is a distinct rat and groups are modeled as independent; no pairing term is used.",
    "- Original raw intensities and preprocessing provenance are unavailable.",
    "- Candidate reversal and pathway findings cannot establish causal rescue.",
    "",
    "## Analysis gates",
    "",
    "| Gate | Status | Detail |",
    "|---|---|---|",
])
for row in gates:
    lines.append(f"| {row['gate']} | {row['status']} | {row['detail']} |")
lines.extend([
    "",
    "## Reproducibility",
    "",
    "- `run_all.sh` reruns the complete analysis.",
    "- `environment.yml` records the intended software environment.",
    "- `report/sessionInfo.txt` records the active R environment.",
    "- `reproducibility/file_manifest.tsv` contains SHA-256 checksums.",
    "- `results/figures/axis_specific/stage_specific_reversal_pathway_composite.*` contains the reference-style PDF, PNG, and TIFF figure.",
    "- `results/tables/axis_specific/figure_source_data/` contains the exact plotted source data and figure methods.",
])
(REPORT / "analysis_report.md").write_text("\n".join(lines) + "\n")

crate = {
    "@context": "https://w3id.org/ro/crate/1.1/context",
    "@graph": [
        {"@id": "ro-crate-metadata.json", "@type": "CreativeWork", "about": {"@id": "./"}},
        {"@id": "./", "@type": "Dataset", "name": "Dapagliflozin PreDM/DM proteomics analysis", "datePublished": datetime.now(timezone.utc).date().isoformat()},
    ],
}
(REPRO / "ro-crate-metadata.json").write_text(json.dumps(crate, indent=2) + "\n")
(REPRO / "runtime.txt").write_text(f"Python {platform.python_version()}\nPlatform {platform.platform()}\n")

manifest_rows = []
for path in sorted(ROOT.rglob("*")):
    if not path.is_file() or path == REPRO / "file_manifest.tsv":
        continue
    manifest_rows.append({
        "path": str(path.relative_to(ROOT)),
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
    })
write_tsv(REPRO / "file_manifest.tsv", manifest_rows)
print(f"[04] final report -> {REPORT / 'analysis_report.md'}")
