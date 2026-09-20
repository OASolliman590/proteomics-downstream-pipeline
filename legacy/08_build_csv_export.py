#!/usr/bin/env python3
"""Build a non-destructive, categorized CSV/figure handoff ZIP."""

from __future__ import annotations

import csv
import hashlib
import os
import shutil
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path


SOURCE = Path(os.environ.get("DAPA_ANALYSIS_ROOT", Path(__file__).resolve().parents[1])).resolve()
EXPORT_PARENT = Path(os.environ.get("DAPA_EXPORT_PARENT", SOURCE.parent)).resolve()
PACKAGE_NAME = os.environ.get("DAPA_EXPORT_NAME", "dapagliflozin_analysis_csv_package_2026-07-20")
PACKAGE = EXPORT_PARENT / PACKAGE_NAME
ZIP_PATH = EXPORT_PARENT / f"{PACKAGE_NAME}.zip"

if PACKAGE.exists() or ZIP_PATH.exists():
    raise SystemExit(f"Refusing to overwrite existing export: {PACKAGE} or {ZIP_PATH}")

for folder in (
    "00_INDEX",
    "01_INPUT_DATA/01_original_workbook",
    "01_INPUT_DATA/02_processed_rat_level",
    "02_QC/tables",
    "02_QC/figures",
    "03_AXIS_SPECIFIC_PRIMARY/01_differential_abundance",
    "03_AXIS_SPECIFIC_PRIMARY/02_three_group_omnibus",
    "03_AXIS_SPECIFIC_PRIMARY/03_reversal",
    "03_AXIS_SPECIFIC_PRIMARY/04_pathways",
    "03_AXIS_SPECIFIC_PRIMARY/05_figure_source_data",
    "03_AXIS_SPECIFIC_PRIMARY/06_figures",
    "04_JOINT_MODEL_SENSITIVITY/01_differential_abundance",
    "04_JOINT_MODEL_SENSITIVITY/02_reversal",
    "04_JOINT_MODEL_SENSITIVITY/03_pathways",
    "04_JOINT_MODEL_SENSITIVITY/04_robustness",
    "04_JOINT_MODEL_SENSITIVITY/05_figures",
    "05_ALL_FIGURES",
    "06_REPORTS",
    "07_PIPELINE_AND_LOGS/pipeline",
    "07_PIPELINE_AND_LOGS/logs",
    "08_REPRODUCIBILITY/native_objects",
):
    (PACKAGE / folder).mkdir(parents=True, exist_ok=False)

catalog: list[dict[str, str | int]] = []


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def description_for(name: str) -> str:
    rules = (
        ("sample_meta", "Rat-level sample metadata and group assignments"),
        ("abundance", "Processed normalized log2 protein abundance matrix by rat"),
        ("protein_annot", "Protein accession, gene symbol, and description mapping"),
        ("sample_qc", "Per-rat missingness, correlation, PCA-distance, and technical-watch metrics"),
        ("protein_coverage", "Per-protein observed-rat counts by experimental group"),
        ("dea_summary", "Differential-abundance counts and minimum P/FDR by contrast"),
        ("axis_dea", "Independent disease-axis limma differential-abundance output"),
        ("omnibus", "Three-group omnibus limma F-test output"),
        ("reversal", "Disease-stage reversal index, class, and evidence annotations"),
        ("disease_axis", "Individual-rat disease-axis scores or treatment tests"),
        ("gsea", "Rank-based pathway enrichment results"),
        ("ora", "Candidate-set over-representation results"),
        ("sensitivity", "Method or threshold robustness results"),
        ("leave_one", "Leave-one-rat/sample stability results"),
        ("contrast_matrix", "Statistical contrast matrix"),
        ("design_matrix", "Limma design matrix"),
        ("pca_scores", "Individual-rat PCA coordinates"),
        ("correlations", "Sample-to-sample abundance correlations"),
        ("analysis_gates", "Scientific and reproducibility gate audit"),
        ("figure", "Exact source data or method definition for a generated figure"),
    )
    lower = name.lower()
    for token, description in rules:
        if token in lower:
            return description
    return "Analysis artifact retained for completeness"


def record(path: Path, category: str, description: str | None = None) -> None:
    catalog.append({
        "category": category,
        "relative_path": path.relative_to(PACKAGE).as_posix(),
        "file_name": path.name,
        "extension": path.suffix.lower().lstrip("."),
        "size_bytes": path.stat().st_size,
        "description": description or description_for(path.name),
    })


def copy_file(src: Path, rel_dest: str, category: str, description: str | None = None) -> Path:
    dest = PACKAGE / rel_dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest)
    record(dest, category, description)
    return dest


def tabular_to_csv(src: Path, rel_dest: str, category: str, description: str | None = None) -> Path:
    dest = PACKAGE / rel_dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    delimiter = "\t" if src.suffix.lower() == ".tsv" else ","
    with src.open("r", encoding="utf-8-sig", newline="") as in_handle, dest.open("w", encoding="utf-8", newline="") as out_handle:
        reader = csv.reader(in_handle, delimiter=delimiter)
        writer = csv.writer(out_handle, lineterminator="\n")
        for row in reader:
            writer.writerow(row)
    record(dest, category, description)
    return dest


def csv_name(path: Path) -> str:
    return f"{path.stem}.csv"


# Original workbook plus processed rat-level inputs.
for src in sorted((SOURCE / "data" / "raw").glob("*")):
    if src.is_file():
        copy_file(src, f"01_INPUT_DATA/01_original_workbook/{src.name}", "Input/original")

for src in sorted((SOURCE / "data" / "processed").glob("*")):
    if not src.is_file():
        continue
    if src.suffix.lower() in {".tsv", ".csv"}:
        tabular_to_csv(src, f"01_INPUT_DATA/02_processed_rat_level/{csv_name(src)}", "Input/processed")
    else:
        copy_file(src, f"08_REPRODUCIBILITY/native_objects/processed_{src.name}", "Reproducibility/native")


def axis_destination(src: Path) -> str:
    name = src.name.lower()
    if "figure_source_data" in src.parts:
        return f"03_AXIS_SPECIFIC_PRIMARY/05_figure_source_data/{csv_name(src)}"
    if "omnibus" in name:
        section = "02_three_group_omnibus"
    elif "reversal" in name or "disease_axis" in name:
        section = "03_reversal"
    elif "gsea" in name or "ora" in name:
        section = "04_pathways"
    else:
        section = "01_differential_abundance"
    return f"03_AXIS_SPECIFIC_PRIMARY/{section}/{csv_name(src)}"


axis_root = SOURCE / "results" / "tables" / "axis_specific"
for src in sorted(axis_root.rglob("*")):
    if not src.is_file():
        continue
    if src.suffix.lower() in {".tsv", ".csv"}:
        tabular_to_csv(src, axis_destination(src), "Primary axis-specific")
    elif src.suffix.lower() in {".md", ".txt"}:
        copy_file(src, f"03_AXIS_SPECIFIC_PRIMARY/05_figure_source_data/{src.name}", "Primary axis-specific/methods")
    else:
        copy_file(src, f"08_REPRODUCIBILITY/native_objects/axis_{src.name}", "Reproducibility/native")


qc_names = {"sample_qc.tsv", "protein_coverage.tsv", "pca_scores.tsv", "sample_correlations.tsv", "analysis_gates.tsv"}
robust_tokens = ("sensitivity", "leave_one", "candidate_leave", "analysis_gates")
reversal_tokens = ("reversal", "disease_axis")
pathway_tokens = ("gsea", "ora", "enrichment")

tables_root = SOURCE / "results" / "tables"
for src in sorted(tables_root.rglob("*")):
    if not src.is_file() or axis_root in src.parents:
        continue
    if src.name in qc_names:
        tabular_to_csv(src, f"02_QC/tables/{csv_name(src)}", "QC")
    elif src.suffix.lower() in {".tsv", ".csv"}:
        lower = src.name.lower()
        if any(token in lower for token in pathway_tokens):
            section = "03_pathways"
        elif any(token in lower for token in reversal_tokens):
            section = "02_reversal"
        elif any(token in lower for token in robust_tokens):
            section = "04_robustness"
        else:
            section = "01_differential_abundance"
        tabular_to_csv(src, f"04_JOINT_MODEL_SENSITIVITY/{section}/{csv_name(src)}", "Joint-model sensitivity")
    else:
        copy_file(src, f"08_REPRODUCIBILITY/native_objects/{src.name}", "Reproducibility/native")


# Figures: preserve original hierarchy twice only where it improves usability:
# a complete figure gallery plus focused axis/QC entry points.
fig_root = SOURCE / "results" / "figures"
for src in sorted(fig_root.rglob("*")):
    if not src.is_file():
        continue
    rel = src.relative_to(fig_root).as_posix()
    copy_file(src, f"05_ALL_FIGURES/{rel}", "Figures/all")
    if "axis_specific" in src.parts:
        copy_file(src, f"03_AXIS_SPECIFIC_PRIMARY/06_figures/{src.name}", "Figures/axis-specific")
    elif "qc" in src.parts:
        copy_file(src, f"02_QC/figures/{src.name}", "Figures/QC")
    else:
        copy_file(src, f"04_JOINT_MODEL_SENSITIVITY/05_figures/{rel}", "Figures/joint sensitivity")


for src in sorted((SOURCE / "report").rglob("*")):
    if src.is_file():
        copy_file(src, f"06_REPORTS/{src.relative_to(SOURCE / 'report').as_posix()}", "Reports")

copy_file(SOURCE / "README.md", "06_REPORTS/PIPELINE_README.md", "Reports")
copy_file(SOURCE / "environment.yml", "08_REPRODUCIBILITY/environment.yml", "Reproducibility/environment")
copy_file(SOURCE / "run_all.sh", "07_PIPELINE_AND_LOGS/run_all.sh", "Pipeline")
for src in sorted((SOURCE / "pipeline").glob("*")):
    if src.is_file() and "__pycache__" not in src.parts:
        copy_file(src, f"07_PIPELINE_AND_LOGS/pipeline/{src.name}", "Pipeline")
for src in sorted((SOURCE / "logs").glob("*")):
    if src.is_file():
        copy_file(src, f"07_PIPELINE_AND_LOGS/logs/{src.name}", "Logs")
for src in sorted((SOURCE / "reproducibility").glob("*")):
    if src.is_file():
        copy_file(src, f"08_REPRODUCIBILITY/{src.name}", "Reproducibility")


# Compact key-findings table for fast orientation.
key_findings = [
    ["section", "axis_or_contrast", "metric", "value", "interpretation"],
    ["Design", "All", "Biological replicates", "4 independent rats/group", "Confirmed by data provider"],
    ["Primary limma", "PreDM+Dapa vs PreDM", "Formal FDR proteins", "0", "No individual protein passed FDR < 0.05 after 3-of-4 coverage gating"],
    ["Primary limma", "DM+Dapa vs DM", "Formal FDR proteins", "0", "No individual protein passed FDR < 0.05 after 3-of-4 coverage gating"],
    ["Reversal", "PreDM", "Full + partial candidates", "20/26", "Exploratory nominal disease set; median RI 0.684"],
    ["Reversal", "DM", "Full + partial candidates", "50/80", "Exploratory nominal disease set; median RI 0.428"],
    ["Pathway", "PreDM+Dapa vs PreDM", "FDR-supported pathways", "3", "EMT down; pentose phosphate and chemokine signaling up"],
    ["Pathway", "DM disease vs CTRL", "Reactome pathways FDR < 0.05", "15", "Predominantly reduced NMDA receptor, ion-channel, and neuronal signaling"],
    ["Pathway", "DM+Dapa vs DM", "Best treatment pathway FDR", "0.0512", "Cell adhesion molecules; supportive, narrowly above 0.05"],
]
key_path = PACKAGE / "00_INDEX" / "KEY_FINDINGS.csv"
with key_path.open("w", newline="", encoding="utf-8") as handle:
    csv.writer(handle, lineterminator="\n").writerows(key_findings)
record(key_path, "Index", "Concise orientation to the principal statistical findings")


# Validate every CSV for parseability and rectangular row structure.
validation_rows = []
for path in sorted(PACKAGE.rglob("*.csv")):
    if path.name in {"CSV_VALIDATION.csv", "CONTENTS.csv", "MANIFEST_SHA256.csv"}:
        continue
    status, detail, n_rows, n_cols = "PASS", "", 0, 0
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            rows = csv.reader(handle)
            header = next(rows, [])
            n_cols = len(header)
            bad_rows = []
            for line_no, row in enumerate(rows, start=2):
                n_rows += 1
                if len(row) != n_cols:
                    bad_rows.append(line_no)
            if not header:
                status, detail = "FAIL", "Missing header"
            elif bad_rows:
                status, detail = "FAIL", f"Non-rectangular rows: {bad_rows[:10]}"
            else:
                detail = "Parsed successfully with consistent column count"
    except Exception as exc:  # pragma: no cover - packaging guard
        status, detail = "FAIL", str(exc)
    validation_rows.append({
        "relative_path": path.relative_to(PACKAGE).as_posix(),
        "status": status,
        "data_rows": n_rows,
        "columns": n_cols,
        "detail": detail,
    })

validation_path = PACKAGE / "00_INDEX" / "CSV_VALIDATION.csv"
with validation_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(validation_rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(validation_rows)
record(validation_path, "Index", "Parse and rectangularity validation for every delivered CSV")

contents_path = PACKAGE / "00_INDEX" / "CONTENTS.csv"
with contents_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(catalog[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(sorted(catalog, key=lambda r: str(r["relative_path"])))
record(contents_path, "Index", "Searchable package inventory with category and description")

readme = f"""# Dapagliflozin PreDM/DM proteomics analysis export

Generated: {datetime.now(timezone.utc).isoformat()}

This is a non-destructive, copy-only handoff of the completed dapagliflozin proteomics analysis. Every numbered abundance column was confirmed to represent one independent rat (n=4/group). The supplied values are processed normalized log2 abundances plus precomputed limma outputs, not raw mass-spectrometry intensities.

## Folder guide

- `00_INDEX`: key findings, searchable contents, CSV validation, and checksum manifest.
- `01_INPUT_DATA`: original workbook and processed rat-level inputs converted to CSV.
- `02_QC`: sample/protein QC tables and figures.
- `03_AXIS_SPECIFIC_PRIMARY`: separate CTRL/PreDM/PreDM+Dapa and CTRL/DM/DM+Dapa limma, omnibus, reversal, pathway, source-data, and figure outputs.
- `04_JOINT_MODEL_SENSITIVITY`: five-group joint-model results, reversal, pathways, and robustness diagnostics.
- `05_ALL_FIGURES`: complete figure gallery preserving the analysis hierarchy.
- `06_REPORTS`: analysis report, methods notes, and software session information.
- `07_PIPELINE_AND_LOGS`: executable scripts and logs.
- `08_REPRODUCIBILITY`: environment, RO-Crate, original analysis manifest, and native R objects.

## Statistical boundaries

- Formal protein discovery requires BH-FDR < 0.05 plus the stated effect-size rule.
- No primary axis-specific individual protein passed FDR < 0.05 after requiring measurements in at least 3/4 rats in every group.
- Reversal-index summaries use exploratory nominal disease sets and must not be called formal protein rescue.
- Pathway GSEA includes FDR-supported results and is labeled separately from protein-level inference.
- Tissue, batch/wave, sex, age, and upstream raw-intensity preprocessing provenance were not supplied.

Start with `00_INDEX/KEY_FINDINGS.csv`, then `06_REPORTS/analysis_report.md`. Figure methods and exact plotted values are under `03_AXIS_SPECIFIC_PRIMARY/05_figure_source_data`.
"""
readme_path = PACKAGE / "README.md"
readme_path.write_text(readme, encoding="utf-8")
record(readme_path, "Documentation", "Package guide, statistical boundaries, and recommended reading order")


# Manifest intentionally excludes itself; all other files, including CONTENTS and README, are hashed.
manifest_rows = []
for path in sorted(PACKAGE.rglob("*")):
    if path.is_file() and path.name != "MANIFEST_SHA256.csv":
        manifest_rows.append({
            "relative_path": path.relative_to(PACKAGE).as_posix(),
            "size_bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
manifest_path = PACKAGE / "00_INDEX" / "MANIFEST_SHA256.csv"
with manifest_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(manifest_rows)

if any(row["status"] != "PASS" for row in validation_rows):
    raise SystemExit("CSV validation failed; package folder retained for inspection and ZIP not created")

with zipfile.ZipFile(ZIP_PATH, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
    for path in sorted(PACKAGE.rglob("*")):
        if path.is_file():
            archive.write(path, arcname=f"{PACKAGE.name}/{path.relative_to(PACKAGE).as_posix()}")

print(f"package={PACKAGE}")
print(f"zip={ZIP_PATH}")
print(f"catalog_files={len(catalog)}")
print(f"csv_files={len(validation_rows)}")
print(f"manifest_files={len(manifest_rows)}")
print(f"zip_bytes={ZIP_PATH.stat().st_size}")
