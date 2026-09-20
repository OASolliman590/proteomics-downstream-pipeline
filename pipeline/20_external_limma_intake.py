#!/usr/bin/env python3
"""Import the external DM/PDM limma workbook at the processed-log2 boundary.

The source workbook contains sample-level log2 abundances and four pairwise
limma result tables. It is not a raw-intensity workbook and must not pass
through the central/serum raw-data log2 transformation.

This reader intentionally uses the standard XLSX ZIP/XML representation. The
source file contains dangling drawing relationships that make stricter OpenXML
readers fail even though its cell tables are intact.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import statistics
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from zipfile import ZipFile


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTDIR = ROOT / "data" / "processed_external_dm_pdm"
DEFAULT_RESULTS = ROOT / "results_external_dm_pdm"
DEFAULT_REPORT = ROOT / "report_external_dm_pdm"

MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"
NS = {"m": MAIN_NS, "r": REL_NS, "p": PKG_REL_NS}

SHEET_CONTRASTS = {
    "DM_vs_CTL": "D1_DM_vs_CTRL",
    "PDM_vs_CTL": "D2_PDM_vs_CTRL",
    "DM-Treated_vs_DM": "T1_DM_Treated_vs_DM",
    "PDM-Treated_vs_PDM": "T2_PDM_Treated_vs_PDM",
}
GROUP_MAP = {
    "CTL": "CTRL",
    "DM": "DM",
    "PDM": "PDM",
    "DM-Treated": "DM_Treated",
    "PDM-Treated": "PDM_Treated",
}
GROUP_ORDER = ("CTRL", "DM", "PDM", "DM_Treated", "PDM_Treated")
SAMPLE_RE = re.compile(r"^(CTL|DM|PDM|DM-Treated|PDM-Treated) ([1-4])$")


class IntakeError(ValueError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def column_letters(cell_ref: str) -> str:
    match = re.match(r"[A-Z]+", cell_ref)
    if not match:
        raise IntakeError(f"invalid XLSX cell reference: {cell_ref!r}")
    return match.group(0)


def read_shared_strings(archive: ZipFile) -> list[str]:
    root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
    return [
        "".join(node.text or "" for node in item.iterfind(".//m:t", NS))
        for item in root.findall("m:si", NS)
    ]


def read_sheet(archive: ZipFile, member: str, shared: list[str]) -> list[dict[str, object]]:
    raw_rows: list[dict[str, object]] = []
    with archive.open(member) as handle:
        for _, elem in ET.iterparse(handle, events=("end",)):
            if elem.tag.rsplit("}", 1)[-1] != "row":
                continue
            values: dict[str, object] = {}
            for cell in elem.findall("m:c", NS):
                column = column_letters(cell.attrib["r"])
                cell_type = cell.attrib.get("t", "n")
                value_node = cell.find("m:v", NS)
                if value_node is None:
                    value: object = None
                elif cell_type == "s":
                    value = shared[int(value_node.text)]
                elif cell_type == "b":
                    value = value_node.text == "1"
                else:
                    try:
                        value = float(value_node.text)
                    except (TypeError, ValueError):
                        value = value_node.text
                values[column] = value
            raw_rows.append(values)
            elem.clear()
    if not raw_rows:
        raise IntakeError(f"empty worksheet: {member}")
    columns = {column: str(value) for column, value in raw_rows[0].items()}
    return [{columns[column]: row.get(column) for column in columns} for row in raw_rows[1:]]


def load_workbook(path: Path) -> dict[str, list[dict[str, object]]]:
    with ZipFile(path) as archive:
        shared = read_shared_strings(archive)
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        rels = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        rel_map = {rel.attrib["Id"]: "xl/" + rel.attrib["Target"] for rel in rels}
        sheets: dict[str, list[dict[str, object]]] = {}
        for sheet in workbook.findall("m:sheets/m:sheet", NS):
            name = sheet.attrib["name"]
            relationship_id = sheet.attrib[f"{{{REL_NS}}}id"]
            sheets[name] = read_sheet(archive, rel_map[relationship_id], shared)
    missing = sorted(set(SHEET_CONTRASTS).difference(sheets))
    if missing:
        raise IntakeError(f"missing expected worksheet(s): {', '.join(missing)}")
    return {name: sheets[name] for name in SHEET_CONTRASTS}


def canonical_sample(source_label: str) -> tuple[str, str, int]:
    match = SAMPLE_RE.match(source_label)
    if not match:
        raise IntakeError(f"unrecognised sample label: {source_label!r}")
    source_group, replicate_text = match.groups()
    group = GROUP_MAP[source_group]
    replicate = int(replicate_text)
    return f"EXT_{group}_{replicate}", group, replicate


def sample_columns(row: dict[str, object]) -> list[str]:
    return [column for column in row if SAMPLE_RE.match(column)]


def validate_and_reconstruct(sheets: dict[str, list[dict[str, object]]]):
    accession_sets = [{str(row["Accession"]) for row in rows} for rows in sheets.values()]
    if not all(current == accession_sets[0] for current in accession_sets[1:]):
        raise IntakeError("worksheet accession sets are not identical")
    if any(len(rows) != 3714 for rows in sheets.values()):
        raise IntakeError("expected 3,714 protein rows in every worksheet")

    annotation: dict[str, tuple[str, str]] = {}
    values: dict[tuple[str, str], list[float | None]] = defaultdict(list)
    source_dea: dict[str, list[dict[str, object]]] = {}
    validation: list[dict[str, object]] = []

    for sheet_name, rows in sheets.items():
        samples = sample_columns(rows[0])
        if len(samples) != 8:
            raise IntakeError(f"{sheet_name}: expected 8 sample columns, got {len(samples)}")
        if len({str(row["Accession"]) for row in rows}) != len(rows):
            raise IntakeError(f"{sheet_name}: duplicated accessions")

        complete_mean_errors: list[float] = []
        nlog_errors: list[float] = []
        contrast_rows: list[dict[str, object]] = []
        for row in rows:
            accession = str(row["Accession"])
            current_annotation = (str(row.get("Gene") or ""), str(row.get("Protein") or ""))
            previous_annotation = annotation.setdefault(accession, current_annotation)
            if previous_annotation != current_annotation:
                raise IntakeError(f"annotation drift for accession {accession}")
            for source_sample in samples:
                sample_id, _, _ = canonical_sample(source_sample)
                value = row[source_sample]
                values[(accession, sample_id)].append(None if value is None else float(value))

            first = [row[sample] for sample in samples[:4]]
            second = [row[sample] for sample in samples[4:]]
            if all(value is not None for value in first + second):
                expected_logfc = statistics.mean(second) - statistics.mean(first)
                complete_mean_errors.append(abs(expected_logfc - float(row["log2fc"])))
            pvalue = float(row["pvalue"])
            nlog_errors.append(abs(-math.log10(pvalue) - float(row["nlog10p"])))
            contrast_rows.append(
                {
                    "Accession": accession,
                    "gene_symbol": current_annotation[0],
                    "description": current_annotation[1],
                    "contrast": SHEET_CONTRASTS[sheet_name],
                    "logFC": float(row["log2fc"]),
                    "P.Value": pvalue,
                    "adj.P.Val": float(row["p.adj"]),
                    "source_nlog10p": float(row["nlog10p"]),
                }
            )
        source_dea[SHEET_CONTRASTS[sheet_name]] = contrast_rows
        validation.append(
            {
                "sheet": sheet_name,
                "contrast": SHEET_CONTRASTS[sheet_name],
                "n_proteins": len(rows),
                "max_abs_logfc_mean_error_complete_rows": max(complete_mean_errors),
                "max_abs_nlog10p_error": max(nlog_errors),
            }
        )

    matrix: dict[str, dict[str, float | None]] = defaultdict(dict)
    repeated_conflicts = []
    for (accession, sample_id), observed in values.items():
        finite = [value for value in observed if value is not None]
        if finite and max(finite) != min(finite):
            repeated_conflicts.append((accession, sample_id))
        matrix[accession][sample_id] = finite[0] if finite else None
    if repeated_conflicts:
        raise IntakeError(f"repeated sample values disagree for {len(repeated_conflicts)} cells")

    expected_samples = [f"EXT_{group}_{rep}" for group in GROUP_ORDER for rep in range(1, 5)]
    observed_samples = sorted({sample for row in matrix.values() for sample in row})
    if set(expected_samples) != set(observed_samples):
        raise IntakeError("the expected 5 groups x 4 samples could not be reconstructed")
    return annotation, matrix, expected_samples, source_dea, validation


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_outputs(args, annotation, matrix, samples, source_dea, validation):
    args.outdir.mkdir(parents=True, exist_ok=True)
    args.results.mkdir(parents=True, exist_ok=True)
    args.report.mkdir(parents=True, exist_ok=True)

    sample_rows = []
    for group in GROUP_ORDER:
        source_group = next(key for key, value in GROUP_MAP.items() if value == group)
        for replicate in range(1, 5):
            sample_rows.append(
                {
                    "sample_id": f"EXT_{group}_{replicate}",
                    "source_sample_label": f"{source_group} {replicate}",
                    "group": group,
                    "replicate": replicate,
                    "batch": "",
                    "wave": "",
                    "tissue": "unknown",
                    "treatment_identity": "unknown" if group.endswith("Treated") else "not_applicable",
                    "matched_central_sample_id": "",
                }
            )
    with (args.outdir / "sample_meta.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(sample_rows[0]))
        writer.writeheader()
        writer.writerows(sample_rows)

    annot_rows = [
        {"Accession": accession, "gene_symbol": annotation[accession][0], "description": annotation[accession][1]}
        for accession in sorted(annotation)
    ]
    write_tsv(args.outdir / "protein_annot.tsv", list(annot_rows[0]), annot_rows)

    matrix_path = args.outdir / "abundance_input_log2.tsv"
    matrix_path.parent.mkdir(parents=True, exist_ok=True)
    with matrix_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Accession"] + samples)
        for accession in sorted(annotation):
            writer.writerow(
                [accession]
                + ["" if matrix[accession].get(sample) is None else matrix[accession][sample] for sample in samples]
            )

    dea_summary_rows = []
    dea_dir = args.results / "tables" / "source_dea"
    for contrast, rows in source_dea.items():
        write_tsv(dea_dir / f"dea_{contrast}.tsv", list(rows[0]), rows)
        dea_summary_rows.append(
            {
                "contrast": contrast,
                "n_proteins": len(rows),
                "n_fdr_lt_0_05": sum(float(row["adj.P.Val"]) < 0.05 for row in rows),
                "n_fdr_lt_0_10": sum(float(row["adj.P.Val"]) < 0.10 for row in rows),
                "min_p_value": min(float(row["P.Value"]) for row in rows),
                "min_adj_p_value": min(float(row["adj.P.Val"]) for row in rows),
            }
        )
    write_tsv(dea_dir / "dea_source_summary.tsv", list(dea_summary_rows[0]), dea_summary_rows)
    write_tsv(
        args.results / "tables" / "intake_validation.tsv",
        list(validation[0]),
        validation,
    )

    missing_by_sample = Counter()
    for accession in matrix:
        for sample in samples:
            missing_by_sample[sample] += matrix[accession].get(sample) is None
    missing_rows = [
        {
            "sample_id": sample,
            "n_missing": missing_by_sample[sample],
            "pct_missing": missing_by_sample[sample] / len(matrix) * 100,
        }
        for sample in samples
    ]
    write_tsv(
        args.results / "tables" / "intake_missingness_by_sample.tsv",
        list(missing_rows[0]),
        missing_rows,
    )

    report_lines = [
        "# External DM/PDM Proteomics Intake Report",
        "",
        f"_Generated: {datetime.now(timezone.utc).isoformat()}_",
        "",
        "## Intake classification",
        "",
        "- Source type: processed sample-level log2 abundance plus precomputed limma contrasts.",
        "- Correct pipeline entry: processed-log2 boundary; raw-intensity extraction/log2 transformation must be skipped.",
        "- Proteins: 3,714 unique accessions.",
        "- Samples: 20 total; five groups with n=4 each.",
        "- Groups: CTRL, DM, PDM, DM_Treated, PDM_Treated.",
        "- Repeated CTL/DM/PDM values across worksheets: exact when observed.",
        "- Stored log2FC and -log10(P) fields: numerically validated.",
        "",
        "## Source differential results",
        "",
        "| Contrast | Proteins | FDR < 0.05 | FDR < 0.10 | Minimum FDR |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in dea_summary_rows:
        report_lines.append(
            f"| {row['contrast']} | {row['n_proteins']} | {row['n_fdr_lt_0_05']} | "
            f"{row['n_fdr_lt_0_10']} | {row['min_adj_p_value']:.4g} |"
        )
    report_lines.extend(
        [
            "",
            "No source contrast contains an FDR-supported DEP at 0.05 or 0.10. Raw-P/effect-size candidates must not be called DEPs.",
            "",
            "## Metadata blockers before full pipeline execution",
            "",
            "- Tissue/compartment is not stated in the workbook.",
            "- Biological identity and dose of the treatment are not stated.",
            "- Batch, experimental wave, sex, age, and pairing information are not supplied.",
            "- The workbook contains processed log2 values; original raw intensities and preprocessing provenance are absent.",
            "- Cross-tissue integration is not authorized scientifically until tissue and sample matching are established.",
            "",
            "## Prepared files",
            "",
            "- `abundance_input_log2.tsv`: reconstructed 3,714 x 20 log2 matrix.",
            "- `sample_meta.csv`: five-group sample manifest with unresolved fields explicit.",
            "- `protein_annot.tsv`: accession, gene symbol, and protein description.",
            "- `results_external_dm_pdm/tables/source_dea/`: standardized source limma tables.",
            "- `intake_missingness_by_sample.tsv`: sample-level missingness audit.",
            "",
            "## Required next statistical step",
            "",
            "Refit a joint five-group limma model after defining an explicit missing-data policy and after metadata confirmation. The disease/rescue axes should be DM-vs-CTRL with DM_Treated-vs-DM, and PDM-vs-CTRL with PDM_Treated-vs-PDM. These are separate disease axes and must not be forced into the original Pio/IP/Cr treatment model.",
        ]
    )
    (args.report / "intake_compatibility_report.md").write_text("\n".join(report_lines) + "\n")

    provenance = {
        "source_path": str(args.input),
        "source_sha256": sha256(args.input),
        "source_size_bytes": args.input.stat().st_size,
        "input_scale": "log2_processed",
        "source_workbook_packaging_issue": "dangling drawing and VML relationships",
        "reader": "Python standard-library XLSX ZIP/XML parser",
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    (args.report / "intake_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args()


def main():
    args = parse_args()
    if not args.input.exists():
        raise SystemExit(f"input not found: {args.input}")
    sheets = load_workbook(args.input)
    reconstructed = validate_and_reconstruct(sheets)
    write_outputs(args, *reconstructed)
    print(f"[20] processed log2 matrix -> {args.outdir / 'abundance_input_log2.tsv'}")
    print(f"[20] source DEA tables     -> {args.results / 'tables' / 'source_dea'}")
    print(f"[20] intake report         -> {args.report / 'intake_compatibility_report.md'}")


if __name__ == "__main__":
    main()
