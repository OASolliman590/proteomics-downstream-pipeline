"""Study manifest and slot input adapters (packet R15d: FR-185, FR-189 to FR-192; ADR 0011).

A study manifest is a versioned JSON document. Every input file is checked against its declared SHA-256 and byte
count before anything is validated or run (E_MANIFEST_HASH). Unknown designs and formats are refused
(E_MANIFEST_DESIGN). A private study's identifier must not carry a path or a person's name. The slot adapters
map the columns of a study's input format through a recorded mapping (its hash is recorded), refuse only the
unmapped column (E_MAPPING_UNMAPPED), apply the MaxQuant contamination and reverse-hit filters, and report
batch confounding before any batch term is fitted.
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

from ..errors import ProteomicsError

MANIFEST_SCHEMA_VERSION = "1.0.0"
SLOTS = ("dda_label_free", "paired_repeated", "tmt_batched", "dia", "multigroup_phenotype_nonhuman", "reference_welch")
DESIGNS = ("two_group", "paired", "repeated", "multi_group")
FORMATS = ("maxquant_proteingroups", "diann", "spectronaut", "canonical_matrix")
REQUIRED_FIELDS = ("study_id", "origin", "slot", "design", "input", "organism", "columns", "licence")
ORIGINS = ("private", "public")
# Columns each input format must map (the names are the format's published headers; the adapter records what it found).
FORMAT_COLUMNS = {
    "maxquant_proteingroups": {"protein_ids": ["Protein IDs"], "reverse": ["Reverse"], "contaminant": ["Potential contaminant"]},
    "diann": {"protein_group": ["Protein.Group"], "genes": ["Genes"]},
    "spectronaut": {"protein_accessions": ["PG.ProteinAccessions"], "genes": ["PG.Genes"]},
    "canonical_matrix": {"feature_id": ["feature_id"]},
}
_NAME_LIKE = re.compile(r"[/\\]|\s")


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_manifest(manifest: dict, manifest_dir: str | Path, *, approvals: dict | None = None) -> dict:
    """Validate the manifest and every input file hash; return the validation record (FR-185, FR-186 step 1)."""
    base = Path(manifest_dir)
    for key in REQUIRED_FIELDS:
        if key not in manifest:
            raise ProteomicsError("E_CONFIG_SCHEMA", f"manifest field {key!r} is required", f"/{key}", exit_code=2)
    if manifest.get("schema_version") != MANIFEST_SCHEMA_VERSION:
        raise ProteomicsError("E_CONFIG_SCHEMA", f"manifest schema_version must be {MANIFEST_SCHEMA_VERSION}", "/schema_version", exit_code=2)
    if manifest["origin"] not in ORIGINS:
        raise ProteomicsError("E_CONFIG_SCHEMA", "origin must be private or public", "/origin", exit_code=2)
    if manifest["slot"] not in SLOTS:
        raise ProteomicsError("E_MANIFEST_DESIGN", f"unknown slot {manifest['slot']!r}", "/slot", exit_code=2)
    if manifest["design"] not in DESIGNS:
        raise ProteomicsError("E_MANIFEST_DESIGN", f"unknown design {manifest['design']!r}", "/design", exit_code=2)
    if manifest["input"].get("format") not in FORMATS:
        raise ProteomicsError("E_MANIFEST_DESIGN", f"unsupported input format {manifest['input'].get('format')!r}", "/input/format", exit_code=2)
    study_id = str(manifest["study_id"])
    if _NAME_LIKE.search(study_id) or not study_id:
        raise ProteomicsError("E_PRIVATE_IDENTIFIER", "study_id must be a local identifier without a path or spaces", "/study_id", exit_code=5)
    if manifest["origin"] == "public":
        approval = (approvals or {}).get(manifest.get("approval_ref"))
        if not manifest.get("approval_ref") or approval is None:
            raise ProteomicsError("E_PUBLIC_DATA_UNAPPROVED", "a public study needs a matching approval record (FR-187)", "/approval_ref", exit_code=2)
    checked = []
    for item in manifest["input"].get("files", []):
        path = (base / item["name"]).resolve()
        if not path.is_file():
            raise ProteomicsError("E_MANIFEST_HASH", f"input file {item['name']!r} is missing", "/input/files", exit_code=5)
        digest = file_sha256(path)
        if digest != item["sha256"] or path.stat().st_size != int(item["bytes"]):
            raise ProteomicsError("E_MANIFEST_HASH", f"input file {item['name']!r} differs from its declared SHA-256 or size; the study is not run", "/input/files", exit_code=5)
        checked.append({"name": item["name"], "sha256": digest, "bytes": path.stat().st_size})
    return {"state": "PASS", "study_id": study_id, "slot": manifest["slot"], "design": manifest["design"], "origin": manifest["origin"],
            "input_format": manifest["input"]["format"], "files_checked": len(checked), "manifest_schema_version": MANIFEST_SCHEMA_VERSION}


def map_columns(input_format: str, header: list[str]) -> dict:
    """Map the required columns of a format through the header; an unmapped column is refused on its own (FR-192)."""
    if input_format not in FORMAT_COLUMNS:
        raise ProteomicsError("E_MANIFEST_DESIGN", f"unsupported input format {input_format!r}", "/input/format", exit_code=2)
    mapping, unmapped = {}, []
    for logical, candidates in FORMAT_COLUMNS[input_format].items():
        found = next((c for c in candidates if c in header), None)
        if found is None:
            unmapped.append(logical)
        else:
            mapping[logical] = found
    if unmapped:
        raise ProteomicsError("E_MAPPING_UNMAPPED", f"required column(s) {unmapped} are not in the {input_format} header; only these are refused", "/input/format", exit_code=2)
    record = {"format": input_format, "mapping": mapping}
    record["mapping_sha256"] = hashlib.sha256(json.dumps(mapping, sort_keys=True).encode("utf-8")).hexdigest()
    return record


def maxquant_filters(rows: list[dict], mapping: dict) -> tuple[list[dict], dict]:
    """The declared MaxQuant filters: reverse hits and potential contaminants are removed and counted (FR-189)."""
    kept, removed_reverse, removed_contaminant = [], 0, 0
    for row in rows:
        if row.get(mapping["reverse"]) == "+":
            removed_reverse += 1
        elif row.get(mapping["contaminant"]) == "+":
            removed_contaminant += 1
        else:
            kept.append(row)
    return kept, {"input_rows": len(rows), "removed_reverse": removed_reverse, "removed_contaminant": removed_contaminant, "kept": len(kept)}


def tmt_batch_confounding(observations: list[dict], group: str, batch: str) -> dict:
    """Group-by-batch cross-tabulation before any batch term is fitted (FR-191). A batch that lacks a group level is
    confounded: it is recorded as an adaptation and its batch term is not fitted."""
    groups = sorted({o[group] for o in observations})
    batches = sorted({o[batch] for o in observations})
    table = {b: {g: sum(1 for o in observations if o[batch] == b and o[group] == g) for g in groups} for b in batches}
    confounded = [b for b in batches if any(table[b][g] == 0 for g in groups)]
    adaptation = None
    if confounded:
        adaptation = {"analysis": "tmt_batch", "item": "batch", "requested": "batch_adjusted", "used": "unadjusted_within_batches", "reason": "batch_confounded_with_group",
                      "batches": confounded}
    return {"groups": groups, "batches": batches, "table": table, "confounded_batches": confounded, "adaptation": adaptation,
            "batch_term_fitted": not confounded}


def write_json(path: Path, value: dict) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8", newline="\n")
    return path
