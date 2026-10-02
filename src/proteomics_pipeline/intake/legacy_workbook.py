"""Read-only wrapper around the recovered historical workbook parser (FR-017).

The recovered ``pipeline/20_external_limma_intake.py`` is imported unchanged
and only its sheet reader and sample-label parser are reused.  Its
``validate_and_reconstruct`` is not called because it hardcodes the
historical study dimensions; the reconstruction below repeats its
agreement checks without those constants.  Contrast identity comes from an
explicit Method A/Method B semantic crosswalk keyed by biological
numerator/denominator, never by D1/T1 position.
"""
from __future__ import annotations

import importlib.util
import json
from functools import lru_cache
from pathlib import Path

from ..errors import CapabilityError
from .wide_long import InputError, RawMatrix, read_text_utf8

REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
HISTORICAL_PARSER = REPOSITORY_ROOT / "pipeline" / "20_external_limma_intake.py"

# Method B names are recorded in docs/PIPELINE.md ("Sample and contrast naming").
METHOD_B_GROUPS = {"CTL": "CTRL", "PDM": "PreDM", "DM": "DM", "PDM-Treated": "PreDM_Dapa", "DM-Treated": "DM_Dapa"}


@lru_cache(maxsize=1)
def historical_parser():
    """Import the recovered parser read-only; unavailable source is a capability error."""
    if not HISTORICAL_PARSER.is_file():
        raise CapabilityError("E_CAPABILITY_NOT_AVAILABLE", "the recovered historical workbook parser is not present in this installation")
    spec = importlib.util.spec_from_file_location("historical_external_limma_intake", HISTORICAL_PARSER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def method_groups(method: str) -> dict[str, str]:
    parser = historical_parser()
    if method == "A":
        return dict(parser.GROUP_MAP)
    if method == "B":
        return dict(METHOD_B_GROUPS)
    raise InputError("E_CROSSWALK_METHOD", f"crosswalk method must be 'A' or 'B', got {method!r}")


def load_crosswalk(path: str | Path) -> dict:
    path = Path(path)
    try:
        crosswalk = json.loads(read_text_utf8(path))
    except json.JSONDecodeError as exc:
        raise InputError("E_CROSSWALK", f"crosswalk is not valid JSON: {exc}", file=str(path)) from exc
    if not isinstance(crosswalk, dict) or set(crosswalk) != {"crosswalk_version", "method", "group_names", "contrasts"}:
        raise InputError("E_CROSSWALK", "crosswalk must contain exactly crosswalk_version, method, group_names and contrasts", file=str(path))
    expected = method_groups(crosswalk["method"])
    if crosswalk["group_names"] != expected:
        raise InputError("E_CROSSWALK_METHOD", f"group_names do not match documented Method {crosswalk['method']} naming {expected}; Method A and B names are not interchangeable", file=str(path))
    return crosswalk


def _sheet_labels(sheet: str) -> tuple[str, str]:
    numerator, _, denominator = sheet.partition("_vs_")
    return numerator, denominator


def validate_crosswalk(crosswalk: dict) -> list[dict]:
    """Check every historical sheet maps to the declared biological contrast."""
    parser = historical_parser()
    names = crosswalk["group_names"]
    contrasts = crosswalk["contrasts"]
    sheets = [item.get("sheet") for item in contrasts]
    if sorted(sheets) != sorted(parser.SHEET_CONTRASTS) or len(set(sheets)) != len(sheets):
        raise InputError("E_CROSSWALK", f"crosswalk must map each historical sheet exactly once: {sorted(parser.SHEET_CONTRASTS)}")
    ids = [item.get("id") for item in contrasts]
    if len(set(ids)) != len(ids) or not all(ids):
        raise InputError("E_CROSSWALK", "crosswalk contrast ids must be unique and nonempty")
    resolved = []
    for item in contrasts:
        numerator_label, denominator_label = _sheet_labels(item["sheet"])
        expected = (names[numerator_label], names[denominator_label])
        declared = (item.get("numerator"), item.get("denominator"))
        if declared != expected:
            raise InputError("E_CONTRAST_IDENTITY", f"contrast {item['id']!r} maps sheet {item['sheet']!r} (biologically {expected[0]} vs {expected[1]}) but declares {declared[0]} vs {declared[1]}; update the semantic crosswalk instead of relabeling results")
        resolved.append({"contrast_id": item["id"], "sheet": item["sheet"], "historical_label": parser.SHEET_CONTRASTS[item["sheet"]],
                         "numerator": expected[0], "denominator": expected[1]})
    return resolved


def import_workbook(workbook: str | Path, crosswalk: dict, *, missing_tokens: list[str]) -> tuple[RawMatrix, list[dict], list[dict], list[dict]]:
    """Return (matrix, features, contrast_crosswalk, source_contrast_rows).

    Source differential statistics are retained as historical provenance
    only; they are never treated as maintained inference results.
    """
    parser = historical_parser()
    resolved = validate_crosswalk(crosswalk)
    try:
        sheets = parser.load_workbook(Path(workbook))
    except parser.IntakeError as exc:
        raise InputError("E_LEGACY_WORKBOOK", str(exc), file=str(workbook)) from exc
    except Exception as exc:  # malformed ZIP/XML from the historical reader
        raise InputError("E_LEGACY_WORKBOOK", f"historical reader failed: {exc}", file=str(workbook)) from exc
    accession_sets = [[str(row["Accession"]) for row in rows] for rows in sheets.values()]
    if any(set(a) != set(accession_sets[0]) for a in accession_sets[1:]):
        raise InputError("E_LEGACY_WORKBOOK", "worksheet accession sets are not identical", file=str(workbook))
    feature_ids: list[str] = []
    for accession in accession_sets[0]:
        if accession in feature_ids:
            raise InputError("E_ID_DUPLICATE", f"duplicated accession {accession!r}", file=str(workbook))
        feature_ids.append(accession)
    values: dict[tuple[str, str], float | None] = {}
    observation_ids: list[str] = []
    annotation: dict[str, tuple[str, str]] = {}
    source_rows: list[dict] = []
    by_sheet = {item["sheet"]: item for item in resolved}
    for sheet_name, rows in sheets.items():
        samples = parser.sample_columns(rows[0])
        if not samples:
            raise InputError("E_LEGACY_WORKBOOK", f"sheet {sheet_name!r} has no recognised sample columns", file=str(workbook))
        for row in rows:
            accession = str(row["Accession"])
            current = (str(row.get("Gene") or ""), str(row.get("Protein") or ""))
            if annotation.setdefault(accession, current) != current:
                raise InputError("E_LEGACY_WORKBOOK", f"annotation drift for accession {accession!r}", file=str(workbook))
            for label in samples:
                sample_id, _, _ = parser.canonical_sample(label)
                if sample_id not in observation_ids:
                    observation_ids.append(sample_id)
                raw = row[label]
                value = None if raw is None or (isinstance(raw, str) and raw in missing_tokens) else float(raw)
                key = (accession, sample_id)
                if key in values and values[key] != value:
                    raise InputError("E_LEGACY_CONFLICT", f"repeated sample value disagrees across sheets for {key!r}", file=str(workbook))
                values[key] = value
            item = by_sheet[sheet_name]
            source_rows.append({"contrast_id": item["contrast_id"], "sheet": sheet_name, "historical_label": item["historical_label"],
                                "numerator": item["numerator"], "denominator": item["denominator"], "feature_id": accession,
                                "source_log2fc": row.get("log2fc"), "source_p_value": row.get("pvalue"), "source_adjusted_p": row.get("p.adj"),
                                "status": "historical_source_value_not_refit"})
    matrix = RawMatrix(feature_ids=feature_ids, observation_ids=observation_ids, values=values)
    for feature_id in feature_ids:
        for observation_id in observation_ids:
            matrix.values.setdefault((feature_id, observation_id), None)
    features = [{"feature_id": f, "accessions": json.dumps([f]), "gene_ids": "[]",
                 "gene_symbols": json.dumps([annotation[f][0]] if annotation[f][0] else [], ensure_ascii=False),
                 "is_decoy": "unknown", "is_contaminant": "unknown", "protein_group_ambiguous": "unknown",
                 "source_profile": "historical_workbook"} for f in feature_ids]
    return matrix, features, resolved, source_rows


def check_observation_groups(crosswalk: dict, observations, group_column: str) -> None:
    """Observation groups must be the crosswalk method's names for each sample label."""
    parser = historical_parser()
    names = crosswalk["group_names"]
    reverse_a = {value: key for key, value in parser.GROUP_MAP.items()}
    for row in observations.rows:
        observation_id = row["observation_id"]
        prefix = observation_id.rsplit("_", 1)[0].removeprefix("EXT_")
        label = reverse_a.get(prefix)
        if label is None:
            raise InputError("E_CROSSWALK_GROUP", f"observation {observation_id!r} is not a historical canonical sample ID", field="observation_id")
        if row[group_column] != names[label]:
            raise InputError("E_CROSSWALK_GROUP", f"observation {observation_id!r} has group {row[group_column]!r}; Method {crosswalk['method']} names workbook label {label!r} as {names[label]!r}", field=group_column)
