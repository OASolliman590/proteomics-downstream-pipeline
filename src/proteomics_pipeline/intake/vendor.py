"""Versioned protein-level vendor export mappings (FR-016).

A mapping profile is an explicit, versioned JSON document under
configs/mappings/.  The adapter never infers a vendor or version from
similar-looking headers: the configuration names the profile and version,
the profile declares every column it understands, and any unexpected header
or sub-protein grain is rejected.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

from .wide_long import InputError, RawMatrix, parse_number, read_text_utf8, read_tsv, reject_subprotein_grain

PROFILE_KEYS = {"profile", "version", "vendor", "qualification", "grain", "documented_header_set", "required_columns",
                "feature_id_column", "accessions", "gene_symbols", "gene_ids", "decoy", "contaminant",
                "intensity_columns", "annotation_columns", "ignored_column_patterns", "declared_source_scale",
                "missing_encoding", "notes"}


def load_profile(path: str | Path, requested: str) -> dict:
    """Load a profile and check it is exactly the requested ``profile@version``."""
    path = Path(path)
    try:
        profile = json.loads(read_text_utf8(path))
    except json.JSONDecodeError as exc:
        raise InputError("E_MAPPING_PROFILE", f"mapping profile is not valid JSON: {exc}", file=str(path)) from exc
    if not isinstance(profile, dict) or set(profile) - PROFILE_KEYS or not {"profile", "version", "grain", "required_columns", "feature_id_column", "intensity_columns"} <= set(profile):
        raise InputError("E_MAPPING_PROFILE", "mapping profile has missing or unknown keys", file=str(path))
    identity = f"{profile['profile']}@{profile['version']}"
    if requested != identity:
        raise InputError("E_MAPPING_VERSION", f"configuration requests profile {requested!r} but the mapping file declares {identity!r}; versions are never inferred", file=str(path))
    if profile["grain"] != "protein":
        raise InputError("E_UNSUPPORTED_SCOPE", f"mapping grain {profile['grain']!r} is not protein-level", file=str(path))
    return profile


def check_scale_declaration(profile: dict, *, source_scale: str, zero_as_missing: bool, missing_tokens: list[str]) -> None:
    declared = profile.get("declared_source_scale")
    if declared is not None and declared != source_scale:
        raise InputError("E_MAPPING_SCALE", f"profile {profile['profile']} documents source scale {declared!r}; configuration declares {source_scale!r}")
    encoding = profile.get("missing_encoding") or {}
    if encoding.get("zero") and not zero_as_missing:
        raise InputError("E_MAPPING_MISSING_ENCODING", f"profile {profile['profile']} encodes missing values as 0; set input.linear_zero_as_missing=true explicitly")
    absent = [token for token in encoding.get("tokens", []) if token not in missing_tokens]
    if absent:
        raise InputError("E_MAPPING_MISSING_ENCODING", f"profile {profile['profile']} documents missing tokens {absent}; declare them in input.missing_values")


def _split(value: str, separator: str) -> list[str]:
    return [item.strip() for item in value.split(separator) if item.strip()] if value else []


def _flag(row: dict[str, str], rule: dict | None) -> str:
    if not rule:
        return "unknown"
    value = row[rule["column"]]
    if value in rule.get("true_values", []):
        return "true"
    if value in rule.get("false_values", []):
        return "false"
    return "unknown"


def map_vendor_table(path: str | Path, profile: dict, *, missing_tokens: list[str], source_columns: dict[str, str] | None) -> tuple[RawMatrix, list[dict]]:
    """Map a vendor protein table to canonical matrix values and feature rows.

    ``source_columns`` optionally maps an exact vendor header to an
    observation_id (from the observation metadata ``source_column`` field).
    """
    path = str(path)
    header, body = read_tsv(path)
    reject_subprotein_grain(header, path)
    missing_required = [name for name in profile["required_columns"] if name not in header]
    if missing_required:
        raise InputError("E_MAPPING_HEADER", f"export lacks required {profile['profile']}@{profile['version']} columns {missing_required}; an explicit mapping is required", file=path, row=1)
    known = set(profile["required_columns"]) | set(profile.get("annotation_columns", []))
    ignored = [re.compile(pattern) for pattern in profile.get("ignored_column_patterns", [])]
    spec = profile["intensity_columns"]
    intensity: list[tuple[str, str]] = []
    for name in header:
        if name in known:
            continue
        if spec.get("mode") == "remaining":
            if any(p.search(name) for p in ignored):
                continue
            observation_id = source_columns.get(name) if source_columns else name
            if source_columns and observation_id is None:
                raise InputError("E_MAPPING_HEADER", f"intensity column {name!r} has no observation metadata source_column entry", file=path, row=1, field=name)
            intensity.append((name, observation_id))
            continue
        match = re.fullmatch(spec["regex"], name)
        if match is None and any(p.search(name) for p in ignored):
            continue
        if match is None:
            raise InputError("E_MAPPING_HEADER", f"unknown header {name!r} is not declared by profile {profile['profile']}@{profile['version']}; it needs an explicit mapping, not a guess", file=path, row=1, field=name)
        observation_id = source_columns.get(name) if source_columns else match.group(1)
        if observation_id is None:
            raise InputError("E_MAPPING_HEADER", f"intensity column {name!r} has no observation metadata source_column entry", file=path, row=1, field=name)
        intensity.append((name, observation_id))
    if not intensity:
        raise InputError("E_MAPPING_HEADER", "no intensity columns were recognised by the profile", file=path, row=1)
    observation_ids = [o for _, o in intensity]
    if len(set(observation_ids)) != len(observation_ids):
        raise InputError("E_ID_DUPLICATE", "two intensity columns map to the same observation_id", file=path, row=1)
    matrix = RawMatrix(feature_ids=[], observation_ids=observation_ids)
    features: list[dict] = []
    id_column = profile["feature_id_column"]
    index = {name: position for position, name in enumerate(header)}
    for number, raw in enumerate(body, start=2):
        row = {name: raw[position] for name, position in index.items()}
        feature_id = row[id_column]
        if not feature_id:
            raise InputError("E_ID_MISSING", "feature identifier is empty", file=path, row=number, field=id_column)
        if feature_id in matrix.source_rows:
            raise InputError("E_ID_DUPLICATE", f"duplicate feature identifier {feature_id!r}", file=path, row=number, field=id_column)
        matrix.source_rows[feature_id] = number
        matrix.feature_ids.append(feature_id)
        for name, observation_id in intensity:
            matrix.values[(feature_id, observation_id)] = parse_number(row[name], set(missing_tokens), file=path, row=number, field=name)
        accessions = _split(row[profile["accessions"]["column"]], profile["accessions"]["separator"]) if profile.get("accessions") else []
        symbols = _split(row[profile["gene_symbols"]["column"]], profile["gene_symbols"]["separator"]) if profile.get("gene_symbols") else []
        gene_ids = _split(row[profile["gene_ids"]["column"]], profile["gene_ids"]["separator"]) if profile.get("gene_ids") else []
        features.append({
            "feature_id": feature_id,
            "accessions": json.dumps(accessions, ensure_ascii=False),
            "gene_ids": json.dumps(gene_ids, ensure_ascii=False),
            "gene_symbols": json.dumps(symbols, ensure_ascii=False),
            "is_decoy": _flag(row, profile.get("decoy")),
            "is_contaminant": _flag(row, profile.get("contaminant")),
            "protein_group_ambiguous": "unknown" if not accessions else ("true" if len(accessions) > 1 else "false"),
            "source_profile": f"{profile['profile']}@{profile['version']}",
        })
    if not matrix.feature_ids:
        raise InputError("E_INPUT_EMPTY", "vendor export has no protein rows", file=path)
    return matrix, features
