"""Preprocessing, missingness and QC stage (packet R03: FR-021 to FR-030).

Python builds the typed stage request from the resolved configuration and
the verified R02 canonical bundle; the R handler ``preprocess_stage`` does
all numerical work.  No model is fit here.
"""
from __future__ import annotations

import os
import shutil
from pathlib import Path

from .provenance import sha256_file

CAPABILITY = "preprocessing"


def capabilities():
    return [{"id": CAPABILITY, "dependencies": ["intake"], "required_r_packages": ["jsonlite", "openssl", "proteomicsCore"]}]


def rscript_executable() -> str:
    return os.environ.get("PROTEOMICS_RSCRIPT") or shutil.which("Rscript") or "Rscript"


def coverage_rules(config: dict) -> list[dict]:
    """One rule per declared model x contrast sharing a design (SM04)."""
    rules = []
    for model in config["models"]:
        coverage = model["coverage"]
        for contrast in config["contrasts"]:
            if contrast["design_id"] != model["design_id"]:
                continue
            rules.append({"model_id": model["id"], "contrast_id": contrast["id"], "required_groups": list(contrast["required_groups"]),
                          "policy": coverage["policy"], "minimum_observed_per_group": coverage.get("minimum_observed_per_group", 2),
                          "minimum_fraction": coverage.get("minimum_fraction", 0.5)})
    return rules


def detection_parameters(config: dict) -> dict:
    enabled = bool(config["preprocessing"].get("detection_enabled", False))
    if not enabled:
        return {"enabled": False}
    primary_design = config["design"]["id"]
    families = [f for f in config.get("multiplicity_families", []) if f["hypothesis_type"] == "detection"]
    if families:
        family = families[0]
        selected = [c for c in config["contrasts"] if c["id"] in family["contrast_ids"]]
        family_id = family["id"]
    else:
        selected = [c for c in config["contrasts"] if c["design_id"] == primary_design and len(c["required_groups"]) == 2]
        family_id = "detection-exploratory"
    return {"enabled": True, "family_id": family_id, "comparisons": [{"contrast_id": c["id"], "required_groups": list(c["required_groups"])} for c in selected]}


def build_parameters(config: dict) -> dict:
    pre = config["preprocessing"]
    design = config["design"]
    return {
        "scope": config["runtime"]["scope"],
        "assay": config["assay"],
        "group_column": design["group_column"],
        "group_levels": list(design["group_levels"]),
        "blocking_mode": design["blocking"]["mode"],
        "normalization": pre["normalization"],
        "normalization_sensitivities": list(pre.get("normalization_sensitivities", [])),
        "reference_features": list(pre.get("reference_features", [])),
        "reference_minimum_fraction": pre.get("reference_minimum_fraction", 1.0),
        "tmt_strategy": pre.get("tmt_strategy"),
        "bridge_channel_column": pre.get("bridge_channel_column"),
        "exclusions": [dict(item) for item in pre.get("exclusions", [])],
        "sensitivities": [dict(item) for item in pre.get("sensitivities", [])],
        "coverage_rules": coverage_rules(config) if config["runtime"]["scope"] == "analysis" else [],
        "detection": detection_parameters(config),
        "watchlist_threshold": 3.5,
    }


def build_request(config: dict, *, config_path: str | Path, bundle_dir: str | Path, run_id: str, output_temp_dir: str | Path,
                  stage_id: str = "preprocessing", plan_hash: str | None = None) -> dict:
    manifest = Path(bundle_dir) / "manifest.json"
    return {"schema_version": "1.2.0", "run_id": run_id, "stage_id": stage_id, "capability": CAPABILITY, "plan_hash": plan_hash,
            "inputs": [{"artifact_id": "canonical_manifest", "path": str(manifest.resolve()), "sha256": sha256_file(manifest)}],
            "output_temp_dir": str(Path(output_temp_dir).resolve()), "config_path": str(Path(config_path).resolve()),
            "parameters": build_parameters(config),
            "rng": {"seed": int(config["runtime"]["seed"]), "kind": "L'Ecuyer-CMRG", "threads": int(config["runtime"]["threads"])}}


def execute(request: dict) -> dict:
    """Run the R handler through the R01 stage protocol (no promotion here)."""
    from .runtime import execute_stage
    return execute_stage(request, rscript=rscript_executable())
