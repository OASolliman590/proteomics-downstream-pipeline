"""Synthetic studies for the R15d framework (V185-V188). Seeded noise with planted structure; nothing here is study data.

A study is a directory with its analysis configuration (from the Welch fixture factory), its input files and a manifest
whose input hashes are computed from those files. Private and public studies differ only in origin and approval_ref.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


W = _load("welch_factory", ROOT / "tests" / "fixtures" / "methods_extension" / "welch" / "factory.py")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def make_study(directory: Path, *, study_id: str = "synthetic_two_group", origin: str = "private", slot: str = "dda_label_free",
               design: str = "two_group", fmt: str = "canonical_matrix", approval_ref: str | None = None, seed: int = 15) -> Path:
    """Writes a synthetic study and its manifest; returns the manifest path."""
    directory.mkdir(parents=True, exist_ok=True)
    values, observations = W.three_group(shifts={0: {"disease": 2.0}, 1: {"treated": -2.0}}, seed=seed)
    contrasts = [W.contrast("disease-control", "disease", "control")]
    cfg = W.study_config(directory, values, observations, groups=["control", "disease", "treated"], contrasts=contrasts, primary="limma", phase=1)
    files = []
    for name in ("abundance.tsv", "observations.tsv", "features.tsv", "provenance.json"):
        path = directory / name
        files.append({"name": name, "sha256": sha256(path), "bytes": path.stat().st_size})
    manifest = {"schema_version": "1.0.0", "study_id": study_id, "origin": origin, "slot": slot, "design": design,
                "input": {"format": fmt, "files": files}, "organism": "rat", "columns": {"subject": None, "group": "group", "batch": None, "phenotype": None},
                "licence": "SYNTHETIC-CC0", "config": cfg.name}
    if origin == "public":
        manifest["approval_ref"] = approval_ref or "approval-001"
    path = directory / "study-manifest.json"
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return path
