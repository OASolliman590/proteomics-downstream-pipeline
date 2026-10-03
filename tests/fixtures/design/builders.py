"""Builders for small synthetic design/limma datasets (R04/R05 tests).

Every dataset is written explicitly from the arguments; nothing is drawn
from study data.  Tests author their own oracle values.
"""
from __future__ import annotations

import copy
import csv
import json
import os
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
EXAMPLES = ROOT / "configs" / "examples"


def write_tsv(path: Path, rows) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for row in rows:
            writer.writerow(["NA" if v is None else v for v in row])
    return path


def read_tsv(path: Path) -> list[dict]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def dataset(directory: Path, values: dict[str, list[float | None]], observations: list[dict], *, extra_columns=()) -> dict:
    """values: feature_id -> list aligned with observations; observations need observation_id and group."""
    directory.mkdir(parents=True, exist_ok=True)
    ids = [o["observation_id"] for o in observations]
    write_tsv(directory / "abundance.tsv", [["feature_id"] + ids] + [[f] + [None if v is None else repr(float(v)) for v in vals] for f, vals in values.items()])
    columns = ["observation_id", "biological_unit_id", "subject_id", "technical_replicate_id", "group", *extra_columns]
    rows = []
    for o in observations:
        rows.append([o.get("observation_id"), o.get("biological_unit_id", o["observation_id"]), o.get("subject_id", "NA"), o.get("technical_replicate_id", "NA"), o["group"],
                     *[o.get(c, "NA") for c in extra_columns]])
    write_tsv(directory / "observations.tsv", [columns] + rows)
    write_tsv(directory / "features.tsv", [["feature_id", "accessions", "gene_ids", "gene_symbols", "is_decoy", "is_contaminant", "protein_group_ambiguous"]] +
              [[f, json.dumps([f]), "[]", "[]", "unknown", "unknown", "false"] for f in values])
    (directory / "provenance.json").write_text(json.dumps({"kind": "synthetic_test_fixture", "biological_data": False}), encoding="utf-8")
    return {"matrix": "abundance.tsv", "observations": "observations.tsv", "features": "features.tsv", "source_provenance": "provenance.json"}


def config(directory: Path, files: dict, *, groups, contrasts, design_overrides=None, model_overrides=None, families=None, mutate=None) -> Path:
    base = json.loads((EXAMPLES / "example-independent.json").read_text(encoding="utf-8"))
    base["input"].update(files)
    base["design"].update({"group_levels": list(groups)})
    base["design"].update(design_overrides or {})
    base["contrasts"] = contrasts
    base["models"][0].update(model_overrides or {})
    primary = [c["id"] for c in contrasts if c["role"] == "primary"]
    secondary = [c["id"] for c in contrasts if c["role"] == "secondary"]
    fams = families if families is not None else ([{"id": "protein-primary", "hypothesis_type": "protein_zero_null", "model_ids": ["limma-main"], "contrast_ids": primary,
                                                   "adjustment": "BH", "denominator": "finite_eligible", "q_cutoff": 0.05, "role": "primary"}] +
                                                 ([{"id": "protein-secondary", "hypothesis_type": "protein_zero_null", "model_ids": ["limma-main"], "contrast_ids": secondary,
                                                    "adjustment": "BH", "denominator": "finite_eligible", "q_cutoff": 0.05, "role": "secondary"}] if secondary else []))
    base["multiplicity_families"] = fams
    if mutate:
        mutate(base)
    path = directory / "analysis.json"
    path.write_text(json.dumps(base, indent=2), encoding="utf-8")
    return path


def contrast(cid, numerator, denominator, role="primary", weights=None, design_id="joint"):
    return {"id": cid, "design_id": design_id, "label": f"{numerator} versus {denominator}", "estimand": "Mean log2 abundance difference",
            "weights": weights or {f"group.{numerator}": 1, f"group.{denominator}": -1}, "role": role, "required_groups": [numerator, denominator]}


def rscript() -> str | None:
    executable = os.environ.get("PROTEOMICS_RSCRIPT") or shutil.which("Rscript")
    if not executable:
        return None
    probe = subprocess.run([executable, "--vanilla", "-e", "quit(status=if(exists('limma_stage', envir=asNamespace('proteomicsCore')) && requireNamespace('limma', quietly=TRUE)) 0L else 1L)"], capture_output=True)
    return executable if probe.returncode == 0 else None


def r_json(code: str, *args: str):
    """Run an independent R oracle snippet that prints JSON."""
    # CI run 36982402784 (Windows): `Rscript -e` with multi-line code lost its arguments ("invalid 'description' argument");
    # the snippet now runs from a UTF-8 script file and every path is passed with forward slashes.
    from proteomics_pipeline.runtime import run_r_code
    result = run_r_code(code, args, rscript=rscript())
    if result.returncode != 0:
        raise RuntimeError(result.stderr)
    lines = [line for line in result.stdout.splitlines() if line.strip()]
    return json.loads(lines[-1])   # packages may print progress lines before the JSON oracle output
