"""Deterministic synthetic pathway fixture (R08).

40 proteins map one-to-one to 40 synthetic rat-namespace genes; six gene sets
of 6-10 members; genes in SET_UP are shifted up in group U, SET_DOWN down, the
rest are null.  Seeded standard-library draws; nothing is study data.
"""
from __future__ import annotations

import csv
import json
import random
from pathlib import Path

SETS = {"SET_UP": range(0, 8), "SET_DOWN": range(8, 16), "SET_MIXED": list(range(0, 4)) + list(range(8, 12)), "SET_NULL1": range(20, 28),
        "SET_NULL2": range(28, 36), "SET_TINY": range(36, 38)}


def build(directory: Path, *, design="independent", seed=7, n_features=40):
    rng = random.Random(seed)
    directory.mkdir(parents=True, exist_ok=True)
    if design == "independent":
        obs = [{"observation_id": f"{g}{i}", "group": g} for g in ("C", "U", "T") for i in range(1, 5)]
    else:
        obs = [{"observation_id": f"S{s}_{v}", "group": v, "subject_id": f"S{s}", "biological_unit_id": f"S{s}_{v}"} for s in range(1, 6) for v in ("pre", "post")]
    subject_effect = {f"S{s}": rng.gauss(0, 0.5) for s in range(1, 6)}
    values = {}
    for f in range(n_features):
        row = []
        for o in obs:
            shift = 0.0
            if o["group"] in ("U", "post"):
                shift = 0.9 if f < 8 else (-0.9 if f < 16 else 0.0)
            base = 15 + 0.1 * f + (subject_effect[o["subject_id"]] if design != "independent" else 0.0)
            row.append(base + shift + rng.gauss(0, 0.3))
        values[f"P{f + 1:03d}"] = row
    with (directory / "mapping_source.tsv").open("w", newline="", encoding="utf-8") as h:
        w = csv.writer(h, delimiter="\t", lineterminator="\n"); w.writerow(["source_id", "id_type", "gene_id", "gene_symbol", "taxonomy_id", "status"])
        for f in range(n_features):
            w.writerow([f"P{f + 1:03d}", "synthetic_accession", f"rat:G{f + 1:03d}", f"Gs{f + 1}", 10116, "current"])
    with (directory / "genesets_source.tsv").open("w", newline="", encoding="utf-8") as h:
        w = csv.writer(h, delimiter="\t", lineterminator="\n"); w.writerow(["set_id", "set_name", "gene_id"])
        for s, idx in SETS.items():
            for i in idx:
                w.writerow([s, f"synthetic {s}", f"rat:G{i + 1:03d}"])
    base = {"version": "synthetic-1", "source": "synthetic pathway fixture", "terms": "synthetic", "id_type": "synthetic", "source_taxonomy_id": 10116, "target_taxonomy_id": 10116}
    (directory / "prepare_mapping.json").write_text(json.dumps({**base, "resource_id": "map", "kind": "mapping", "files": [{"name": "mapping.tsv", "source": "mapping_source.tsv"}]}), encoding="utf-8")
    (directory / "prepare_sets.json").write_text(json.dumps({**base, "resource_id": "sets", "kind": "gene_sets", "files": [{"name": "gene_sets.tsv", "source": "genesets_source.tsv"}]}), encoding="utf-8")
    return values, obs
