#!/usr/bin/env python3
"""Performance benchmark (packet R11, V108).

Generates the declared workload (default 20000 features x 100 observations,
five groups giving eight contrasts, 2000 gene sets), runs QC + limma + CAMERA
+ reports through the real CLI workflow and records wall time, per-stage
time and peak memory (largest child process and orchestrator).  Targets:
<= 30 minutes and <= 8 GiB on documented four-core hardware.  The workload is
never downsized to pass; ``--scale`` exists only for smoke runs and is
labelled as such.  Usage:
  R_LIBS_USER=$PWD/.r-lib python scripts/maintained/benchmark.py [--scale 1.0] [--output DIR]
"""
from __future__ import annotations

import argparse
import json
import platform
import random
import resource
import subprocess
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))


def hardware() -> dict:
    def sysctl(key):
        try:
            return subprocess.run(["sysctl", "-n", key], capture_output=True, text=True, encoding="utf-8").stdout.strip()
        except OSError:
            return None
    return {"platform": platform.platform(), "machine": platform.machine(), "cpu": sysctl("machdep.cpu.brand_string"), "logical_cpus": sysctl("hw.ncpu"),
            "memory_bytes": sysctl("hw.memsize"), "python": platform.python_version()}


def build(directory: Path, n_features: int, n_obs: int, n_sets: int, seed: int = 20261002) -> Path:
    from proteomics_pipeline import resources
    rng = random.Random(seed)
    groups = ["G0", "G1", "G2", "G3", "G4"]
    per = n_obs // len(groups)
    obs = [f"{g}_{i}" for g in groups for i in range(per)]
    directory.mkdir(parents=True)
    with (directory / "abundance.tsv").open("w", encoding="utf-8") as h:
        h.write("feature_id\t" + "\t".join(obs) + "\n")
        for f in range(n_features):
            shift = [0.0, 0.8 if f % 20 == 0 else 0.0, -0.8 if f % 25 == 0 else 0.0, 0.0, 0.5 if f % 30 == 0 else 0.0]
            vals = []
            for gi, g in enumerate(groups):
                for _ in range(per):
                    v = 18 + (f % 50) * 0.05 + shift[gi] + rng.gauss(0, 0.4)
                    vals.append("NA" if rng.random() < 0.0005 else repr(round(v, 6)))
            h.write(f"P{f:06d}\t" + "\t".join(vals) + "\n")
    with (directory / "observations.tsv").open("w", encoding="utf-8") as h:
        h.write("observation_id\tbiological_unit_id\tsubject_id\ttechnical_replicate_id\tgroup\n")
        for o in obs:
            h.write(f"{o}\t{o}\tNA\tNA\t{o.split('_')[0]}\n")
    with (directory / "features.tsv").open("w", encoding="utf-8") as h:
        h.write("feature_id\taccessions\tgene_ids\tgene_symbols\tis_decoy\tis_contaminant\tprotein_group_ambiguous\n")
        for f in range(n_features):
            h.write(f'P{f:06d}\t"[""P{f:06d}""]"\t[]\t[]\tunknown\tunknown\tfalse\n')
    (directory / "provenance.json").write_text(json.dumps({"kind": "synthetic_benchmark"}), encoding="utf-8")
    with (directory / "mapping_source.tsv").open("w", encoding="utf-8") as h:
        h.write("source_id\tid_type\tgene_id\tgene_symbol\ttaxonomy_id\tstatus\n")
        for f in range(n_features):
            h.write(f"P{f:06d}\tsynthetic\tG{f:06d}\tg{f}\t10116\tcurrent\n")
    with (directory / "sets_source.tsv").open("w", encoding="utf-8") as h:
        h.write("set_id\tset_name\tgene_id\n")
        for s in range(n_sets):
            for g in rng.sample(range(n_features), 20):
                h.write(f"SET{s:05d}\tsynthetic\tG{g:06d}\n")
    base = {"version": "bench-1", "source": "synthetic benchmark", "terms": "synthetic", "id_type": "synthetic", "source_taxonomy_id": 10116, "target_taxonomy_id": 10116}
    (directory / "prep_map.json").write_text(json.dumps({**base, "resource_id": "map", "kind": "mapping", "files": [{"name": "mapping.tsv", "source": "mapping_source.tsv"}]}), encoding="utf-8")
    (directory / "prep_sets.json").write_text(json.dumps({**base, "resource_id": "sets", "kind": "gene_sets", "files": [{"name": "gene_sets.tsv", "source": "sets_source.tsv"}]}), encoding="utf-8")
    snaps = {k: resources.prepare(directory / f"prep_{k}.json", directory / "snap" / k) for k in ("map", "sets")}
    config = json.loads((ROOT / "configs" / "examples" / "example-independent.json").read_text(encoding="utf-8"))
    config["input"].update({"matrix": "abundance.tsv", "observations": "observations.tsv", "features": "features.tsv", "source_provenance": "provenance.json"})
    config["runtime"]["phase"] = 2
    config["design"]["group_levels"] = groups
    pairs = [("G1", "G0"), ("G2", "G0"), ("G3", "G0"), ("G4", "G0"), ("G2", "G1"), ("G3", "G2"), ("G4", "G3"), ("G4", "G1")]
    config["contrasts"] = [{"id": f"{a}-{b}", "design_id": "joint", "label": f"{a} vs {b}", "estimand": "mean log2 difference", "weights": {f"group.{a}": 1, f"group.{b}": -1},
                            "role": "primary" if i < 4 else "secondary", "required_groups": [a, b]} for i, (a, b) in enumerate(pairs)]
    fam = lambda fid, h, ids, role, coll=None: {"id": fid, "hypothesis_type": h, "model_ids": ["limma-main"], "contrast_ids": ids, "adjustment": "BH", "denominator": "finite_eligible",
                                                "q_cutoff": 0.05, "role": role, **({"collection_ids": coll} if coll else {})}
    ids = [c["id"] for c in config["contrasts"]]
    config["multiplicity_families"] = [fam("primary", "protein_zero_null", ids[:4], "primary"), fam("secondary", "protein_zero_null", ids[4:], "secondary"),
                                       fam("camera", "competitive_enrichment", ids, "secondary", ["sets"])]
    config["resources"] = [{"id": rid, "kind": kind, "path": str(Path(snaps[k]["snapshot"]) / "manifest.json"), "sha256": snaps[k]["manifest_sha256"], "version": "bench-1",
                            "source": "synthetic benchmark", "source_taxonomy_id": 10116, "target_taxonomy_id": 10116, "terms": "synthetic"} for k, rid, kind in (("map", "map", "mapping"), ("sets", "sets", "gene_sets"))]
    config["pathways"] = {"enabled": True, "methods": ["camera"], "mapping_resource_id": "map", "gene_set_resource_ids": ["sets"], "representative_rule": "coverage_median_stable_id",
                          "multi_gene_policy": "exclude", "min_size": 10, "max_size": 500}
    path = directory / "analysis.json"; path.write_text(json.dumps(config), encoding="utf-8")
    return path


def main(argv=None):
    parser = argparse.ArgumentParser(); parser.add_argument("--scale", type=float, default=1.0); parser.add_argument("--output", type=Path, default=ROOT / "docs" / "validation" / "benchmarks")
    args = parser.parse_args(argv)
    n_features, n_obs, n_sets = int(20000 * args.scale), 100, int(2000 * args.scale)
    with tempfile.TemporaryDirectory(prefix="proteomics-benchmark-") as scratch:
        scratch = Path(scratch)
        t0 = time.time(); config = build(scratch / "data", n_features, n_obs, n_sets); build_s = time.time() - t0
        start = time.time()
        proc = subprocess.run([sys.executable, "-m", "proteomics_pipeline", "run", "--config", str(config), "--output", str(scratch / "run"), "--json"], capture_output=True, text=True, encoding="utf-8")
        wall = time.time() - start
        children = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        payload = json.loads(proc.stdout) if proc.stdout.strip().startswith("{") else {"raw": proc.stdout[-2000:], "stderr": proc.stderr[-2000:]}
        stages = []
        status_path = scratch / "run" / "run_status.json"
        if status_path.is_file():
            for s in json.loads(status_path.read_text(encoding="utf-8"))["stages"]:
                if s.get("result_path"):
                    r = json.loads((scratch / "run" / s["result_path"]).read_text(encoding="utf-8"))
                    if r.get("started_at") and r.get("finished_at"):
                        stages.append({"stage_id": s["stage_id"], "state": s["state"], "seconds": (datetime.fromisoformat(r["finished_at"].replace("Z", "+00:00")) - datetime.fromisoformat(r["started_at"].replace("Z", "+00:00"))).total_seconds()})
    peak = children if sys.platform == "darwin" else children * 1024
    record = {"workload": {"features": n_features, "observations": n_obs, "contrasts": 8, "gene_sets": n_sets, "scale": args.scale,
                           "label": "release workload" if args.scale == 1.0 else "smoke (downsized; not evidence for the target)"},
              "hardware": hardware(), "exit_code": proc.returncode, "run_state": payload.get("state"), "wall_seconds": round(wall, 1), "data_generation_seconds": round(build_s, 1),
              "peak_child_rss_bytes": peak, "stages": stages, "targets": {"wall_seconds_max": 1800, "memory_bytes_max": 8 * 1024 ** 3},
              "assessment": ("PASS" if proc.returncode == 0 and wall <= 1800 and peak <= 8 * 1024 ** 3 else "FAIL") if args.scale == 1.0 else "NOT_RUN (smoke scale)"}
    args.output.mkdir(parents=True, exist_ok=True)
    name = "benchmark.json" if args.scale == 1.0 else f"benchmark_smoke_{args.scale}.json"
    (args.output / name).write_text(json.dumps(record, indent=2), encoding="utf-8")
    print(json.dumps({k: record[k] for k in ("exit_code", "run_state", "wall_seconds", "peak_child_rss_bytes", "assessment")}))
    return 0 if record["assessment"] in ("PASS", "NOT_RUN (smoke scale)") else 1


if __name__ == "__main__":
    raise SystemExit(main())
