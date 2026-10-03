"""R07 acceptance: V061-V070 (local snapshots, explicit preparation, species-aware mapping, representatives,
ambiguous groups, finite gene matrix and model, set eligibility, ORA universe, regression)."""
from __future__ import annotations

import copy
import importlib.util
import json
import os
import socket
import subprocess
import sys
from pathlib import Path

import pytest

from proteomics_pipeline import resources, workflow
from proteomics_pipeline.errors import ProteomicsError

ROOT = Path(__file__).resolve().parents[2]
FIX = ROOT / "tests" / "fixtures" / "resources"
EXAMPLES = ROOT / "configs" / "examples"
spec = importlib.util.spec_from_file_location("design_builders", ROOT / "tests" / "fixtures" / "design" / "builders.py")
B = importlib.util.module_from_spec(spec); spec.loader.exec_module(B)
needs_r = pytest.mark.skipif(B.rscript() is None, reason="NOT_RUN: Rscript/proteomicsCore/limma unavailable")


def snapshots(tmp_path):
    out = {}
    for name in ("mapping", "genesets", "orthology", "human_genesets"):
        out[name] = resources.prepare(FIX / f"prepare_{name}.json", tmp_path / "snapshots" / name)
    return out


def resource_entry(rid, kind, prepared, source_tax=10116, target_tax=10116):
    return {"id": rid, "kind": kind, "path": str(Path(prepared["snapshot"]) / "manifest.json"), "sha256": prepared["manifest_sha256"], "version": "synthetic-1",
            "source": "synthetic test fixture (not a real database)", "source_taxonomy_id": source_tax, "target_taxonomy_id": target_tax, "terms": "synthetic; no redistribution restrictions"}


def config(tmp_path, snaps, *, gene_sets=("sets",), mutate=None):
    raw = json.loads((EXAMPLES / "example-independent.json").read_text(encoding="utf-8"))
    for key in ("matrix", "observations", "features", "source_provenance"):
        raw["input"][key] = str((EXAMPLES / raw["input"][key]).resolve())
    raw["runtime"]["phase"] = 2
    raw["resources"] = [resource_entry("map", "mapping", snaps["mapping"]), resource_entry("sets", "gene_sets", snaps["genesets"]),
                        resource_entry("orth", "orthology", snaps["orthology"], 9606, 10116), resource_entry("hsets", "gene_sets", snaps["human_genesets"], 9606, 10116)]
    raw["pathways"] = {"enabled": False, "mapping_resource_id": "map", "gene_set_resource_ids": list(gene_sets), "min_size": 2, "max_size": 3}
    if mutate:
        mutate(raw)
    tmp_path.mkdir(parents=True, exist_ok=True)
    path = tmp_path / "analysis.json"; path.write_text(json.dumps(raw), encoding="utf-8")
    return path


# ----------------------------------------------------------------------------- V061
def test_v061_snapshots_verified_before_analysis(tmp_path):
    snaps = snapshots(tmp_path)
    path = config(tmp_path, snaps)
    cfg = workflow.load(path)
    verified = resources.verify_resources(cfg, path.parent)
    assert {v["resource_id"] for v in verified} == {"map", "sets", "orth", "hsets"}
    manifest = json.loads((Path(snaps["mapping"]["snapshot"]) / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["terms"].startswith("synthetic") and manifest["source_taxonomy_id"] == 10116
    import hashlib
    assert manifest["files"][0]["sha256"] == hashlib.sha256((Path(snaps["mapping"]["snapshot"]) / "mapping.tsv").read_bytes()).hexdigest()


@pytest.mark.parametrize("mutation,code", [("byte", "E_RESOURCE_HASH"), ("missing", "E_RESOURCE_MISSING"), ("zero", "E_RESOURCE_HASH"), ("inner", "E_RESOURCE_HASH")])
def test_v061_negative_mutation_missing_and_placeholder(tmp_path, mutation, code):
    snaps = snapshots(tmp_path)
    def mutate(raw):
        if mutation == "zero":
            raw["resources"][0]["sha256"] = "0" * 64
        if mutation == "missing":
            raw["resources"][0]["path"] = str(tmp_path / "nowhere" / "manifest.json")
    path = config(tmp_path, snaps, mutate=mutate)
    if mutation == "byte":
        manifest = Path(snaps["mapping"]["snapshot"]) / "manifest.json"; manifest.write_bytes(manifest.read_bytes() + b" ")
    if mutation == "inner":
        inner = Path(snaps["mapping"]["snapshot"]) / "mapping.tsv"; inner.write_text(inner.read_text(encoding="utf-8").replace("rat:G1", "rat:G9"), encoding="utf-8")
    with pytest.raises(ProteomicsError) as error:
        resources.verify_resources(workflow.load(path), path.parent)
    assert error.value.code == code


# ----------------------------------------------------------------------------- V062
def test_v062_explicit_preparation_then_offline_run(tmp_path, monkeypatch):
    out = subprocess.run([sys.executable, "-m", "proteomics_pipeline", "resources", "prepare", "--manifest", str(FIX / "prepare_mapping.json"), "--output", str(tmp_path / "snap"), "--json"],
                         capture_output=True, text=True, encoding="utf-8")
    assert out.returncode == 0, out.stdout + out.stderr
    payload = json.loads(out.stdout)
    assert payload["resource_id"] == "map" and (tmp_path / "snap" / "manifest.json").is_file()
    remote = json.loads((FIX / "prepare_mapping.json").read_text(encoding="utf-8")); remote["files"][0]["source"] = "https://example.invalid/mapping.tsv"
    (tmp_path / "remote.json").write_text(json.dumps(remote), encoding="utf-8")
    with pytest.raises(ProteomicsError) as error:
        resources.prepare(tmp_path / "remote.json", tmp_path / "snap2")
    assert error.value.code == "E_RESOURCE_SOURCE_UNSUPPORTED"


@needs_r
def test_v062_analysis_without_snapshot_fails_and_never_downloads(tmp_path, monkeypatch):
    snaps = snapshots(tmp_path)
    path = config(tmp_path, snaps, mutate=lambda raw: raw["resources"][0].update({"path": str(tmp_path / "absent" / "manifest.json")}))
    payload, code = workflow.run_command(path, tmp_path / "run")
    assert code == 2 and payload["error"]["code"] == "E_RESOURCE_MISSING"
    def no_network(*_a, **_k):
        raise AssertionError("analysis attempted a network connection")
    monkeypatch.setattr(socket, "create_connection", no_network)
    monkeypatch.setenv("http_proxy", "http://127.0.0.1:9"); monkeypatch.setenv("https_proxy", "http://127.0.0.1:9")
    good = config(tmp_path / "g", snapshots(tmp_path / "g"))
    payload, code = workflow.run_command(good, tmp_path / "g" / "run")
    assert code == 0, payload
    assert any(s["stage_id"] == "resources" and s["state"] == "COMPLETED" for s in payload["stages"])


# ----------------------------------------------------------------------------- V063 / V064
def test_v063_negative_cross_species_join_refused(tmp_path):
    snaps = snapshots(tmp_path)
    path = config(tmp_path, snaps, mutate=lambda raw: raw.update({"organism": {"taxonomy_id": 9606, "scientific_name": "Homo sapiens"}}))
    with pytest.raises(ProteomicsError) as error:
        resources.verify_resources(workflow.load(path), path.parent)
    assert error.value.code == "E_RESOURCE_TAXONOMY"


def test_v064_negative_projection_without_evidence(tmp_path):
    spec_ = json.loads((FIX / "prepare_orthology.json").read_text(encoding="utf-8")); spec_.pop("projection")
    (tmp_path / "p.json").write_text(json.dumps(spec_), encoding="utf-8"); (tmp_path / "orthology_source.tsv").write_text((FIX / "orthology_source.tsv").read_text(encoding="utf-8"), encoding="utf-8")
    with pytest.raises(ProteomicsError) as error:
        resources.prepare(tmp_path / "p.json", tmp_path / "snap")
    assert error.value.code == "E_RESOURCE_TAXONOMY"


@needs_r
def test_v063_to_v069_end_to_end_mapping_and_gene_matrix(tmp_path):
    snaps = snapshots(tmp_path)
    path = config(tmp_path, snaps, gene_sets=("sets", "hsets"))
    payload, code = workflow.run_command(path, tmp_path / "run")
    assert code == 0, payload
    res = tmp_path / "run" / "resources"
    mapping = {r["feature_id"]: r for r in B.read_tsv(res / "gene_mapping.tsv")}
    # V063: only the rat namespace is used; retired and unmapped keep reasons
    assert mapping["P01"]["gene_id"] == "rat:G1" and mapping["P04"]["mapping_reason"] == "retired_identifier" and mapping["P06"]["mapping_reason"] == "unmapped"
    # V066: the two-gene group is excluded; same-gene proteins compete for one representative
    assert mapping["P05"]["mapping_reason"] == "multi_gene_group" and mapping["P05"]["representative_state"] == "excluded"
    assert {mapping["P02"]["representative_state"], mapping["P03"]["representative_state"]} == {"representative", "candidate_not_selected"}
    # V065 oracle: coverage desc, median desc, feature_id asc (P02, P03 both fully observed; higher median wins)
    canonical = {r["feature_id"]: r for r in B.read_tsv(tmp_path / "run" / "inputs" / "matrix.tsv")}
    def median(fid):
        vals = sorted(float(v) for k, v in canonical[fid].items() if k != "feature_id" and v != "NA"); n = len(vals)
        return vals[n // 2] if n % 2 else (vals[n // 2 - 1] + vals[n // 2]) / 2
    winner = max(("P02", "P03"), key=lambda f: (12 - sum(v == "NA" for v in canonical[f].values()), median(f), [-ord(c) for c in f]))
    assert mapping[winner]["representative_state"] == "representative"
    # V067: finite gene matrix and consistent gene model
    loss = {r["gene_id"]: r for r in B.read_tsv(res / "gene_universe_loss.tsv")}
    assert loss["rat:G7"]["in_finite_universe"] == "false" and loss["rat:G7"]["reason"] == "incomplete_in_required_observations"
    finite = B.read_tsv(res / "gene_matrix_finite.tsv")
    assert [r["feature_id"] for r in finite] == ["rat:G1", "rat:G2"]
    oracle = B.r_json("m <- as.matrix(read.delim(commandArgs(TRUE)[1], row.names=1, check.names=FALSE)); g <- factor(sub('[0-9]+$','',colnames(m)), levels=c('C','U','T')); "
                      "X <- model.matrix(~ relevel(g, ref='C')); fit <- limma::eBayes(limma::lmFit(m, X), trend=TRUE, robust=TRUE); "
                      "cat(jsonlite::toJSON(list(t=unname(fit$t[,2]), p=unname(fit$p.value[,2])), digits=NA))", str(res / "gene_matrix_finite.tsv"))
    genes = [r for r in B.read_tsv(res / "gene_zero_null.tsv") if r["contrast_id"] == "disease-control"]
    for row, t, p in zip(genes, oracle["t"], oracle["p"]):
        assert abs(float(row["statistic"]) - t) <= 1e-8 and abs(float(row["p_value"]) - p) <= 1e-8 + 1e-6 * p
    # V068: eligibility on set ∩ universe with bounds 2-3 before testing
    summary = {(r["resource_id"], r["set_id"]): r for r in B.read_tsv(res / "gene_set_summary.tsv")}
    assert summary[("sets", "SET_A")]["eligible_overlap_size"] == "2" and summary[("sets", "SET_A")]["eligible"] == "true"
    assert summary[("sets", "SET_B")]["eligible_overlap_size"] == "0" and summary[("sets", "SET_B")]["eligible"] == "false"
    # V064: projected collection keeps evidence and label; ambiguous orthologs excluded
    member = [r for r in B.read_tsv(res / "gene_set_membership.tsv") if r["resource_id"] == "hsets"]
    assert {r["collection_label"] for r in member} == {"ortholog_projected"}
    assert next(r for r in member if r["source_gene_id"] == "hum:G7")["projection_status"] == "ambiguous_ortholog_excluded"


@needs_r
def test_v070_repeatable_and_label_permutation_invariant(tmp_path):
    snaps = snapshots(tmp_path)
    workflow.run_command(config(tmp_path / "a", snaps), tmp_path / "a" / "run")
    workflow.run_command(config(tmp_path / "b", snaps), tmp_path / "b" / "run")
    for name in ("gene_mapping.tsv", "gene_matrix_finite.tsv", "gene_set_summary.tsv"):
        assert (tmp_path / "a" / "run" / "resources" / name).read_bytes() == (tmp_path / "b" / "run" / "resources" / name).read_bytes()
    # permute group labels (C<->T): representatives are label-independent
    obs = B.read_tsv(EXAMPLES / "fixtures" / "independent-observations.tsv")
    swapped = tmp_path / "swapped.tsv"
    swap = {"control": "treated", "treated": "control"}
    B.write_tsv(swapped, [list(obs[0].keys())] + [[swap.get(o[k], o[k]) if k == "group" else o[k] for k in obs[0]] for o in obs])
    workflow.run_command(config(tmp_path / "c", snaps, mutate=lambda raw: raw["input"].update({"observations": str(swapped)})), tmp_path / "c" / "run")
    def reps(d):
        return {r["feature_id"]: r["representative_state"] for r in B.read_tsv(d / "run" / "resources" / "gene_mapping.tsv")}
    assert reps(tmp_path / "a") == reps(tmp_path / "c")
    corrupt = Path(snaps["genesets"]["snapshot"]) / "gene_sets.tsv"; corrupt.write_text(corrupt.read_text(encoding="utf-8") + "SET_X\tx\trat:G1\n", encoding="utf-8")
    payload, code = workflow.run_command(config(tmp_path / "d", snaps), tmp_path / "d" / "run")
    assert code == 2 and payload["error"]["code"] == "E_RESOURCE_HASH"
