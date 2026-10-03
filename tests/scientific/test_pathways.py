"""R08 acceptance through the real pipeline: V071-V075, V077, V078, V080 (R-level V076/V079 in test-pathways.R)."""
from __future__ import annotations

import copy
import importlib.util
import json
from pathlib import Path

import pytest

from proteomics_pipeline import pathway_service, resources, workflow
from proteomics_pipeline.errors import ProteomicsError

ROOT = Path(__file__).resolve().parents[2]


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path); module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module); return module


B = _load("design_builders", ROOT / "tests" / "fixtures" / "design" / "builders.py")
PF = _load("pathway_fixture", ROOT / "tests" / "fixtures" / "pathways" / "make_pathway_fixture.py")


def _fgsea():
    if B.rscript() is None:
        return False
    try:
        return bool(B.r_json("cat(jsonlite::toJSON(requireNamespace('fgsea', quietly=TRUE), auto_unbox=TRUE))"))
    except RuntimeError:
        return False


pytestmark = pytest.mark.skipif(not _fgsea(), reason="NOT_RUN: Rscript/proteomicsCore/limma/fgsea unavailable")


def fam(fid, htype, contrasts):
    return {"id": fid, "hypothesis_type": htype, "model_ids": ["limma-main"], "contrast_ids": contrasts, "collection_ids": ["sets"], "adjustment": "BH",
            "denominator": "finite_eligible", "q_cutoff": 0.05, "role": "secondary"}


def make(tmp_path, *, design="independent", methods=("camera", "roast", "fgsea", "ora"), nrot=199, profile="smoke_test", mutate=None):
    d = tmp_path / "data"
    values, obs = PF.build(d, design=design)
    files = B.dataset(d, values, obs)
    snaps = {k: resources.prepare(d / f"prepare_{k}.json", d / "snap" / k) for k in ("mapping", "sets")}
    if design == "independent":
        contrasts = [B.contrast("U-C", "U", "C"), B.contrast("T-C", "T", "C", role="secondary")]
        groups, over = ["C", "U", "T"], {}
    else:
        contrasts = [B.contrast("post-pre", "post", "pre")]
        groups = ["pre", "post"]
        over = {"blocking": {"mode": design, "subject_column": "subject_id"}}
    cids = [c["id"] for c in contrasts]
    def m(c):
        c["runtime"]["phase"] = 2; c["runtime"]["execution_profile"] = profile
        c["resources"] = [{"id": k2, "kind": kind, "path": str(Path(snaps[k]["snapshot"]) / "manifest.json"), "sha256": snaps[k]["manifest_sha256"], "version": "synthetic-1",
                           "source": "synthetic pathway fixture", "source_taxonomy_id": 10116, "target_taxonomy_id": 10116, "terms": "synthetic"}
                          for k, k2, kind in (("mapping", "map", "mapping"), ("sets", "sets", "gene_sets"))]
        c["pathways"] = {"enabled": True, "methods": list(methods), "mapping_resource_id": "map", "gene_set_resource_ids": ["sets"], "min_size": 3, "max_size": 20, "representative_rule": "coverage_median_stable_id", "multi_gene_policy": "exclude",
                         "nrot": nrot, "midp": False, "rank": "native_zero_null", "eps": 0, "ora_foreground": "family_q", "ora_threshold": 0.05}
        for method in methods:
            for h in pathway_service.METHOD_TYPES[method]:
                c["multiplicity_families"].append(fam(f"{h}-fam", h, cids))
        if mutate:
            mutate(c)
    path = B.config(d, files, groups=groups, contrasts=contrasts, design_overrides=over, mutate=m)
    payload, code = workflow.run_command(path, tmp_path / "run")
    return payload, code, tmp_path / "run"


def table(run, name):
    return B.read_tsv(run / "pathways" / f"{name}.tsv")


def test_v072_v074_v075_independent_camera_roast_fgsea_match_direct_calls(tmp_path):
    payload, code, run = make(tmp_path)
    assert code == 0, payload
    seed = json.loads((run / "plan.json").read_text(encoding="utf-8"))["runtime"]["seed"]
    oracle = B.r_json("""a <- commandArgs(TRUE); G <- as.matrix(read.delim(a[1], row.names=1, check.names=FALSE)); mem <- read.delim(a[2]); sm <- read.delim(a[3])
      X <- as.matrix(read.delim(a[4], row.names=1, check.names=FALSE)); el <- sm$set_id[sm$eligible == "true"]
      mem <- mem[mem$in_universe == "true" & mem$set_id %in% el,]; sets <- split(mem$gene_id, mem$set_id); idx <- lapply(sets, function(g) match(g, rownames(G)))
      w <- c(-1, 1, 0); cam <- limma::camera(G, idx, design=X, contrast=w, inter.gene.cor=NA, sort=FALSE)
      set.seed(as.integer(a[5]), kind="L'Ecuyer-CMRG"); ro <- limma::mroast(G, idx, design=X, contrast=w, set.statistic="mean", nrot=199, approx.zscore=TRUE, midp=FALSE, sort="none")
      cat(jsonlite::toJSON(list(sets=rownames(cam), cam_p=cam$PValue, cam_cor=cam$Correlation, ro_sets=rownames(ro), ro_p=ro$PValue, ro_mixed=ro$PValue.Mixed), digits=NA))""",
      str(run / "resources" / "gene_matrix_finite.tsv"), str(run / "resources" / "gene_set_membership.tsv"), str(run / "resources" / "gene_set_summary.tsv"),
      str(run / "designs" / "designs" / "joint" / "design_matrix.tsv"), str(seed))
    camera = {r["set_id"]: r for r in table(run, "competitive_enrichment") if r["contrast_id"] == "U-C"}
    for s, p, cor in zip(oracle["sets"], oracle["cam_p"], oracle["cam_cor"]):
        assert abs(float(camera[s]["p_value"]) - p) <= 1e-8 + 1e-6 * p and abs(float(camera[s]["correlation"]) - cor) <= 1e-10
    assert float(camera["SET_UP"]["p_value"]) < 0.05 and camera["SET_UP"]["direction"] == "Up" and "SET_TINY" not in camera   # tiny set ineligible before testing
    directional = {r["set_id"]: r for r in table(run, "self_contained_directional") if r["contrast_id"] == "U-C"}
    mixed = {r["set_id"]: r for r in table(run, "self_contained_mixed") if r["contrast_id"] == "U-C"}
    for s, p, pm in zip(oracle["ro_sets"], oracle["ro_p"], oracle["ro_mixed"]):
        assert float(directional[s]["p_value"]) == pytest.approx(p) and float(mixed[s]["p_value"]) == pytest.approx(pm)
        assert directional[s]["resolution_label"] == "smoke_test_only"
    # V074: central q recomputed per family from the correct raw P
    for name in ("self_contained_directional", "self_contained_mixed"):
        rows = table(run, name); ps = [float(r["p_value"]) for r in rows]; n = len(ps)
        order = sorted(range(n), key=lambda i: ps[i]); bh = [0.0] * n; running = 1.0
        for rank in range(n - 1, -1, -1):
            i = order[rank]; running = min(running, ps[i] * n / (rank + 1)); bh[i] = min(running, 1.0)
        assert all(abs(float(r["q_value"]) - q) <= 1e-12 for r, q in zip(rows, bh)) and {r["family_id"] for r in rows} == {f"{name}-fam"}
    fg = table(run, "preranked_gene_set")
    assert {r["null_type"] for r in fg} == {"preranked_gene_set_exploratory"} and all(r["exploratory"] == "true" for r in fg)
    ranks = B.r_json("""a <- commandArgs(TRUE); dea <- read.delim(a[1]); gm <- read.delim(a[2]); mem <- read.delim(a[3]); sm <- read.delim(a[4])
      reps <- gm[gm$representative_state == "representative",]; d <- dea[dea$contrast_id == "U-C" & dea$eligibility == "tested",]
      st <- setNames(d$statistic[match(reps$feature_id, d$feature_id)], reps$gene_id); st <- st[is.finite(st)]; st <- st[order(-st, names(st), method="radix")]
      el <- sm$set_id[sm$eligible == "true"]; m <- mem[mem$set_id %in% el & !is.na(mem$gene_id),]; pw <- lapply(split(m$gene_id, m$set_id), unique)
      set.seed(as.integer(a[5]), kind="L'Ecuyer-CMRG"); r <- fgsea::fgsea(pw, st, minSize=1, maxSize=length(st)-1, eps=0, nproc=1)
      cat(jsonlite::toJSON(list(sets=r$pathway, p=r$pval, nes=r$NES), digits=NA))""", str(run / "dea" / "zero_null.tsv"), str(run / "resources" / "gene_mapping.tsv"),
      str(run / "resources" / "gene_set_membership.tsv"), str(run / "resources" / "gene_set_summary.tsv"), str(seed))
    fgu = {r["set_id"]: r for r in fg if r["contrast_id"] == "U-C"}
    for s, p, nes in zip(ranks["sets"], ranks["p"], ranks["nes"]):
        assert abs(float(fgu[s]["p_value"]) - p) <= 1e-10 and abs(float(fgu[s]["nes"]) - nes) <= 1e-10
    ora = table(run, "ora_up")
    up = next(r for r in ora if r["contrast_id"] == "U-C" and r["set_id"] == "SET_UP")
    assert int(up["k"]) >= 1 and up["status"] == "tested" and all(r["N"] == ora[0]["N"] for r in ora if r["contrast_id"] == "U-C")
    data = json.loads((run / "report" / "report_data.json").read_text(encoding="utf-8"))
    assert data["sections"]["pathways"]["state"] in ("COMPLETED", "NOT_RUN")


def test_v071_v073_blocked_designs_dispatch_to_roast(tmp_path):
    payload, code, run = make(tmp_path / "camera", design="fixed_subject", methods=("camera",))
    assert code == 2 and payload["error"]["code"] == "E_CAMERA_BLOCKED_DESIGN"
    payload, code, run = make(tmp_path / "dc", design="duplicate_correlation", methods=("roast",))
    assert code == 0, payload
    plan = json.loads((run / "plan.json").read_text(encoding="utf-8")); rho = plan["blocking"][0]["consensus_correlation"]; seed = plan["runtime"]["seed"]
    oracle = B.r_json("""a <- commandArgs(TRUE); G <- as.matrix(read.delim(a[1], row.names=1, check.names=FALSE)); mem <- read.delim(a[2]); sm <- read.delim(a[3])
      X <- as.matrix(read.delim(a[4], row.names=1, check.names=FALSE)); el <- sm$set_id[sm$eligible == "true"]; mem <- mem[mem$in_universe == "true" & mem$set_id %in% el,]
      idx <- lapply(split(mem$gene_id, mem$set_id), function(g) match(g, rownames(G))); block <- sub("_.*$", "", colnames(G))
      set.seed(as.integer(a[6]), kind="L'Ecuyer-CMRG"); ro <- limma::mroast(G, idx, design=X, contrast=c(-1,1), block=block, correlation=as.numeric(a[5]), set.statistic="mean", nrot=199, approx.zscore=TRUE, midp=FALSE, sort="none")
      cat(jsonlite::toJSON(list(sets=rownames(ro), p=ro$PValue), digits=NA))""", str(run / "resources" / "gene_matrix_finite.tsv"), str(run / "resources" / "gene_set_membership.tsv"),
      str(run / "resources" / "gene_set_summary.tsv"), str(run / "designs" / "designs" / "joint" / "design_matrix.tsv"), repr(rho), str(seed))
    rows = {r["set_id"]: r for r in table(run, "self_contained_directional")}
    for s, p in zip(oracle["sets"], oracle["p"]):
        assert float(rows[s]["p_value"]) == pytest.approx(p)
    payload, code, run = make(tmp_path / "prod", methods=("roast",), profile="production")
    assert code == 2 and payload["error"]["code"] == "E_ROTATION_RESOLUTION"


def test_v078_alternative_engine_primary_needs_named_linear_sensitivity():
    cfg = json.loads((B.EXAMPLES / "example-independent.json").read_text(encoding="utf-8"))
    cfg["runtime"]["phase"] = 2; cfg["primary_engine"] = "proda"
    cfg["models"][0].update({"engine": "proda", "coverage": {"policy": "native_dropout", "mask_source": "original_observed"}})
    cfg["pathways"] = {"enabled": True, "methods": ["camera"], "mapping_resource_id": "map", "gene_set_resource_ids": ["sets"]}
    cfg["multiplicity_families"].append(fam("cam", "competitive_enrichment", ["disease-control"]))
    with pytest.raises(ProteomicsError) as error:
        pathway_service.plan_checks(cfg)
    assert error.value.code == "E_PATHWAY_ENGINE_MISMATCH"
    lin = dict(cfg["models"][0]); lin.update({"id": "linear-sens", "engine": "limma", "role": "sensitivity", "execution_requirement": "optional", "trend": True, "robust": True,
                                              "coverage": {"policy": "available_case", "minimum_observed_per_group": 2, "minimum_fraction": 0.5, "mask_source": "original_observed"}})
    cfg["models"].append(lin); cfg["pathways"]["linear_sensitivity_model_id"] = "linear-sens"
    pathway_service.plan_checks(cfg)
    model, role = pathway_service.gene_model(cfg)
    assert model["id"] == "linear-sens" and role == "linear_sensitivity"


def test_v077_v080_families_and_empty_foreground(tmp_path):
    payload, code, run = make(tmp_path, methods=("camera", "ora"), mutate=lambda c: c["pathways"].update({"ora_threshold": 1e-300}))
    assert code == 0, payload
    for name in ("ora_up", "ora_down"):
        rows = table(run, name)
        assert rows and all(r["status"] == "empty_foreground" and float(r["p_value"]) == 1.0 for r in rows)          # all eligible sets kept, p = 1
    cam = table(run, "competitive_enrichment")
    assert {r["family_id"] for r in cam} == {"competitive_enrichment-fam"} and not any(r["hypothesis_type"] != "competitive_enrichment" for r in cam)
    assert len({(r["contrast_id"], r["set_id"]) for r in cam}) == len(cam)
