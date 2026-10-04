"""R14e acceptance through the real Python -> R pipeline: V161-V165 (SM39, SM40 with SM14).

Oracles are written here: the declared minimum of units, the planted module labels with an adjusted Rand index
computed in the test, prcomp() eigengenes, a direct limma fit on eigengenes, the snapshot manifest hash, and a
degree-preserving sampling loop over the measured mapped universe with the recorded seed.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

from proteomics_pipeline import resources, workflow
from proteomics_pipeline.errors import ProteomicsError
from proteomics_pipeline.reproduction import offline

ROOT = Path(__file__).resolve().parents[2]


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path); module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module); return module


F = _load("post_de_factory", ROOT / "tests" / "fixtures" / "post_de" / "design_factory.py")
B = F.B


def _ready():
    if B.rscript() is None:
        return False
    try:
        return bool(B.r_json("cat(jsonlite::toJSON(requireNamespace('dynamicTreeCut', quietly=TRUE), auto_unbox=TRUE))"))
    except RuntimeError:
        return False


pytestmark = pytest.mark.skipif(not _ready(), reason="NOT_RUN: Rscript/proteomicsCore/dynamicTreeCut unavailable")


def coab_run(tmp_path, *, n_per_group=30, coabundance=True, extra=None, expect=0):
    values, obs = F.coabundance_design(n_per_group=n_per_group)
    files = F.write_dataset(tmp_path / "data", values, obs)
    networks = {"enabled": True}
    if coabundance:
        networks["coabundance"] = {"enabled": True, "min_units": 20, "rule": "wgcna_signed", "min_module_size": 8, "bootstrap": 20, "seed": 5, "traits": ["group"], **(extra or {})}
    else:
        networks["ppi"] = {"enabled": False}
    path = F.config(tmp_path / "data", files, groups=["A", "B"], contrasts=[B.contrast("B-A", "B", "A")], post_de={"enabled": True, "networks": networks})
    payload, code = workflow.run_command(path, tmp_path / "run")
    if expect is not None:
        assert code == expect, payload
    return payload, code, tmp_path / "run"


@pytest.fixture(scope="module")
def modules(tmp_path_factory):
    return coab_run(tmp_path_factory.mktemp("coab"))[2]


def _net(out):
    return out / "post_de" / "networks"


# ----------------------------------------------------------------------------- V161
def test_v161_small_cohort_refused_large_cohort_runs(tmp_path, modules):
    payload, code, out = coab_run(tmp_path, n_per_group=6)
    stage = next(s for s in payload["stages"] if s["stage_id"] == "post_de_networks")
    assert stage["state"] == "INAPPLICABLE" and stage["reason_code"] == "E_COABUNDANCE_SMALL_N" and not _net(out).exists()
    assert (_net(modules) / "modules.tsv").is_file()


def test_v161_repeated_design_is_refused_with_its_design_reason(tmp_path):
    values, obs = F.coabundance_design(n_per_group=30)
    for i, o in enumerate(obs):
        o["subject_id"] = f"S{i % 30 + 1}"                                        # every subject observed in both groups (paired)
    files = F.write_dataset(tmp_path / "data", values, obs)
    networks = {"enabled": True, "coabundance": {"enabled": True, "min_units": 20, "min_module_size": 8, "bootstrap": 10, "seed": 5}}
    path = F.config(tmp_path / "data", files, groups=["A", "B"], contrasts=[B.contrast("B-A", "B", "A")], post_de={"enabled": True, "networks": networks},
                    design_overrides={"blocking": {"mode": "duplicate_correlation", "subject_column": "subject_id"}})
    payload, code = workflow.run_command(path, tmp_path / "run")
    stage = next(s for s in payload["stages"] if s["stage_id"] == "post_de_networks")
    assert code == 0 and stage["state"] == "INAPPLICABLE" and stage["reason_code"] == "E_COABUNDANCE_DESIGN_UNSUPPORTED" and not _net(tmp_path / "run").exists()


def test_v161_negative_coabundance_is_never_on_by_default(tmp_path, ppi):
    with pytest.raises(ProteomicsError) as error:                                 # an enabled networks block with no declared analysis is refused
        coab_run(tmp_path, coabundance=False, expect=None)
    assert error.value.code == "E_NETWORKS_NOT_DECLARED" and not (tmp_path / "run" / "post_de").exists()
    out = ppi[0]                                                                  # a PPI-only declaration builds no modules
    assert (_net(out) / "connectivity.tsv").is_file() and not (_net(out) / "modules.tsv").exists()
    assert json.loads((_net(out) / "eligibility.json").read_text(encoding="utf-8"))["summary"].get("coabundance") is None


# ----------------------------------------------------------------------------- V162
def test_v162_modules_recover_planted_structure_with_unit_bootstrap(modules):
    rows = F.read_tsv(_net(modules) / "modules.tsv")
    found = {r["feature_id"]: r["module"] for r in rows}
    planted = F.planted_modules()
    from collections import Counter
    pairs = Counter((planted[f], found[f]) for f in planted)
    comb = lambda n: n * (n - 1) / 2
    n = len(planted); s = sum(comb(v) for v in pairs.values())
    a = sum(comb(v) for v in Counter(planted.values()).values()); b = sum(comb(v) for v in Counter(found.values()).values())
    expected = a * b / comb(n); ari = (s - expected) / (0.5 * (a + b) - expected)
    assert ari >= 0.9, ari
    eig = F.read_tsv(_net(modules) / "eigengenes.tsv")
    oracle = B.r_json("""a <- commandArgs(TRUE); Y <- as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE, encoding = 'UTF-8'))
      m <- read.delim(a[2], colClasses = 'character', encoding = 'UTF-8'); out <- list()
      for (mod in setdiff(unique(m$module), 'unassigned')) { Z <- scale(t(Y[m$feature_id[m$module == mod], ])); out[[mod]] <- unname(prcomp(Z, center = FALSE)$x[, 1]) }
      cat(jsonlite::toJSON(out, digits = NA))""", modules / "preprocessing" / "primary" / "matrix.tsv", _net(modules) / "modules.tsv")
    for mod, values in oracle.items():
        mine = [float(r[mod]) for r in eig]
        same = all(abs(x - y) <= 1e-8 for x, y in zip(mine, values)); flipped = all(abs(x + y) <= 1e-8 for x, y in zip(mine, values))
        assert same or flipped, mod                                                # first principal component up to sign
    stability = F.read_tsv(_net(modules) / "module_stability.tsv")
    assert {r["seed"] for r in stability} == {"5"} and {r["bootstrap"] for r in stability} == {"20"}
    assert {r["resampling"] for r in stability} == {"biological units with replacement"}          # never features (R-level hand-loop oracle in test-post-de-networks.R)
    assert {r["failed_draws"] for r in stability} == {"0"}
    assert all(float(r["mean_jaccard"]) > 0.5 for r in stability)
    beta = F.read_tsv(_net(modules) / "soft_threshold.tsv")
    chosen = {r["soft_threshold"] for r in rows}
    assert len(chosen) == 1 and any(r["beta"] == next(iter(chosen)) for r in beta)


# ----------------------------------------------------------------------------- V163
def test_v163_module_trait_equals_direct_limma_on_eigengenes(modules):
    rows = F.read_tsv(_net(modules) / "module_trait.tsv")
    assert {r["claim_label"] for r in rows} == {"module_level"} and {r["family_id"] for r in rows} == {"module_trait__group"}
    assert {r["moderation"] for r in rows} == {"limma eBayes trend=FALSE (eigengene means are 0), robust=TRUE"}
    oracle = B.r_json("""a <- commandArgs(TRUE); e <- read.delim(a[1], check.names = FALSE, encoding = 'UTF-8'); o <- read.delim(a[2], colClasses = 'character', encoding = 'UTF-8')
      o <- o[match(e$observation_id, o$observation_id), ]; E <- t(as.matrix(e[, grep('^M', names(e))])); colnames(E) <- e$observation_id
      eb <- limma::eBayes(limma::lmFit(E, model.matrix(~ factor(o$group, levels = c('A', 'B')))), trend = FALSE, robust = TRUE)
      cat(jsonlite::toJSON(list(m = rownames(E), t = unname(eb$t[, 2]), p = unname(eb$p.value[, 2]), q = p.adjust(eb$p.value[, 2], 'BH')), digits = NA))""",
                      _net(modules) / "eigengenes.tsv", modules / "preprocessing" / "primary" / "observations.tsv")
    by = {r["module"]: r for r in rows}
    for m, t, p, q in zip(oracle["m"], oracle["t"], oracle["p"], oracle["q"]):
        assert abs(float(by[m]["statistic"]) - t) <= 1e-8 and abs(float(by[m]["p_value"]) - p) <= 1e-10 and abs(float(by[m]["q_value"]) - q) <= 1e-10
    planted = F.planted_modules()
    module_of_m1 = next(r["module"] for r in F.read_tsv(_net(modules) / "modules.tsv") if r["feature_id"] == "F01")
    assert float(by[module_of_m1]["q_value"]) < 0.05                              # the planted group-associated module
    assert {r["family_id"] for r in F.read_tsv(modules / "dea" / "zero_null.tsv")} == {"protein-primary"}   # never pooled into protein families


def test_v163_phenotype_trait_equals_direct_limma_with_group_adjustment(tmp_path):
    import random
    values, obs = F.coabundance_design(n_per_group=30)
    rng = random.Random(5)
    for i, o in enumerate(obs):                                                   # a phenotype tracking module 3 (F25-F36)
        o["score"] = repr(sum(values[f"F{f:02d}"][i] for f in range(25, 37)) / 12 + rng.gauss(0, 0.3))
    files = F.write_dataset(tmp_path / "data", values, obs, extra_columns=("score",))
    networks = {"enabled": True, "coabundance": {"enabled": True, "min_units": 20, "rule": "wgcna_signed", "min_module_size": 8, "bootstrap": 10, "seed": 5, "traits": ["group", "score"]}}
    path = F.config(tmp_path / "data", files, groups=["A", "B"], contrasts=[B.contrast("B-A", "B", "A")], post_de={"enabled": True, "networks": networks})
    payload, code = workflow.run_command(path, tmp_path / "run")
    assert code == 0, payload
    out = tmp_path / "run"
    rows = [r for r in F.read_tsv(_net(out) / "module_trait.tsv") if r["trait"] == "score"]
    assert {r["family_id"] for r in rows} == {"module_trait__score"} and {r["claim_label"] for r in rows} == {"module_level"}
    oracle = B.r_json("""a <- commandArgs(TRUE); e <- read.delim(a[1], check.names = FALSE, encoding = 'UTF-8'); o <- read.delim(a[2], colClasses = 'character', encoding = 'UTF-8')
      o <- o[match(e$observation_id, o$observation_id), ]; E <- t(as.matrix(e[, grep('^M', names(e))])); colnames(E) <- e$observation_id
      s <- as.numeric(o$score); X <- model.matrix(~ factor(o$group, levels = c('A', 'B')) + I(s - mean(s)))
      eb <- limma::eBayes(limma::lmFit(E, X), trend = FALSE, robust = TRUE)
      cat(jsonlite::toJSON(list(m = rownames(E), t = unname(eb$t[, 3]), p = unname(eb$p.value[, 3]), q = p.adjust(eb$p.value[, 3], 'BH')), digits = NA))""",
                      _net(out) / "eigengenes.tsv", out / "preprocessing" / "primary" / "observations.tsv")
    by = {r["module"]: r for r in rows}
    for m, t, p, q in zip(oracle["m"], oracle["t"], oracle["p"], oracle["q"]):
        assert abs(float(by[m]["statistic"]) - t) <= 1e-8 and abs(float(by[m]["p_value"]) - p) <= 1e-10 and abs(float(by[m]["q_value"]) - q) <= 1e-10
    module_of_m3 = next(r["module"] for r in F.read_tsv(_net(out) / "modules.tsv") if r["feature_id"] == "F25")
    assert float(by[module_of_m3]["q_value"]) < 0.05


# ----------------------------------------------------------------------------- V164 / V165
def ppi_setup(tmp_path, *, tamper=False, url=False, null_universe=None):
    values, obs = F.biomarker_design(n_pos=6, n_neg=12, n_features=40, up=8, down=0, shift=2.0, seed=73)
    for o in [o for o in obs if o["group"] == "N"][6:]:
        o["group"] = "Q"                                                          # a third group without any effect: Q-N is a null contrast
    features = list(values)
    snap = F.write_ppi_snapshot(tmp_path / "sources", features, cluster=[f"F{i:02d}" for i in range(1, 9)])
    prepared = {k: resources.prepare(snap[f"{k}_manifest"], tmp_path / "snapshots" / k) for k in ("mapping", "ppi")}
    files = F.write_dataset(tmp_path / "data", values, obs)
    def resource(rid, kind, prep, extra=None):
        entry = {"id": rid, "kind": kind, "path": str(Path(prep["snapshot"]) / "manifest.json"), "sha256": prep["manifest_sha256"], "version": "synthetic-ppi-1",
                 "source": "synthetic test fixture (not a real database)", "source_taxonomy_id": 10116, "target_taxonomy_id": 10116, "terms": "synthetic; no redistribution restrictions"}
        entry.update(extra or {}); return entry
    ppi_entry = resource("ppi", "ppi", prepared["ppi"], {"score_type": "combined_score"})
    if tamper:
        ppi_entry["sha256"] = "1" * 64
    if url:
        ppi_entry["path"] = "https://example.invalid/ppi/manifest.json"
    ppi = {"enabled": True, "snapshot_id": "ppi", "mapping_resource_id": "map", "min_score": 300, "sets": ["dep", "rand"], "null_draws": 199, "seed": 9}
    if null_universe:
        ppi["null_universe"] = null_universe
    post = {"enabled": True, "sets": {"enabled": True, "definitions": [{"id": "dep", "rule": F.leaf("protein-primary", "P-N")},
                                                                         {"id": "rand", "rule": F.leaf("protein-secondary", "Q-N", criterion="exploratory_raw_p", threshold=0.5)}]},
            "networks": {"enabled": True, "ppi": ppi}}
    def mutate(c):
        c["resources"] = [resource("map", "mapping", prepared["mapping"]), ppi_entry]
    path = F.config(tmp_path / "data", files, groups=["P", "N", "Q"], contrasts=[B.contrast("P-N", "P", "N"), B.contrast("Q-N", "Q", "N", role="secondary")], post_de=post, mutate=mutate)
    return path, snap, prepared


@pytest.fixture(scope="module")
def ppi(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("ppi")
    path, snap, prepared = ppi_setup(tmp)
    with offline():                                                               # any network call during analysis would fail
        payload, code = workflow.run_command(path, tmp / "run")
    assert code == 0, payload
    return tmp / "run", snap, prepared


def test_v164_hashed_offline_snapshot_loads_with_its_provenance(ppi):
    out, snap, prepared = ppi
    manifest = Path(prepared["ppi"]["snapshot"]) / "manifest.json"
    assert hashlib.sha256(manifest.read_bytes()).hexdigest() == prepared["ppi"]["manifest_sha256"]
    summary = json.loads((_net(out) / "eligibility.json").read_text(encoding="utf-8"))["summary"]["ppi"]
    assert summary["source"] == "synthetic test fixture (not a real database)" and json.loads(manifest.read_text(encoding="utf-8"))["source"] == summary["source"]
    assert summary["release"] == "synthetic-ppi-1" and summary["species"] == 10116 and summary["score_type"] == "combined_score"
    assert summary["snapshot_sha256"] == prepared["ppi"]["manifest_sha256"]
    assert json.loads(manifest.read_text(encoding="utf-8"))["kind"] == "ppi"


def test_v164_negative_hash_mismatch_and_download_are_refused(tmp_path):
    path, _, _ = ppi_setup(tmp_path / "hash", tamper=True)
    payload, code = workflow.run_command(path, tmp_path / "hash" / "run")
    assert code == 2 and payload["error"]["code"] == "E_RESOURCE_HASH"
    path, _, _ = ppi_setup(tmp_path / "url", url=True)
    payload, code = workflow.run_command(path, tmp_path / "url" / "run")
    assert code == 2 and payload["error"]["code"] == "E_RESOURCE_DOWNLOAD"


def test_v165_connectivity_null_from_the_measured_mapped_universe(ppi):
    """Observed edges, set sizes and the universe are counted independently here; the null sampler itself is checked
    against exact enumeration on a tiny graph in test-post-de-networks.R (review 2026-10-05, minor 5)."""
    out, snap, _ = ppi
    rows = {r["set_id"]: r for r in F.read_tsv(_net(out) / "connectivity.tsv")}
    dep, rest = rows["dep"], rows["rand"]
    assert float(dep["p_value"]) == 1 / 200 and float(rest["p_value"]) > 0.1      # planted cluster at the floor; random-like set non-small
    oracle = B.r_json("""a <- commandArgs(TRUE); e <- read.delim(a[1], colClasses = 'character', encoding = 'UTF-8'); e <- e[as.numeric(e$score) >= 300 & e$gene_a != e$gene_b, ]
      e <- unique(data.frame(a = pmin(e$gene_a, e$gene_b), b = pmax(e$gene_a, e$gene_b)))
      mp <- read.delim(a[2], colClasses = 'character', encoding = 'UTF-8'); mem <- read.delim(a[3], colClasses = 'character', encoding = 'UTF-8')
      universe <- intersect(unique(mp$gene_id), unique(c(e$a, e$b))); acc <- setNames(mp$gene_id, mp$source_id)
      res <- list(universe = length(universe))
      for (s in c('dep', 'rand')) { genes <- intersect(unique(na.omit(acc[mem$accession_key[mem[[s]] == 'true']])), universe)
        res[[s]] <- list(obs = sum(e$a %in% genes & e$b %in% genes), n = length(genes)) }
      cat(jsonlite::toJSON(res, auto_unbox = TRUE, digits = NA))""",
                      Path(snap["ppi_manifest"]).parent / "ppi_source.tsv", Path(snap["mapping_manifest"]).parent / "mapping_source.tsv", _accession_membership(out, snap))
    for s in ("dep", "rand"):
        assert int(rows[s]["observed_edges"]) == oracle[s]["obs"] and int(rows[s]["n_set_genes_in_universe"]) == oracle[s]["n"]
        assert abs(float(rows[s]["p_value"]) - (int(rows[s]["k"]) + 1) / (int(rows[s]["null_draws"]) + 1)) <= 1e-15
        assert rows[s]["seed"] == "9"
    assert (dep["stream"], rest["stream"]) == ("1", "2")                          # one derived L'Ecuyer stream per set, never a shared seed
    assert int(dep["universe_size"]) == oracle["universe"] and dep["universe"] == "measured mapped genes in the snapshot"
    assert int(dep["universe_size"]) < len(snap["genes"]) + len(snap["unmeasured"])                 # not the whole snapshot
    hubs = F.read_tsv(_net(out) / "hub_degree.tsv")
    assert {h["claim_label"] for h in hubs} == {"descriptive"}


def _accession_membership(out, snap):
    """Membership of the R14a sets with each feature's accession (the oracle maps accessions itself)."""
    rows = F.read_tsv(out / "post_de" / "sets" / "membership.tsv")
    accession = {f: F.SPECIES["rat"]["accession"](i) for i, f in enumerate(snap["genes"], start=1)}
    path = out.parent / "membership_with_accessions.tsv"
    B.write_tsv(path, [["feature_id", "accession_key", "dep", "rand"]] + [[r["feature_id"], accession[r["feature_id"]], r["dep"], r["rand"]] for r in rows])
    return path


def test_v165_negative_genome_wide_null_is_refused(tmp_path):
    path, _, _ = ppi_setup(tmp_path, null_universe="snapshot")
    with pytest.raises(ProteomicsError) as error:
        workflow.run_command(path, tmp_path / "run")
    assert error.value.code == "E_NETWORK_UNIVERSE"


def test_r14e_report_section(modules, ppi):
    for out in (modules, ppi[0]):
        data = json.loads((out / "report" / "report_data.json").read_text(encoding="utf-8"))
        section = data["sections"]["post_de_networks"]
        assert section["state"] == "COMPLETED"
        for figure in section["figures"]:
            assert (out / figure["source"]).is_file() and (out / figure["src"]).is_file()
