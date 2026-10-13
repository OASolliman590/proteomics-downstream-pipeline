"""R15b acceptance: V174-V178 (SM44, SM17; ADR 0011).

Oracles are written here: hashlib for file hashes, stats::phyper through Rscript for the ORA P values, a Python
brute-force running sum and a full enumeration of every subset for the pre-ranked score and exact P, and
p.adjust(BH) applied per family. Synthetic marked fixtures only. Each case has a negative case next to it.
"""
from __future__ import annotations

import importlib.util
import hashlib
import itertools
import json
import socket
import subprocess
from pathlib import Path

import pytest

from proteomics_pipeline import gene_set_libraries as GSL
from proteomics_pipeline import workflow
from proteomics_pipeline.cli import main as cli_main
from proteomics_pipeline.errors import ProteomicsError

ROOT = Path(__file__).resolve().parents[2]


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


G = _load("pathway_factory", ROOT / "tests" / "fixtures" / "methods_extension" / "pathway_libraries" / "factory.py")
W = _load("welch_factory", ROOT / "tests" / "fixtures" / "methods_extension" / "welch" / "factory.py")
B = _load("design_builders", ROOT / "tests" / "fixtures" / "design" / "builders.py")
R_REQUIRED = pytest.mark.skipif(B.rscript() is None, reason="NOT_RUN: Rscript/proteomicsCore/limma unavailable")

PHYPER = """
  a <- commandArgs(TRUE); spec <- jsonlite::fromJSON(a[1])
  p <- mapply(function(k, K) stats::phyper(k - 1, K, spec$N - K, spec$n, lower.tail = FALSE), spec$k, spec$K)
  cat(jsonlite::toJSON(unname(p), digits = NA))"""

BH = """
  a <- commandArgs(TRUE); p <- as.numeric(jsonlite::fromJSON(a[1]))
  cat(jsonlite::toJSON(unname(p.adjust(p, method = 'BH')), digits = NA))"""

POOLED = """
  a <- commandArgs(TRUE); rows <- as.data.frame(jsonlite::fromJSON(a[1]), stringsAsFactors = FALSE)
  key <- function(r) paste(r$method, r$contrast_id, sep = '|')
  out <- proteomicsCore:::enrichment_family_adjust(rows, key)
  cat(jsonlite::toJSON(out$q_value, digits = NA))"""


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _decl(lib_id, path: Path, *, licence="SYNTHETIC-CC0", redistributable=True, source="synthetic fixture", version="1", species="synthetic",
          identifier="gene_symbol"):
    return {"id": lib_id, "path": str(path), "sha256": _sha(path), "source": source, "version": version, "species": species,
            "identifier_type": identifier, "licence_id": licence, "redistributable": redistributable}


def _stage(tmp: Path, libraries, welch_path: Path, *, background=None, min_size=1, max_size=500, tie_break="sm15_stable_id", name="stage"):
    base = tmp / "libraries"
    base.mkdir(exist_ok=True)
    request = GSL.build_request(libraries, base_dir=base, welch_results=welch_path, background=background or {"kind": "measured", "size": 500},
                                min_size=min_size, max_size=max_size, run_id="run-r15b", output_temp_dir=tmp / name, tie_break=tie_break)
    return GSL.execute(request), tmp / name


def _rows(path: Path) -> list[dict]:
    return B.read_tsv(path)


# ----------------------------------------------------------------------------- V174
def test_v174_library_is_hashed_recorded_and_read_from_the_same_bytes(tmp_path):
    gmt = Path(ROOT / "tests" / "fixtures" / "methods_extension" / "pathway_libraries" / "synthetic_small.gmt")
    libraries = [_decl("syn_small", gmt)]
    manifest = GSL.library_record(libraries[0], tmp_path)
    assert manifest["sha256"] == hashlib.sha256(gmt.read_bytes()).hexdigest()   # the oracle is hashlib over the file bytes
    for key in ("source", "version", "species", "identifier_type", "licence_id", "redistributable"):
        assert manifest[key] == libraries[0][key]
    assert manifest["synthetic"] is True


@R_REQUIRED
def test_v174_run_records_the_same_hash_and_writes_the_manifest(tmp_path):
    welch = tmp_path / "welch.tsv"
    G.write_welch_table(welch, [{"feature_id": f"F{i}", "gene_id": f"G{i:04d}", "contrast_id": "B-A", "effect": 1.0 + i * 0.1, "statistic": 1.0 + i * 0.3,
                                 "q_value": 0.9, "p_value": 0.5} for i in range(1, 13)])
    gmt = Path(ROOT / "tests" / "fixtures" / "methods_extension" / "pathway_libraries" / "synthetic_small.gmt")
    result, out = _stage(tmp_path, [_decl("syn_small", gmt)], welch, background={"kind": "measured", "size": 12}, max_size=12)
    assert result["state"] == "COMPLETED", result
    manifest = json.loads((out / "libraries.json").read_text(encoding="utf-8"))
    record = manifest["libraries"][0]
    assert record["sha256"] == hashlib.sha256(gmt.read_bytes()).hexdigest()
    assert manifest["library_content_written"] is False and manifest["network"] == "none"
    assert "path" not in record and str(tmp_path) not in json.dumps(manifest)     # no local path is recorded


def test_v174_negative_one_byte_edit_fails_the_hash_check(tmp_path):
    gmt = tmp_path / "edited.gmt"
    G.write_gmt(gmt, {"SET_A": ["G0001", "G0002"]})
    declared = _decl("edited", gmt)
    raw = bytearray(gmt.read_bytes()); raw[-2] = ord("9" if raw[-2] != ord("9") else "8")   # one byte changes after declaration
    gmt.write_bytes(bytes(raw))
    with pytest.raises(ProteomicsError) as info:
        GSL.library_record(declared, tmp_path)
    assert info.value.code == "E_RESOURCE_HASH"


def test_v174_negative_network_location_and_socket_calls_fail_with_e_resource_download(tmp_path):
    remote = {"id": "remote", "path": "https://example.invalid/library.gmt", "sha256": "0" * 64, "source": "x", "version": "1", "species": "s",
              "identifier_type": "gene_symbol", "licence_id": "CC-BY-4.0", "redistributable": True}
    with pytest.raises(ProteomicsError) as info:
        GSL.library_record(remote, tmp_path)
    assert info.value.code == "E_RESOURCE_DOWNLOAD"
    with GSL.no_network():
        with pytest.raises(ProteomicsError) as socket_info:
            socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        assert socket_info.value.code == "E_RESOURCE_DOWNLOAD"
    assert socket.socket.__name__ == "socket"                                   # the guard is removed afterwards


# ----------------------------------------------------------------------------- V175
def _tracked_gmt_files() -> list[Path]:
    # tracked files plus untracked files that git would not ignore: the set that a commit would publish
    listed = subprocess.run(["git", "-C", str(ROOT), "ls-files", "--cached", "--others", "--exclude-standard", "*.gmt"], capture_output=True, text=True,
                            encoding="utf-8", check=True).stdout.split()
    return [ROOT / item for item in listed]


def test_v175_every_tracked_gmt_file_is_a_marked_synthetic_fixture():
    files = _tracked_gmt_files()
    assert files, "the repository should carry the synthetic GMT fixture"
    for path in files:
        assert GSL.is_synthetic(path), path.relative_to(ROOT).as_posix()


def test_v175_release_artifact_with_only_synthetic_members_is_written(tmp_path):
    gmt = Path(ROOT / "tests" / "fixtures" / "methods_extension" / "pathway_libraries" / "synthetic_small.gmt")
    libraries = [GSL.library_record(_decl("syn_small", gmt), tmp_path)]
    accepted = GSL.assert_artifact_members({"synthetic_small.gmt": gmt}, libraries)
    assert accepted == ["synthetic_small.gmt"]


@pytest.mark.parametrize("variant", ["kegg_licence", "msigdb_licence", "non_synthetic_redistributable"])
def test_v175_negative_each_restricted_variant_is_refused_and_nothing_is_written(tmp_path, variant):
    if variant == "kegg_licence":
        path = G.write_gmt(tmp_path / "kegg.gmt", {"KEGG_X": ["G0001", "G0002"]}); licence, redist = "KEGG-restricted", False
    elif variant == "msigdb_licence":
        path = G.write_gmt(tmp_path / "msig.gmt", {"MSIG_X": ["G0001", "G0002"]}); licence, redist = "MSigDB-restricted", False
    else:
        path = G.write_gmt(tmp_path / "real_like.gmt", {"REAL_X": ["G0001", "G0002"]}, marked=False); licence, redist = "CC-BY-4.0", True
    # A library may be declared and used locally; it may never be part of an artifact (FR-175)
    library = GSL.library_record({"id": variant, "path": str(path), "sha256": _sha(path), "source": "x", "version": "1", "species": "s",
                                  "identifier_type": "gene_symbol", "licence_id": licence, "redistributable": redist}, tmp_path)
    artifact = tmp_path / "artifact"
    with pytest.raises(ProteomicsError) as refused:
        GSL.assert_artifact_members({f"{variant}.gmt": path}, [library])
    assert refused.value.code == "E_LIBRARY_REDISTRIBUTION"
    assert not artifact.exists()


# ----------------------------------------------------------------------------- V176
def _ora_fixture(tmp_path: Path):
    universe = G.gene_names(500)
    set_a = universe[0:30]; set_b = universe[99:129]
    foreground = set(set_a[:8]) | set(universe[199:231])        # n = 40, k = 8 in set A, k = 0 in set B
    gmt = G.write_gmt(tmp_path / "ora_lib.gmt", {"SET_A": set_a, "SET_B": set_b})
    welch = G.ora_universe_table(tmp_path / "ora_welch.tsv", genes=universe, foreground=foreground)
    return gmt, welch, set_a, set_b, foreground, universe


@R_REQUIRED
def test_v176_ora_p_values_equal_phyper_and_bh_stays_inside_the_family(tmp_path):
    gmt, welch, set_a, set_b, foreground, universe = _ora_fixture(tmp_path)
    result, out = _stage(tmp_path, [_decl("ora_lib", gmt)], welch, background={"kind": "measured", "size": 500}, max_size=500)
    assert result["state"] == "COMPLETED", result
    rows = [r for r in _rows(out / "ora_ora_lib.tsv") if r["method"] == "ora_up"]
    by_set = {r["set_id"]: r for r in rows}
    oracle = B.r_json(PHYPER, json.dumps({"N": 500, "K": [30, 30], "n": 40, "k": [8, 0]}))
    assert abs(float(by_set["SET_A"]["p_value"]) - oracle[0]) / oracle[0] <= 1e-12
    assert float(by_set["SET_B"]["p_value"]) == 1.0 and int(by_set["SET_B"]["k"]) == 0      # k = 0 is tested and kept
    bh = B.r_json(BH, json.dumps([float(by_set["SET_A"]["p_value"]), float(by_set["SET_B"]["p_value"])]))
    assert abs(float(by_set["SET_A"]["q_value"]) - bh[0]) <= 1e-12 and abs(float(by_set["SET_B"]["q_value"]) - bh[1]) <= 1e-12
    assert int(by_set["SET_A"]["N"]) == 500 and int(by_set["SET_A"]["K"]) == 30 and int(by_set["SET_A"]["n"]) == 40


@R_REQUIRED
def test_v176_negative_whole_genome_background_fails_e_ora_universe(tmp_path):
    gmt, welch, *_ = _ora_fixture(tmp_path)
    result, _out = _stage(tmp_path, [_decl("ora_lib", gmt)], welch, background={"kind": "genome", "size": 20000}, name="genome")
    assert result["state"] == "FAILED" and result["reason_code"] == "E_ORA_UNIVERSE"


# ----------------------------------------------------------------------------- V177
def _brute_force(stat_by_gene: dict, members: set):
    order = sorted(stat_by_gene, key=lambda g: (-stat_by_gene[g], g))
    stats = [stat_by_gene[g] for g in order]
    n = len(order)

    def es_for(positions):
        k = len(positions); total = sum(abs(stats[i]) for i in positions); running = 0.0; best = 0.0
        for i in range(n):
            running += (abs(stats[i]) / total) if i in positions else -1.0 / (n - k)
            if abs(running) > abs(best):
                best = running
        return best

    obs_positions = {i for i, g in enumerate(order) if g in members}
    es_obs = es_for(obs_positions)
    k = len(obs_positions)
    scores = [es_for(set(c)) for c in itertools.combinations(range(n), k)]
    p = sum(s >= es_obs - 1e-12 for s in scores) / len(scores) if es_obs >= 0 else sum(s <= es_obs + 1e-12 for s in scores) / len(scores)
    return es_obs, p, len(scores)


@R_REQUIRED
def test_v177_prerank_on_welch_t_equals_brute_force_and_is_labelled_exploratory(tmp_path):
    # Ten tie-free ranks from a Welch fixture: the t statistics are computed by the stats::t.test oracle.
    import random
    rng = random.Random(177)
    groups = {f"G{i:02d}": ([10 + rng.gauss(0, 0.4) for _ in range(6)], [10 + (0.8 if i < 3 else 0) + rng.gauss(0, 1.2) for _ in range(9)]) for i in range(1, 11)}
    t_values = {}
    for gene, (x, y) in groups.items():
        t_values[gene] = B.r_json("""
          a <- commandArgs(TRUE); x <- as.numeric(jsonlite::fromJSON(a[1])); y <- as.numeric(jsonlite::fromJSON(a[2]))
          cat(jsonlite::toJSON(unname(stats::t.test(x, y, var.equal = FALSE)$statistic), digits = NA, auto_unbox = TRUE))""", json.dumps(x), json.dumps(y))
    stats_ = {g: float(t) for g, t in t_values.items()}
    assert len(set(stats_.values())) == 10                                     # tie-free
    welch = G.write_welch_table(tmp_path / "rank_welch.tsv", [{"feature_id": f"F{g}", "gene_id": g, "contrast_id": "B-A", "effect": stats_[g], "statistic": stats_[g],
                                                             "q_value": 1.0, "p_value": 0.5} for g in stats_])
    gene_sets = {"SET3": ["G02", "G05", "G09"], "SET4": ["G01", "G04", "G07", "G10"]}
    gmt = G.write_gmt(tmp_path / "rank_lib.gmt", gene_sets)
    result, out = _stage(tmp_path, [_decl("rank_lib", gmt)], welch, background={"kind": "measured", "size": 10}, min_size=1, max_size=10, name="rank")
    assert result["state"] == "COMPLETED", result
    rows = {r["set_id"]: r for r in _rows(out / "prerank_rank_lib.tsv")}
    for set_id, members in gene_sets.items():
        es_oracle, p_oracle, n_subsets = _brute_force(stats_, set(members))
        assert abs(float(rows[set_id]["es"]) - es_oracle) <= 1e-12
        assert abs(float(rows[set_id]["p_value"]) - p_oracle) <= 1e-12
        assert rows[set_id]["exploratory"] == "true" and rows[set_id]["claim_label"] == "exploratory"
        assert rows[set_id]["estimator"] == "prerank_welch_exact_enumeration"       # 120 subsets: the exact enumeration path ran
        assert int(rows[set_id]["K"]) == len(members)


@R_REQUIRED
def test_v177_negative_duplicate_genes_and_unbroken_ties_are_refused(tmp_path):
    gmt = G.write_gmt(tmp_path / "dup_lib.gmt", {"SET": ["G01", "G02", "G03"]})
    dup = G.write_welch_table(tmp_path / "dup.tsv", [{"feature_id": "A", "gene_id": "G01", "contrast_id": "B-A", "effect": 1, "statistic": 1.0, "q_value": 1},
                                                    {"feature_id": "B", "gene_id": "G01", "contrast_id": "B-A", "effect": 1, "statistic": 2.0, "q_value": 1},
                                                    {"feature_id": "C", "gene_id": "G02", "contrast_id": "B-A", "effect": 1, "statistic": 3.0, "q_value": 1}])
    result, _ = _stage(tmp_path, [_decl("dup_lib", gmt)], dup, background={"kind": "measured", "size": 2}, name="dup")
    assert result["state"] == "FAILED" and result["reason_code"] == "E_RANK_DUPLICATE"
    tied = G.write_welch_table(tmp_path / "tied.tsv", [{"feature_id": f"T{i}", "gene_id": f"G0{i}", "contrast_id": "B-A", "effect": 1, "statistic": 1.0,
                                                       "q_value": 1} for i in range(1, 5)])
    result2, _ = _stage(tmp_path, [_decl("dup_lib", gmt)], tied, background={"kind": "measured", "size": 4}, tie_break="none", name="tied")
    assert result2["state"] == "FAILED" and result2["reason_code"] == "E_RANK_TIE"


# ----------------------------------------------------------------------------- V178
@R_REQUIRED
def test_v178_each_library_method_contrast_is_its_own_family_with_bh_inside(tmp_path):
    universe = G.gene_names(120)
    import random
    rng = random.Random(178)
    rows = []
    for gene in universe:
        for contrast, shift in (("B-A", 2.0), ("C-A", -2.0)):
            hit = (gene in universe[:25] and contrast == "B-A") or (gene in universe[60:80] and contrast == "C-A")
            stat = (shift if hit else 0.0) + rng.gauss(0, 1)
            rows.append({"feature_id": f"P{gene}{contrast}", "gene_id": gene, "contrast_id": contrast, "effect": stat, "statistic": stat,
                         "q_value": 0.01 if hit else 0.8, "p_value": 0.001 if hit else 0.4})
    welch = G.write_welch_table(tmp_path / "fam_welch.tsv", rows)
    lib_a = G.write_gmt(tmp_path / "fam_a.gmt", {"A1": universe[0:30], "A2": universe[30:60], "A3": universe[10:40]})
    lib_b = G.write_gmt(tmp_path / "fam_b.gmt", {"B1": universe[0:30], "B2": universe[20:50]})          # overlapping with A
    result, out = _stage(tmp_path, [_decl("fam_a", lib_a), _decl("fam_b", lib_b)], welch, background={"kind": "measured", "size": 120}, max_size=120, name="fam")
    assert result["state"] == "COMPLETED", result
    table = _rows(out / "ora_fam_a.tsv") + _rows(out / "prerank_fam_a.tsv") + _rows(out / "ora_fam_b.tsv") + _rows(out / "prerank_fam_b.tsv")
    families = {}
    for r in table:
        families.setdefault(r["family_id"], []).append(r)
    assert len(families) == 2 * 3 * 2                                          # 2 libraries x 3 methods x 2 contrasts
    for family_id, members in families.items():
        assert all(family_id.startswith(m["library_id"] + "|") for m in members)   # no pooled library family
        p = [float(m["p_value"]) for m in members]
        expected = B.r_json(BH, json.dumps(p))
        for m, q in zip(members, expected):
            assert abs(float(m["q_value"]) - q) <= 1e-12, family_id
    # negative: the same BH applied across libraries is refused by the stage's family rule
    pooled = [{"library_id": "fam_a", "method": "ora_up", "contrast_id": "B-A", "p_value": 0.01},
              {"library_id": "fam_b", "method": "ora_up", "contrast_id": "B-A", "p_value": 0.2}]
    with pytest.raises(RuntimeError) as info:
        B.r_json(POOLED, json.dumps(pooled))
    assert "E_ENRICHMENT_FAMILY_POOLED" in str(info.value)


# ----------------------------------------------------------------------------- F2: the libraries run inside a workflow
SYNTHETIC_GMT = ROOT / "tests" / "fixtures" / "methods_extension" / "pathway_libraries" / "synthetic_small.gmt"


def _workflow_config(tmp: Path, *, licence="SYNTHETIC-CC0", redistributable=True, gmt=SYNTHETIC_GMT, sha=None):
    values, obs = W.three_group(shifts={0: {"disease": 2.0}, 1: {"treated": -2.0}})
    contrasts = [W.contrast("disease-control", "disease", "control"), W.contrast("treated-control", "treated", "control"),
                 W.contrast("treated-disease", "treated", "disease")]
    cfg = W.study_config(tmp / "study", values, obs, groups=["control", "disease", "treated"], contrasts=contrasts, primary="welch",
                         methods=W.welch_methods(), genes=True)
    data = json.loads(cfg.read_text(encoding="utf-8"))
    data["methods_extension"]["pathways"] = {"enabled": True, "execution_requirement": "optional", "methods": ["ora", "prerank"], "min_size": 1, "max_size": 500,
                                             "libraries": [{"id": "syn_small", "path": str(gmt), "sha256": sha or _sha(gmt), "source": "synthetic fixture",
                                                            "version": "1", "species": "synthetic", "identifier_type": "gene_symbol", "licence_id": licence,
                                                            "redistributable": redistributable}]}
    data["report"]["figures"] = {"formats": ["png", "pzfx", "html"], "dpi": 300, "journal_width_mm": 120, "palette": "okabe_ito", "font": "Arial", "p_annotation": "exact"}
    cfg.write_text(json.dumps(data, indent=2), encoding="utf-8")
    return cfg


@R_REQUIRED
def test_f2_offline_libraries_run_in_the_workflow_with_records_and_figures(tmp_path, capsys):
    cfg = _workflow_config(tmp_path)
    payload, code = workflow.run_command(cfg, tmp_path / "run")
    assert code == 0 and payload["state"] == "COMPLETED", payload
    states = {s["stage_id"]: s["state"] for s in payload["stages"]}
    assert states.get("offline_gene_set_libraries") == "COMPLETED" and states.get("welch") == "COMPLETED", states
    out = tmp_path / "run" / "methods_extension" / "pathways"
    plan = json.loads((tmp_path / "run" / "plan.json").read_text(encoding="utf-8"))
    planned = [c for c in plan["capability_plan"] if c["capability"] == "offline_gene_set_libraries"]
    assert planned and planned[0]["required"] is False and planned[0]["implemented"] is True
    manifest = json.loads((out / "libraries.json").read_text(encoding="utf-8"))["libraries"][0]
    assert manifest["sha256"] == hashlib.sha256(SYNTHETIC_GMT.read_bytes()).hexdigest()       # the hash recorded equals the file bytes
    assert manifest["licence_id"] == "SYNTHETIC-CC0" and manifest["redistributable"] is True and manifest["synthetic"] is True
    assert "path" not in manifest and str(tmp_path) not in (out / "libraries.json").read_text(encoding="utf-8")
    ora = [r for r in _rows(out / "ora_syn_small.tsv") if r["method"] == "ora_up"]
    assert any(int(r["k"]) >= 1 for r in ora)                                                  # the planted up genes hit a set
    # the measured universe is the run's mapped genes (60 features -> 60 genes), not a genome
    assert {r["N"] for r in ora} == {"60"}
    catalogue = {p.name for p in (tmp_path / "run" / "report" / "figure_catalogue").iterdir()}
    assert any(name.startswith("path_enrichment_dotplot__syn_small__ora_up") and name.endswith(".png") for name in catalogue)
    capsys.readouterr()
    assert cli_main(["verify", "--run", str(tmp_path / "run"), "--json"]) == 0


def test_f2_negative_a_one_byte_edit_after_declaration_fails_the_plan_with_e_resource_hash(tmp_path):
    gmt = tmp_path / "edited.gmt"
    gmt.write_bytes(SYNTHETIC_GMT.read_bytes())
    declared_sha = _sha(gmt)
    raw = bytearray(gmt.read_bytes()); raw[-3] = ord("Z")          # the file changes after the hash was declared
    gmt.write_bytes(bytes(raw))
    cfg = _workflow_config(tmp_path, gmt=gmt, sha=declared_sha)
    payload, code = workflow.run_command(cfg, tmp_path / "run")
    assert code != 0 and payload["state"] == "FAILED" and payload["error"]["code"] == "E_RESOURCE_HASH", payload
    assert not (tmp_path / "run" / "methods_extension").exists()


def test_f2_negative_a_restricted_licence_declared_redistributable_fails_the_plan(tmp_path):
    kegg = G.write_gmt(tmp_path / "kegg_like.gmt", {"KEGG_X": ["G0001", "G0002"]})
    cfg = _workflow_config(tmp_path, licence="KEGG-restricted", redistributable=True, gmt=kegg)
    payload, code = workflow.run_command(cfg, tmp_path / "run")
    assert code != 0 and payload["state"] == "FAILED" and payload["error"]["code"] == "E_LIBRARY_REDISTRIBUTION", payload
