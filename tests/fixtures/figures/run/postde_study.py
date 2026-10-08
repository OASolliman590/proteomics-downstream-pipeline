"""A synthetic post-DE study for the run-level figure catalogue (R16e, D-81). Built at test time under the caller's directory
(the resource snapshots hold absolute paths, so nothing is committed).

It reuses the synthetic post-DE factory (tests/fixtures/post_de/design_factory.py): the two-group study with planted DEPs,
phenotype-linked features, a co-abundance module and a PPI snapshot. The post-DE block is cheap: a small biomarker budget
(outer 3 folds, one repeat, B 19), the association and the sensitivity modules off, the co-abundance module on, and a gene-set
enrichment on the planted sets (gene-set resource of the pipeline; the synthetic sets are not a GMT file).
"""
from __future__ import annotations

import copy
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


F = _load("post_de_factory", ROOT / "tests" / "fixtures" / "post_de" / "design_factory.py")
SETS = {"SET_UP": range(1, 6), "SET_SCORE": range(6, 11), "SET_COMODULE": range(11, 21), "SET_NULL": range(21, 31)}


def _prepare_gene_sets(directory: Path, species: str) -> dict:
    from proteomics_pipeline import resources
    gene = F.SPECIES[species]["gene"]
    rows = ["set_id\tset_name\tgene_id"] + [f"{s}\tsynthetic {s}\t{gene(i)}" for s, idx in SETS.items() for i in idx]
    (directory / "genesets_source.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8", newline="\n")
    tax = F.SPECIES[species]["taxonomy_id"]
    base = {"version": "synthetic-1", "source": "synthetic gene-set fixture (not a real database)", "terms": "synthetic; no redistribution restrictions",
            "id_type": "synthetic", "source_taxonomy_id": tax, "target_taxonomy_id": tax}
    manifest = directory / "prepare_sets.json"
    manifest.write_text(json.dumps({**base, "resource_id": "sets", "kind": "gene_sets", "files": [{"name": "gene_sets.tsv", "source": "genesets_source.tsv"}]}),
                        encoding="utf-8", newline="\n")
    prep = resources.prepare(manifest, directory / "snapshots" / "sets")
    return {"id": "sets", "kind": "gene_sets", "path": str(Path(prep["snapshot"]) / "manifest.json"), "sha256": prep["manifest_sha256"],
            "version": base["version"], "source": base["source"], "source_taxonomy_id": tax, "target_taxonomy_id": tax, "terms": base["terms"]}


def study_config(directory: Path, *, species: str = "rat") -> Path:
    """The configuration file of the study, written to ``directory`` with its data and resources."""
    from proteomics_pipeline import resources
    cell = F.generality_cell("two_group")
    data = directory / "data"
    files = F.write_dataset(data, cell["values"], cell["observations"], species=cell["species"], extra_columns=cell["extra_columns"])
    snap = F.write_ppi_snapshot(directory / "sources", list(cell["values"]), cluster=[f"F{i:02d}" for i in range(1, 6)], species=cell["species"])
    resource_entries = []
    for kind in ("mapping", "ppi"):
        prep = resources.prepare(snap[f"{kind}_manifest"], directory / "snapshots" / kind)
        entry = {"id": "map" if kind == "mapping" else "ppi", "kind": kind, "path": str(Path(prep["snapshot"]) / "manifest.json"), "sha256": prep["manifest_sha256"],
                 "version": "synthetic-ppi-1", "source": "synthetic test fixture (not a real database)", "source_taxonomy_id": F.SPECIES[cell["species"]]["taxonomy_id"],
                 "target_taxonomy_id": F.SPECIES[cell["species"]]["taxonomy_id"], "terms": "synthetic; no redistribution restrictions"}
        if kind == "ppi":
            entry["score_type"] = "combined_score"
        resource_entries.append(entry)
    resource_entries.append(_prepare_gene_sets(directory / "sources", cell["species"]))
    block = F.generality_post_de(cell, ppi=cell["species"])
    block["sensitivity"]["enabled"] = False
    block["association"]["enabled"] = False
    block["biomarker"]["cv"] = {"outer": {"k": 3, "repeats": 1}, "inner": {"k": 3}}
    block["biomarker"]["permutation"] = {"B": 19}
    block["biomarker"]["single_feature"] = {"bootstrap": 100}
    block["networks"]["coabundance"] = {"enabled": True, "min_units": 20, "min_module_size": 5, "bootstrap": 10, "seed": 4, "traits": ["group"]}

    from proteomics_pipeline import pathway_service
    contrast_ids = [c["id"] for c in cell["contrasts"]]

    def change(config: dict) -> None:
        config["runtime"]["threads"] = 2
        config["report"]["figures"] = {"formats": ["png", "pzfx", "html"], "dpi": 300, "journal_width_mm": 120, "palette": "okabe_ito", "font": "Arial", "p_annotation": "exact"}
        config["resources"] = resource_entries
        for hypothesis in pathway_service.METHOD_TYPES["ora"]:
            config["multiplicity_families"].append({"id": f"{hypothesis}-fam", "hypothesis_type": hypothesis, "model_ids": ["limma-main"],
                                                    "contrast_ids": contrast_ids, "collection_ids": ["sets"], "adjustment": "BH",
                                                    "denominator": "finite_eligible", "q_cutoff": 0.05, "role": "secondary"})
        config["pathways"] = {"enabled": True, "methods": ["ora"], "mapping_resource_id": "map", "gene_set_resource_ids": ["sets"], "min_size": 3,
                              "max_size": 20, "representative_rule": "coverage_median_stable_id", "multi_gene_policy": "exclude",
                              "gene_aggregation_sensitivity": "off", "nrot": 10,
                              "midp": False, "rank": "native_zero_null", "eps": 0, "ora_foreground": "family_q", "ora_threshold": 0.05}

    return F.config(data, files, groups=cell["groups"], contrasts=cell["contrasts"], post_de=block, design_overrides=cell["design_overrides"],
                    mutate=change, species=cell["species"])
