"""Figure registry (packet R16a, FR-197, V197).

Each catalogue figure is registered once with its figure-source table. The table is copied to
``sources/<figure_id>.tsv`` and its SHA-256 is recorded. An output is produced only through ``produce``, which first
requires the registered source (E_FIGURE_SOURCE_MISSING otherwise, with nothing written) and records the source hash that
the output was generated from. ``registry.json`` also holds the resolved style of the run.
"""
from __future__ import annotations

import json
from collections.abc import Callable, Mapping
from pathlib import Path
from typing import Any

from ..errors import ConfigurationError
from ..provenance import sha256_file
from .config import FIGURE_FORMATS

# Catalogue of contract figures.md, "Figure IDs": id -> (packet, area, source family, claim or method label).
CATALOGUE: dict[str, tuple[str, str, str, str]] = {
    "qc_intensity_distributions": ("R16e", "QC", "observed log2 intensity (SM06)", "descriptive"),
    "qc_missingness_heatmap": ("R16e", "QC", "observed mask (SM04)", "descriptive"),
    "qc_sample_correlation_heatmap": ("R16e", "QC", "SM06 coverage rule", "descriptive"),
    "qc_pca_ellipses": ("R16e", "QC", "PCA on observed cells", "descriptive"),
    "diff_volcano_labelled": ("R16e", "Differential", "zero-null family (SM12)", "family q"),
    "diff_ma_plot": ("R16e", "Differential", "contrast means", "descriptive"),
    "diff_protein_dotplot_brackets": ("R16e", "Differential", "family P (FR-198)", "family P"),
    "diff_dep_heatmap_clustered": ("R16e", "Differential", "declared DEP set", "descriptive"),
    "diff_upset_venn": ("R16e", "Differential", "SM31 exclusive regions (Venn for k <= 3)", "descriptive"),
    "bio_roc_ci_band": ("R16e", "Biomarker", "SM34", "claim label"),
    "bio_panel_forest": ("R16e", "Biomarker", "SM34", "claim label"),
    "bio_cv_auc_vs_permutation": ("R16e", "Biomarker", "SM36", "claim label"),
    "path_enrichment_dotplot": ("R16e", "Pathways", "SM18", "family q"),
    "net_module_heatmap": ("R16e", "Networks", "SM39", "module_level"),
    "net_network_graph": ("R16e", "Networks", "SM40", "descriptive"),
    "perm_pcoa_centroids": ("R16e", "PERMANOVA", "SM26", "descriptive"),
    "perm_permdisp": ("R16e", "PERMANOVA", "SM29", "family P"),
}
SOURCE_DIR = "sources"
REGISTRY_FILE = "registry.json"


class FigureRegistry:
    """Registry of one run's figure catalogue, written under ``<run>/report/figure_catalogue/``."""

    def __init__(self, directory: str | Path, style: Mapping[str, Any]):
        self.directory = Path(directory)
        self.style = dict(style)
        self.figures: dict[str, dict[str, Any]] = {}

    def register(self, figure_id: str, source_table: str | Path) -> dict[str, Any]:
        """Register a catalogue figure with its source table (copied into the catalogue, hashed)."""
        if figure_id not in CATALOGUE:
            raise ConfigurationError("E_FIGURE_SOURCE_MISSING", f"{figure_id!r} is not a catalogue figure", "/figure_id")
        source = Path(source_table)
        if not source.is_file():
            raise ConfigurationError("E_FIGURE_SOURCE_MISSING", f"figure {figure_id!r} has no source table", "/source")
        digest = sha256_file(source)
        destination = self.directory / SOURCE_DIR / f"{figure_id}.tsv"
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(source.read_bytes())
        record = {"state": "REGISTERED", "reason": None, "source_path": f"{SOURCE_DIR}/{figure_id}.tsv", "source_sha256": digest, "outputs": {}}
        self.figures[figure_id] = record
        return dict(record)

    def require(self, figure_id: str) -> dict[str, Any]:
        """The registered record of a figure; refused before any output is written when it is missing or changed."""
        record = self.figures.get(figure_id)
        if record is None:
            raise ConfigurationError("E_FIGURE_SOURCE_MISSING", f"no registered source table for {figure_id!r}; no output is produced", "/figure_id")
        current = self.directory / record["source_path"]
        if not current.is_file() or sha256_file(current) != record["source_sha256"]:
            raise ConfigurationError("E_FIGURE_SOURCE_MISSING", f"the source table of {figure_id!r} is missing or changed since registration; no output is produced", "/source")
        return dict(record)

    def produce(self, figure_id: str, kind: str, render: Callable[[Path, Path], None]) -> Path:
        """Render one output of a registered figure from its registered source table and record its hashes."""
        if kind not in FIGURE_FORMATS:
            raise ConfigurationError("E_FIGURE_FORMAT", f"unknown figure format {kind!r}", "/formats")
        record = self.require(figure_id)
        output = self.directory / f"{figure_id}.{kind}"
        render(self.directory / record["source_path"], output)
        if not output.is_file():
            raise ConfigurationError("E_FIGURE_OUTPUT_MISSING", f"the {kind} output of {figure_id!r} was not written", "/outputs")
        self.figures[figure_id]["outputs"][kind] = {"path": output.name, "sha256": sha256_file(output), "source_sha256": record["source_sha256"]}
        return output

    def record_missing_output(self, figure_id: str, kind: str, reason: str) -> None:
        """An output that cannot be produced is NOT_RUN with its reason; the other outputs still run (E_FIGURE_OUTPUT_MISSING)."""
        self.figures.setdefault(figure_id, {"state": "NOT_REGISTERED", "reason": None, "source_path": None, "source_sha256": None, "outputs": {}})
        self.figures[figure_id]["outputs"][kind] = {"state": "NOT_RUN", "code": "E_FIGURE_OUTPUT_MISSING", "reason": reason}

    def write(self) -> Path:
        """Write registry.json (style, figures, outputs with SHA-256). UTF-8 with LF line endings."""
        self.directory.mkdir(parents=True, exist_ok=True)
        document = {"schema": "figure-catalogue-registry/1", "style": self.style, "font_substitution": None,
                    "catalogue": {key: {"packet": v[0], "area": v[1], "source_family": v[2], "claim_label": v[3]} for key, v in CATALOGUE.items()},
                    "figures": {key: self.figures[key] for key in sorted(self.figures)}}
        path = self.directory / REGISTRY_FILE
        path.write_text(json.dumps(document, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
        return path
