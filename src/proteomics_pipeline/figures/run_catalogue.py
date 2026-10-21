"""The figure catalogue of a run (contracts/figures.md, Output layout; D-72 retrofit; `proteomics run` with report.figures).

``stage`` is called by the workflow after the report stages. It resolves ``report.figures`` (typed refusals make only this stage
INAPPLICABLE), discovers the registered figure-source tables that the producers of the run wrote (``figure_sources/<stem>.tsv``
and ``heatmap_source.tsv``), applies the retrofit rules to them, writes the three outputs of each figure into the catalogue
directory, and writes ``registry.json`` with the resolved style, the font-substitution record and the retrofit notes.
"""
from __future__ import annotations

import json
import re
import tempfile
from pathlib import Path
from typing import Any, Mapping

from ..errors import ConfigurationError
from ..provenance import is_os_metadata
from ..runtime import run_r_code
from .config import resolve_figure_config
from .catalogue import build_figure
from .registry import FigureRegistry
from . import run_sources
from .retrofit import build_retrofit, rule_for

ASSOCIATION_STEM = "association_heatmap"


def font_substitution() -> str | None:
    """The font status of the R16c theme on this machine: None when Arial is available, else the recorded substitution."""
    code = (
        "suppressPackageStartupMessages(library(proteomicsCore))\n"
        "status <- figure_font_status()\n"
        "cat(jsonlite::toJSON(list(substitution = status$substitution), auto_unbox = TRUE, null = 'null'))\n"
    )
    result = run_r_code(code)
    if result.returncode != 0:
        raise ConfigurationError("E_DEPENDENCY_UNAVAILABLE", f"the font status could not be read from proteomicsCore: {result.stderr[-200:]}")
    return json.loads(result.stdout.strip().splitlines()[-1])["substitution"]


def discover(run_root: Path) -> tuple[list[dict[str, Any]], list[dict[str, str]]]:
    """Every retrofittable registered table of the run, in path order, and the tables that no rule covers."""
    found: list[dict[str, Any]] = []
    skipped: list[dict[str, str]] = []
    used: set[str] = set()
    for path in sorted(run_root.rglob("*.tsv")):
        if is_os_metadata(path):
            continue
        relative = path.relative_to(run_root).as_posix()
        if relative.startswith("report/figure_catalogue/"):
            continue
        if path.parent.name == "figure_sources":
            stem = path.stem
        elif path.name == "heatmap_source.tsv":
            stem = ASSOCIATION_STEM
        else:
            continue
        rule = rule_for(stem)
        if rule is None:
            skipped.append({"path": relative, "reason": "no retrofit rule for this table"})
            continue
        if rule.skip:
            skipped.append({"path": relative, "reason": rule.note or "companion table, drawn with its figure"})
            continue
        companion = None
        if rule.companion is not None:
            name = rule.companion.removeprefix("parent:")
            folder = path.parent.parent if rule.companion.startswith("parent:") else path.parent
            candidate = folder / f"{name}.tsv"
            companion = candidate if candidate.is_file() else None
        figure_id = stem
        if figure_id in used:
            prefix = re.sub(r"[^A-Za-z0-9]+", "_", str(path.parent.parent.relative_to(run_root))).strip("_")
            figure_id = f"{prefix}_{stem}"
        used.add(figure_id)
        found.append({"figure_id": figure_id, "path": path, "rule": rule, "companion": companion, "relative": relative})
    return found, skipped


def build(run_root: Path, figures_block: Mapping[str, Any] | None, catalogue_dir: Path) -> dict[str, Any]:
    """Write the catalogue of the run into ``catalogue_dir``. Typed refusals of the block propagate (nothing is written)."""
    style = resolve_figure_config(figures_block)
    registry = FigureRegistry(catalogue_dir, style)
    registry.font_substitution = font_substitution()
    found, skipped = discover(run_root)
    registry.retrofit = {"producers": {"R10b": "report figures of the volcano and response-axes tables", "R14": "post-DE figure-source tables",
                                       "R13": "not retrofitted: the PERMANOVA producer writes no figure-source table"},
                         "not_retrofitted": skipped}
    results = []
    for item in found:
        try:
            result = build_retrofit(registry, item["figure_id"], item["path"], item["rule"], style["formats"], style["journal_width_mm"],
                                    companion=item["companion"])
        except ConfigurationError as error:
            if item["figure_id"] not in registry.figures:
                registry.figures[item["figure_id"]] = {"state": "NOT_RUN", "reason": error.message, "source_path": None, "source_sha256": None,
                                                       "outputs": {}, "kind": item["rule"].kind, "companions": {}}
            registry.figures[item["figure_id"]]["state"] = "NOT_RUN"
            registry.figures[item["figure_id"]]["reason"] = error.message
            result = {"figure_id": item["figure_id"], "kind": item["rule"].kind, "outputs": {},
                      "refusals": [{"code": error.code, "reason": error.message}]}
        result["source"] = item["relative"]
        results.append(result)
    with tempfile.TemporaryDirectory(prefix="catalogue-sources-") as work:
        jobs, run_skipped = run_sources.build_jobs(run_root, Path(work))
        registry.retrofit["run_sources"] = run_skipped
        for job in jobs:
            try:
                result = build_figure(run_root, job["figure_id"], job["source"], registry=registry, output_id=job["output_id"],
                                      formats=style["formats"], width_mm=style["journal_width_mm"], **job["options"])
            except ConfigurationError as error:
                registry.figures.setdefault(job["output_id"], {"state": "NOT_RUN", "reason": error.message, "source_path": None,
                                                              "source_sha256": None, "outputs": {}, "kind": job["figure_id"], "companions": {}})
                registry.figures[job["output_id"]]["state"] = "NOT_RUN"
                registry.figures[job["output_id"]]["reason"] = error.message
                result = {"figure_id": job["output_id"], "kind": job["figure_id"], "outputs": {},
                          "refusals": [{"code": error.code, "reason": error.message}]}
            result["source"] = job["output_id"]
            results.append(result)
    registry.write()
    return {"figures": results, "not_retrofitted": skipped, "formats": style["formats"]}


def stage(run_root: Path, figures_block: Mapping[str, Any] | None, temp: Path) -> tuple[str, str | None, str]:
    """The workflow stage: (state, reason_code, message). An invalid block is INAPPLICABLE with its typed reason."""
    try:
        summary = build(run_root, figures_block, temp)
    except ConfigurationError as error:
        return "INAPPLICABLE", error.code, error.message
    except Exception as error:                      # an unexpected defect fails this optional stage only
        return "FAILED", "E_FIGURE_OUTPUT_MISSING", f"{type(error).__name__}: {error}"
    produced = sum(len(r["outputs"]) for r in summary["figures"])
    missing = sum(len(r["refusals"]) for r in summary["figures"])
    return "COMPLETED", None, f"{len(summary['figures'])} figure(s) registered; {produced} output(s) written; {missing} output(s) not run"
