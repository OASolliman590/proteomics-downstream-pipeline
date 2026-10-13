"""R16a acceptance: V195 (report.figures configuration) and V197 (figure registry), FR-195 and FR-197.

Oracles are written here: the configuration defaults and the figure IDs are copied from
specs/017-visualization/contracts/figures.md, and SHA-256 values are computed with hashlib, not with the code under test.
Negative cases use the same loader and registry with one changed input.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from proteomics_pipeline.config import load_config
from proteomics_pipeline.errors import ConfigurationError
from proteomics_pipeline.figures import FigureRegistry, resolve_figure_config

ROOT = Path(__file__).resolve().parents[2]
DEMO_CONFIG = ROOT / "tests" / "fixtures" / "reports" / "demo-phase2.json"

# Oracle: the configuration table of contracts/figures.md (defaults, written out here).
CONTRACT_DEFAULTS = {
    "formats": ["png", "pzfx", "html"],
    "dpi": 300,
    "journal_width_mm": 120,
    "palette": "okabe_ito",
    "group_colours": None,
    "font": "Arial",
    "p_annotation": "exact",
}
# Oracle: the registered catalogue figure IDs used by this test (contracts/figures.md, "Figure IDs").
REGISTERED_IDS = ["qc_intensity_distributions", "diff_volcano_labelled", "bio_roc_ci_band"]
UNREGISTERED_ID = "diff_ma_plot"


def _config_with_figures(tmp_path: Path, figures: dict | None) -> Path:
    """A copy of the demo report configuration; ``figures`` is written as report.figures when it is not None."""
    raw = json.loads(DEMO_CONFIG.read_text(encoding="utf-8"))
    if figures is not None:
        raw["report"]["figures"] = figures
    path = tmp_path / "config.json"
    path.write_text(json.dumps(raw, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    return path


def _resolved_from_file(path: Path) -> dict:
    report = load_config(path)["report"]
    return resolve_figure_config(report.get("figures"))


def test_v195_defaults_without_block_and_with_explicit_defaults(tmp_path):
    # (a) no report.figures block
    absent = _resolved_from_file(_write(tmp_path / "a", None))
    # (b) explicit defaults
    explicit = _resolved_from_file(_write(tmp_path / "b", {"formats": ["png", "pzfx", "html"], "dpi": 300, "journal_width_mm": 120,
                                                          "palette": "okabe_ito", "font": "Arial", "p_annotation": "exact"}))
    assert absent == CONTRACT_DEFAULTS
    assert explicit == CONTRACT_DEFAULTS
    # the resolved values are written to the run settings (the catalogue registry of the run)
    run_dir = tmp_path / "run" / "report" / "figure_catalogue"
    FigureRegistry(run_dir, absent).write()
    settings = json.loads((run_dir / "registry.json").read_text(encoding="utf-8"))["style"]
    assert settings == CONTRACT_DEFAULTS


def test_v195_report_figure_formats_is_not_merged_into_figures(tmp_path):
    raw = json.loads(DEMO_CONFIG.read_text(encoding="utf-8"))
    assert raw["report"]["figure_formats"] == ["png", "pdf"]
    # the legacy key keeps its meaning and never changes the figures defaults
    assert _resolved_from_file(_write(tmp_path, None))["formats"] == CONTRACT_DEFAULTS["formats"]


def test_v195_negative_dpi_150_fails_e_figure_dpi(tmp_path):
    path = _write(tmp_path, {"dpi": 150})
    with pytest.raises(ConfigurationError) as caught:
        _resolved_from_file(path)
    assert caught.value.code == "E_FIGURE_DPI"
    assert not (tmp_path / "run").exists(), "no figure or registry is written on refusal"


def test_v195_negative_jpg_format_fails_e_figure_format(tmp_path):
    path = _write(tmp_path, {"formats": ["jpg"]})
    with pytest.raises(ConfigurationError) as caught:
        _resolved_from_file(path)
    assert caught.value.code == "E_FIGURE_FORMAT"
    assert not (tmp_path / "run").exists()


def test_v195_negative_unknown_key_fails_with_the_unknown_key_error(tmp_path):
    path = _write(tmp_path, {"dpi": 300, "colour_scheme": "dark"})
    with pytest.raises(ConfigurationError) as caught:       # the configuration loader refuses it before any figure logic
        load_config(path)
    assert caught.value.code == "E_CONFIG_SCHEMA"
    with pytest.raises(ConfigurationError) as direct:       # and the figure resolver refuses it as well
        resolve_figure_config({"colour_scheme": "dark"})
    assert direct.value.code == "E_CONFIG_SCHEMA"
    assert not (tmp_path / "run").exists()


@pytest.mark.parametrize(("block", "code"), [
    ({"formats": []}, "E_FIGURE_FORMAT"),
    ({"formats": ["png", "png"]}, "E_FIGURE_FORMAT"),
    ({"journal_width_mm": 100}, "E_CONFIG_SCHEMA"),
    ({"font": "Helvetica"}, "E_CONFIG_SCHEMA"),
    ({"p_annotation": "letters"}, "E_CONFIG_SCHEMA"),
    ({"palette": ["#GG0000"]}, "E_FIGURE_PALETTE"),
    ({"palette": ["#0072B2", "#0072B2"]}, "E_FIGURE_PALETTE"),
    ({"group_colours": {"Control": "blue"}}, "E_FIGURE_PALETTE"),
])
def test_v195_other_refusals_are_typed(block, code):
    with pytest.raises(ConfigurationError) as caught:
        resolve_figure_config(block)
    assert caught.value.code == code


def _write(directory: Path, figures: dict | None) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    return _config_with_figures(directory, figures)


def _synthetic_table(path: Path, seed: int) -> Path:
    """A small synthetic figure-source table (deterministic, no study data)."""
    rows = ["feature\tvalue"] + [f"F{i:03d}\t{((i * 37 + seed) % 101) / 10:.1f}" for i in range(1, 21)]
    path.write_text("\n".join(rows) + "\n", encoding="utf-8", newline="\n")
    return path


def test_v197_registry_lists_source_path_and_sha256_and_generates_outputs_from_it(tmp_path):
    inputs = [_synthetic_table(tmp_path / f"input_{i}.tsv", seed=i) for i in range(3)]
    oracle = {fid: hashlib.sha256(path.read_bytes()).hexdigest() for fid, path in zip(REGISTERED_IDS, inputs)}
    run_dir = tmp_path / "run" / "report" / "figure_catalogue"
    registry = FigureRegistry(run_dir, CONTRACT_DEFAULTS)
    for fid, path in zip(REGISTERED_IDS, inputs):
        registry.register(fid, path)

    def render(source: Path, output: Path) -> None:
        output.write_bytes(source.read_bytes())            # every output is rendered from the registered source only

    for fid in REGISTERED_IDS:
        for kind in CONTRACT_DEFAULTS["formats"]:
            registry.produce(fid, kind, render)
    registry.write()
    document = json.loads((run_dir / "registry.json").read_text(encoding="utf-8"))
    for fid in REGISTERED_IDS:
        record = document["figures"][fid]
        assert record["source_path"] == f"sources/{fid}.tsv"
        assert record["source_sha256"] == oracle[fid]
        for kind in CONTRACT_DEFAULTS["formats"]:
            assert record["outputs"][kind]["source_sha256"] == oracle[fid]
            assert record["outputs"][kind]["sha256"] == oracle[fid]    # the output is a copy of its registered source


def test_v197_negative_request_without_source_table_fails_and_writes_nothing(tmp_path):
    inputs = [_synthetic_table(tmp_path / f"input_{i}.tsv", seed=i) for i in range(3)]
    run_dir = tmp_path / "run" / "report" / "figure_catalogue"
    registry = FigureRegistry(run_dir, CONTRACT_DEFAULTS)
    for fid, path in zip(REGISTERED_IDS, inputs):
        registry.register(fid, path)
    calls = []

    def render(source: Path, output: Path) -> None:
        calls.append(output)
        output.write_bytes(source.read_bytes())

    before = sorted(p.name for p in run_dir.iterdir())
    with pytest.raises(ConfigurationError) as caught:
        registry.produce(UNREGISTERED_ID, "png", render)
    assert caught.value.code == "E_FIGURE_SOURCE_MISSING"
    assert calls == [], "the renderer is never reached without a registered source"
    assert sorted(p.name for p in run_dir.iterdir()) == before
    assert not (run_dir / f"{UNREGISTERED_ID}.png").exists()


def test_v197_negative_missing_or_changed_source_is_refused(tmp_path):
    source = _synthetic_table(tmp_path / "input.tsv", seed=7)
    run_dir = tmp_path / "run" / "report" / "figure_catalogue"
    registry = FigureRegistry(run_dir, CONTRACT_DEFAULTS)
    with pytest.raises(ConfigurationError) as absent:
        registry.register("qc_pca_ellipses", tmp_path / "does-not-exist.tsv")
    assert absent.value.code == "E_FIGURE_SOURCE_MISSING"
    registry.register("qc_pca_ellipses", source)
    (run_dir / "sources" / "qc_pca_ellipses.tsv").write_text("feature\tvalue\nF001\t9.9\n", encoding="utf-8", newline="\n")
    with pytest.raises(ConfigurationError) as changed:
        registry.produce("qc_pca_ellipses", "html", lambda s, o: o.write_bytes(s.read_bytes()))
    assert changed.value.code == "E_FIGURE_SOURCE_MISSING"
    assert not (run_dir / "qc_pca_ellipses.html").exists()
