"""R10b acceptance: V095-V100 (complete views, pathway/response views, figures, methods, robustness, documented workflow).
Also re-exercises the R10a interface read-only (stub report still produced and unchanged in kind)."""
from __future__ import annotations

import html.parser
import importlib.util
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from proteomics_pipeline import workflow
from proteomics_pipeline.reporting import assemble as assembler
from proteomics_pipeline.reporting import full

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("design_builders", ROOT / "tests" / "fixtures" / "design" / "builders.py")
B = importlib.util.module_from_spec(spec); spec.loader.exec_module(B)
pytestmark = pytest.mark.skipif(B.rscript() is None, reason="NOT_RUN: Rscript/proteomicsCore/limma unavailable")
DEMO = ROOT / "tests" / "fixtures" / "reports" / "demo-phase2.json"


class Rows(html.parser.HTMLParser):
    def __init__(self):
        super().__init__(); self.captions = []; self.counts = []; self._in_caption = False; self._rows = 0; self.urls = []
    def handle_starttag(self, tag, attrs):
        if tag == "table": self._rows = 0
        if tag == "caption": self._in_caption = True; self.captions.append("")
        if tag == "tr": self._rows += 1
        for k, v in attrs:
            if k in ("href", "src"): self.urls.append(v)
    def handle_endtag(self, tag):
        if tag == "caption": self._in_caption = False
        if tag == "table": self.counts.append(self._rows - 1)
    def handle_data(self, data):
        if self._in_caption: self.captions[-1] += data


def run_example(tmp_path, name="example-independent.json", mutate=None, formats=("png", "pdf", "svg")):
    raw = json.loads((B.EXAMPLES / name).read_text(encoding="utf-8"))
    for key in ("matrix", "observations", "features", "source_provenance"):
        raw["input"][key] = str((B.EXAMPLES / raw["input"][key]).resolve())
    raw["report"]["figure_formats"] = list(formats)
    if mutate:
        mutate(raw)
    tmp_path.mkdir(parents=True, exist_ok=True)
    path = tmp_path / "analysis.json"; path.write_text(json.dumps(raw), encoding="utf-8")
    payload, code = workflow.run_command(path, tmp_path / "run")
    return payload, code, tmp_path / "run"


def test_v095_every_planned_endpoint_rendered(tmp_path):
    payload, code, run = run_example(tmp_path, mutate=lambda c: c["models"][0]["coverage"].update({"minimum_observed_per_group": 4, "minimum_fraction": 1.0}))
    assert code == 0, payload
    rows = B.read_tsv(run / "dea" / "zero_null.tsv")
    assert {r["contrast_id"] for r in rows} == {"disease-control", "treated-disease", "treated-control"}
    assert any(r["eligibility"] != "tested" for r in rows) and any(r["role"] == "secondary" for r in rows)
    page = (run / "report-full" / "index.html").read_text(encoding="utf-8")
    parser = Rows(); parser.feed(page)
    i = next(i for i, c in enumerate(parser.captions) if c.startswith("Complete zero_null table"))
    assert parser.counts[i] == len(rows) == 24
    for r in rows:
        assert f"<td>{r['feature_id']}</td>" in page
    families = B.read_tsv(run / "dea" / "families.tsv")
    data = json.loads((run / "report-full" / "report_data.json").read_text(encoding="utf-8"))
    assert [f["n_planned"] for f in data["dea"]["families"]] == [f["n_planned"] for f in families]   # totals, never display-filtered counts


def test_v096_pathway_response_and_multivariate_views(tmp_path):
    payload, code = workflow.run_command(DEMO, tmp_path / "run")
    assert code == 0, payload
    data = json.loads((tmp_path / "run" / "report-full" / "report_data.json").read_text(encoding="utf-8"))
    response = data["sections"]["response"]
    desc = next(t for t in response["tables"] if t["caption"].startswith("Descriptive response"))
    assert "P" not in [c.upper() for c in desc["columns"]] and not any("p_value" in c for c in desc["columns"])
    assert any("not statistical equivalence" in n for n in response["notes"])
    assert data["sections"]["permanova"]["state"] == "COMPLETED"
    page = (tmp_path / "run" / "report-full" / "index.html").read_text(encoding="utf-8").lower()
    for phrase in ("confirms the", "validated mechanism", "rescues", "proves"):
        assert phrase not in page


def test_v096_pathway_section_labels(tmp_path):
    from proteomics_pipeline.reporting.pathway_data import NULL_LABELS
    assert "not sample-level inference" in NULL_LABELS["preranked_gene_set"] and "exploratory" in NULL_LABELS["ora_up"]


def test_v097_figures_match_their_source_tables(tmp_path):
    payload, code, run = run_example(tmp_path)
    figs = json.loads((run / "report-full" / "figures.json").read_text(encoding="utf-8"))
    assert figs and all(f["source"] for f in figs)
    rows = {(r["contrast_id"], r["feature_id"]): r for r in B.read_tsv(run / "dea" / "zero_null.tsv")}
    for fig in figs:
        src = B.read_tsv(run / "report-full" / fig["source"])
        contrast = fig["stem"].split("_limma-main_")[1]
        import math
        for r in src:
            d = rows[(contrast, r["feature_id"])]
            assert float(r["effect_log2"]) == float(d["effect"]) and abs(float(r["neg_log10_p"]) + math.log10(float(d["p_value"]))) <= 1e-12
        svg = (run / "report-full" / f"figures/{fig['stem']}.svg").read_text(encoding="utf-8")
        assert svg.count("<circle") == len({(r["effect_log2"], r["neg_log10_p"]) for r in src}) and "*" not in svg   # identical coordinates share one mark (item 10)
        for f in fig["files"]:
            head = (run / "report-full" / f).read_bytes()[:8]
            assert head.startswith(b"%PDF") if f.endswith(".pdf") else (head == b"\x89PNG\r\n\x1a\n" if f.endswith(".png") else head.startswith(b"<svg"))


def test_v098_methods_follow_the_executed_run(tmp_path):
    _, _, zero = run_example(tmp_path / "a")
    _, _, treat = run_example(tmp_path / "b", name="example-effect-threshold.json")
    m1 = (zero / "report-full" / "methods_full.md").read_text(encoding="utf-8"); m2 = (treat / "report-full" / "methods_full.md").read_text(encoding="utf-8")
    assert "protein_zero_null" in m1 and "protein_treat" in m2 and m1 != m2
    for text in (m1, m2):
        for engine in ("DEqMS", "proDA", "fgsea", "CAMERA", "PERMANOVA"):
            assert engine not in text                                          # not executed, not mentioned
        assert "Unknown: study.tissue" in text


def test_v099_null_partial_failed_and_cancelled_reports_are_truthful(tmp_path):
    payload, code, qc = run_example(tmp_path / "qc", mutate=lambda c: c["runtime"].update({"scope": "qc_only"}))
    page = (qc / "report-full" / "index.html").read_text(encoding="utf-8")
    assert code == 0 and "no differential rows exist for this run" in page
    status = json.loads((qc / "run_status.json").read_text(encoding="utf-8"))
    snapshot = {k: status[k] for k in ("run_id", "plan_hash", "exit_code", "reason_code", "requested_phase", "implemented_capabilities")}
    snapshot.update(state="CANCELLED", reason_code="E_CHILD_TIMEOUT", exit_code=6,
                    stages=[dict(s, message="") for s in status["stages"]] + [{"stage_id": "limma", "capability": "limma", "required": True, "state": "CANCELLED", "reason_code": "E_CHILD_TIMEOUT", "result_path": None, "message": "interrupted"}])
    out = tmp_path / "cancelled"; out.mkdir()
    full.build_report(qc, snapshot, out, [])
    text = (out / "index.html").read_text(encoding="utf-8")
    assert "state-CANCELLED" in text and "COMPLETED — exit code 0" not in text
    values = {f"F{i}": [1e300 if j % 2 == 0 else -1e300 for j in range(12)] for i in range(4)}
    obs = [{"observation_id": f"{g}{i}", "group": g} for g in ("C", "U", "T") for i in range(4)]
    files = B.dataset(tmp_path / "crash" / "data", values, obs)
    path = B.config(tmp_path / "crash" / "data", files, groups=["C", "U", "T"], contrasts=[B.contrast("d", "U", "C"), B.contrast("t", "T", "U")])
    payload, code = workflow.run_command(path, tmp_path / "crash" / "run")
    crash = (tmp_path / "crash" / "run" / "report-full" / "index.html").read_text(encoding="utf-8")
    assert code == 4 and "state-FAILED" in crash and "E_ENGINE_FAILED" in crash


def _documented_block(guide: str, fence: str) -> list[str]:
    block = guide.split("## Workflow", 1)[1].split(f"```{fence}", 1)[1].split("```", 1)[0]
    return [line.strip() for line in block.splitlines() if line.strip()]


def test_v100_documented_commands_run_verbatim(tmp_path):
    import sysconfig
    guide = (ROOT / "docs" / "user-guide" / "usage.md").read_text(encoding="utf-8")
    bash = [line for line in _documented_block(guide, "bash") if line.startswith(".venv/bin/proteomics")]
    # review follow-up 2026-10-03: the PowerShell variant is documented and checked token-for-token against the bash one
    powershell = _documented_block(guide, "powershell")
    assert powershell[0] == '$env:R_LIBS_USER = "$PWD\\.r-lib"'
    ps_commands = [line for line in powershell[1:]]
    assert len(bash) == len(ps_commands) == 7
    for posix, windows in zip(bash, ps_commands):
        p_args, w_args = posix.split(), windows.split()          # no quoting in either block; PowerShell is not parsed with POSIX shlex
        assert p_args[0] == ".venv/bin/proteomics" and w_args[0] == ".venv\\Scripts\\proteomics.exe" and p_args[1:] == w_args[1:]
    env = {**os.environ}
    help_text = subprocess.run([sys.executable, "-m", "proteomics_pipeline", "--help"], capture_output=True, text=True, encoding="utf-8").stdout
    # The documented executable is the installed console script; it lives in sysconfig's scripts directory
    # (<venv>/bin on POSIX, <venv>\\Scripts on Windows), so the block for this platform runs through it.
    console = shutil.which("proteomics", path=os.pathsep.join([sysconfig.get_path("scripts"), str(Path(sys.executable).parent)]))
    executable = [console] if console else [sys.executable, "-m", "proteomics_pipeline"]
    commands = ps_commands if os.name == "nt" else bash
    for command in commands:
        # CI run 36982402784: a Windows tmp path has backslashes, which POSIX shlex treats as escapes; substitute a forward-slash path after splitting
        argv = [arg.replace("runs/", f"{tmp_path.as_posix()}/runs/") for arg in command.split()]
        assert argv[1] in help_text                                            # only implemented commands are documented
        result = subprocess.run([*executable, *argv[1:]], cwd=ROOT, capture_output=True, text=True, env=env, encoding="utf-8")
        assert result.returncode == 0, (command, result.stdout[-500:], result.stderr[-500:])
    assert (tmp_path / "runs" / "comparison" / "index.html").is_file() and (tmp_path / "runs" / "example" / "report" / "index.html").is_file()
