"""R10a acceptance: V091-V094 (typed ReportData, one-run integration, offline HTML, QC/inclusion sections)."""
from __future__ import annotations

import html.parser
import importlib.util
import json
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pytest

from proteomics_pipeline import workflow
from proteomics_pipeline.reporting import assemble as assembler

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("design_builders", ROOT / "tests" / "fixtures" / "design" / "builders.py")
B = importlib.util.module_from_spec(spec); spec.loader.exec_module(B)
pytestmark = pytest.mark.skipif(B.rscript() is None, reason="NOT_RUN: Rscript/proteomicsCore/limma unavailable")
OBS = [{"observation_id": f"{g}{i}", "group": g} for g in ("C", "U", "T") for i in range(4)]
CONTRASTS = [B.contrast("d", "U", "C"), B.contrast("t", "T", "U"), B.contrast("r", "T", "C", role="secondary")]


def zero_discovery(tmp_path, mutate=None, values=None):
    offsets = [-0.3, -0.1, 0.1, 0.3]
    values = values or {f"F{f + 1}": [10 + f + offsets[(i + f) % 4] for i in range(12)] for f in range(6)}
    files = B.dataset(tmp_path / "data", values, OBS)
    return B.config(tmp_path / "data", files, groups=["C", "U", "T"], contrasts=[dict(c) for c in CONTRASTS], mutate=mutate)


def failing(tmp_path):
    values = {f"F{f + 1}": [1e300 if i % 2 == 0 else -1e300 for i in range(12)] for f in range(4)}
    return zero_discovery(tmp_path, values=values)


def run_cli(*args, cwd=None):
    return subprocess.run([sys.executable, "-m", "proteomics_pipeline", *args], capture_output=True, text=True, cwd=cwd, env=os.environ.copy())


class Parser(html.parser.HTMLParser):
    def __init__(self):
        super().__init__(); self.tags = []; self.urls = []; self.captions = 0; self.th_scope = 0; self.th = 0; self.tables = 0
    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs); self.tags.append((tag, attrs))
        for key in ("href", "src"):
            if key in attrs:
                self.urls.append(attrs[key])
        self.captions += tag == "caption"; self.tables += tag == "table"; self.th += tag == "th"; self.th_scope += tag == "th" and "scope" in attrs


# ----------------------------------------------------------------------------- V091
def test_v091_zero_discoveries_reported_as_zero_from_actual_tables(tmp_path):
    payload, code = workflow.run_command(zero_discovery(tmp_path), tmp_path / "run")
    assert code == 0 and payload["state"] == "COMPLETED"
    data = json.loads((tmp_path / "run" / "report" / "report_data.json").read_text())
    assert data["dea"]["state"] == "COMPLETED"
    rows = B.read_tsv(tmp_path / "run" / "dea" / "zero_null.tsv")
    summary = {(c["contrast_id"], c["hypothesis_type"]): c for c in data["dea"]["contrast_summaries"]}
    for contrast in ("d", "t", "r"):
        actual = [r for r in rows if r["contrast_id"] == contrast]
        assert summary[(contrast, "protein_zero_null")]["n_planned"] == len(actual) == 6
        assert summary[(contrast, "protein_zero_null")]["n_q_at_or_below_cutoff"] == 0
    assert [f["rejection_count"] for f in data["dea"]["families"]] == ["0", "0"]
    assert all(s["verified"] for s in data["stages"])
    page = (tmp_path / "run" / "report" / "index.html").read_text()
    assert "PASS" not in page
    lowered = page.lower()
    import re
    for match in re.finditer("no biological effect", lowered):                     # only ever negated
        assert "not evidence of" in lowered[max(0, match.start() - 40):match.start()]


def test_v091_missing_dea_stays_unknown_not_zero(tmp_path):
    path = zero_discovery(tmp_path, mutate=lambda c: c["runtime"].update({"scope": "qc_only"}))
    payload, code = workflow.run_command(path, tmp_path / "run")
    assert code == 0
    data = json.loads((tmp_path / "run" / "report" / "report_data.json").read_text())
    assert data["dea"]["state"] == "NOT_REQUESTED" and data["dea"]["values"] is None and "families" not in data["dea"]
    assert data["design"]["state"] == "NOT_REQUESTED"


def test_v091_failed_required_stage_remains_failed(tmp_path):
    payload, code = workflow.run_command(failing(tmp_path), tmp_path / "run")
    assert payload["state"] == "FAILED" and code == 4
    limma = next(s for s in payload["stages"] if s["stage_id"] == "limma")
    assert limma["state"] == "FAILED" and limma["reason_code"] == "E_ENGINE_FAILED"
    data = json.loads((tmp_path / "run" / "report" / "report_data.json").read_text())
    assert data["run"]["state"] == "FAILED" and data["dea"]["state"] == "FAILED" and data["dea"]["values"] is None
    assert (tmp_path / "run" / "logs" / "limma-failed").is_dir() and not (tmp_path / "run" / "dea").exists()


def test_v091_negative_unverified_artifact_is_not_completed(tmp_path):
    workflow.run_command(zero_discovery(tmp_path), tmp_path / "run")
    status = json.loads((tmp_path / "run" / "run_status.json").read_text())
    (tmp_path / "run" / "dea" / "families.tsv").write_text("family_id\trejection_count\nforged\t99\n")
    snapshot = {"run_id": status["run_id"], "plan_hash": status["plan_hash"], "state": status["state"], "exit_code": status["exit_code"], "reason_code": status["reason_code"],
                "requested_phase": 1, "implemented_capabilities": status["implemented_capabilities"], "stages": [dict(s, message="") for s in status["stages"]]}
    data = assembler.assemble(tmp_path / "run", snapshot)
    limma = next(s for s in data["stages"] if s["stage_id"] == "limma")
    assert limma["display_state"] == "UNVERIFIED" and limma["verified"] is False
    assert data["dea"]["state"] == "UNVERIFIED" and data["dea"]["values"] is None
    assert run_cli("verify", "--run", str(tmp_path / "run"), "--json").returncode == 5


# ----------------------------------------------------------------------------- V092
def test_v092_one_command_runs_the_phase1_dag_with_plan_before_fit(tmp_path):
    result = run_cli("run", "--config", str(zero_discovery(tmp_path)), "--output", str(tmp_path / "run"), "--json")
    assert result.returncode == 0, result.stdout + result.stderr
    payload = json.loads(result.stdout)
    # Amendment A-2026-10-01-06: R10b adds an optional full-report stage after the thin report; the Phase 1 DAG order is unchanged.
    assert [s["stage_id"] for s in payload["stages"]][:5] == ["intake", "preprocessing", "design", "limma", "report"]
    assert [s["stage_id"] for s in payload["stages"]][5:] in ([], ["report_full"])
    run = tmp_path / "run"
    limma = json.loads((run / "dea" / "stage-result.json").read_text())
    started = datetime.fromisoformat(limma["started_at"].replace("Z", "+00:00")).timestamp()
    assert (run / "plan.json").stat().st_mtime <= started + 1.0
    assert limma["plan_hash"] == json.loads((run / "plan.json").read_text())["plan_hash"] == payload["plan_hash"]
    assert limma["exit_code"] == 0 and limma["session_info_path"] == "sessionInfo.txt"


def test_v092_genuine_r_failure_stops_dependants_with_partial_report(tmp_path):
    result = run_cli("run", "--config", str(failing(tmp_path)), "--output", str(tmp_path / "run"), "--json")
    assert result.returncode == 4
    payload = json.loads(result.stdout)
    assert payload["state"] == "FAILED" and payload["report"] == "report/index.html"
    page = (tmp_path / "run" / "report" / "index.html").read_text()
    assert "state-FAILED" in page and "E_ENGINE_FAILED" in page


def test_v092_negative_phase1_score_request_is_rejected(tmp_path):
    path = zero_discovery(tmp_path, mutate=lambda c: c.update({"score_test": "independent_exact"}))
    result = run_cli("run", "--config", str(path), "--output", str(tmp_path / "run"), "--json")
    assert result.returncode == 2 and json.loads(result.stdout)["errors"][0]["code"] == "E_PHASE_CAPABILITY"
    assert not (tmp_path / "run").exists() or not any((tmp_path / "run").rglob("*score*"))


# ----------------------------------------------------------------------------- V093
def test_v093_offline_html_without_scripts_or_remote_resources(tmp_path):
    workflow.run_command(zero_discovery(tmp_path), tmp_path / "run")
    report = tmp_path / "run" / "report" / "index.html"
    text = report.read_text()
    parser = Parser(); parser.feed(text)
    tags = [t for t, _ in parser.tags]
    assert "script" not in tags and "iframe" not in tags and "style" in tags
    assert not any(t == "link" and a.get("rel") == "stylesheet" for t, a in parser.tags)
    for url in parser.urls:
        assert not url.startswith(("http:", "https:", "//", "data:", "file:")), url
        if url.startswith("#"):
            assert f'id="{url[1:]}"' in text
            continue
        assert (report.parent / url).resolve().is_file(), url
    assert parser.tables >= 6 and parser.captions == parser.tables and parser.th == parser.th_scope > 0
    assert str(tmp_path) not in text and "/private/" not in text and "/Users/" not in text
    assert '<nav aria-label="Report sections">' in text


# ----------------------------------------------------------------------------- V094
def test_v094_qc_inclusion_sections_show_biological_vs_technical_n(tmp_path):
    obs = []
    values = {f"F{f}": [] for f in range(1, 6)}
    for g in ("C", "U"):
        for unit in range(1, 4):
            injections = 4 if (g == "C" and unit == 1) else 1
            for inj in range(1, injections + 1):
                oid = f"{g}{unit}" + (f"_inj{inj}" if injections > 1 else "")
                obs.append({"observation_id": oid, "biological_unit_id": f"{g}{unit}", "technical_replicate_id": f"i{inj}" if injections > 1 else "NA", "group": g})
                for f in range(1, 6):
                    values[f"F{f}"].append(9 + f + (0.4 if g == "U" else 0) + 0.07 * ((unit * 3 + f + inj) % 4))
    files = B.dataset(tmp_path / "data", values, obs)
    reason = "prespecified: sample flagged by the synthetic study protocol"
    def mutate(c):
        c["preprocessing"]["technical_replicates"] = {"mode": "aggregate", "method": "mean_log2", "biological_unit_column": "biological_unit_id"}
        c["preprocessing"]["exclusions"] = [{"observation_id": "U3", "reason": reason}]
        c["input"]["normalization_state"] = "unknown"; c["study"]["tissue"] = None
    path = B.config(tmp_path / "data", files, groups=["C", "U"], contrasts=[B.contrast("U-C", "U", "C")], mutate=mutate)
    payload, code = workflow.run_command(path, tmp_path / "run")
    assert code == 0, payload
    data = json.loads((tmp_path / "run" / "report" / "report_data.json").read_text())
    groups = {g["group"]: g for g in data["inputs"]["values"]["groups"]}
    assert groups["C"]["n_injections"] == 6 and groups["C"]["n_biological_units"] == 3        # 4 injections -> 1 specimen
    qc_n = {r["group"]: r for r in data["qc"]["tables"]["sample_n"]["rows"]}
    assert qc_n["C"]["n_biological_units"] == "3" and qc_n["C"]["n_injections"] == "6" and qc_n["U"]["n_biological_units"] == "2"
    exclusions = data["qc"]["tables"]["exclusions"]["rows"]
    assert exclusions == [{"observation_id": "U3", "reason": reason}]
    analysed = {r["observation_id"] for r in data["qc"]["tables"]["distributions"]["rows"]}
    assert analysed | {"U3"} == set(data["inputs"]["canonical_observation_ids"]) and "U3" not in analysed   # nobody silently disappears
    gaps = {g["field"] for g in data["inputs"]["values"]["provenance_gaps"]}
    assert {"study.tissue", "input.normalization_state"} <= gaps
    page = (tmp_path / "run" / "report" / "index.html").read_text()
    assert reason in page and "Biological versus technical n by group" in page and "flag only; no sample removed" in page
    for source in (data["qc"]["tables"]["sample_n"]["source"], data["qc"]["tables"]["exclusions"]["source"], data["inputs"]["sources"]["sample_n"]):
        assert (tmp_path / "run" / source).is_file() and f'href="../{source}"' in page
