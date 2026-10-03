import json, os, shutil, subprocess, sys, sysconfig
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[3]; EXAMPLES=ROOT/"configs"/"examples"
def run_cli(*args):
    return subprocess.run([sys.executable,"-m","proteomics_pipeline",*args],cwd=ROOT,text=True,capture_output=True,env={**os.environ,"PYTHONPATH":str(ROOT/"src")}, encoding="utf-8")
def test_help_and_version_are_available_without_r():
    h=run_cli("--help"); v=run_cli("--version")
    assert h.returncode==0 and v.returncode==0; assert "doctor" in h.stdout and "validate" in h.stdout; assert v.stdout.strip()=="0.1.0.dev0"
@pytest.mark.parametrize("name",["example-independent.json","example-paired.json","example-effect-threshold.json"])
def test_schema_only_examples_resolve(name):
    result=run_cli("validate","--config",str(EXAMPLES/name),"--schema-only","--json")
    assert result.returncode==0,result.stderr; payload=json.loads(result.stdout); assert payload["valid"] is True; assert payload["config"]["primary_engine"]=="limma"; assert payload["config"]["response_mode"]=="descriptive_only"; assert payload["config"]["score_test"]=="off"
def test_unknown_fields_have_json_pointer_and_typed_error(tmp_path):
    source=json.loads((EXAMPLES/"example-independent.json").read_text(encoding="utf-8")); source["study"]["typo"]=True; path=tmp_path/"bad.json"; path.write_text(json.dumps(source), encoding="utf-8")
    result=run_cli("validate","--config",str(path),"--schema-only","--json"); assert result.returncode==2; error=json.loads(result.stdout)["errors"][0]; assert error["code"]=="E_CONFIG_SCHEMA"; assert error["pointer"]=="/study"
def test_full_validation_and_plan_are_explicitly_unavailable(tmp_path):
    # Amendment A-2026-10-01-03: R02-R05 now provide real validate/plan/run handlers, so the
    # still-unimplemented case is a required capability without a handler (a Phase 2 DEqMS model, packet R06).
    source=json.loads((EXAMPLES/"example-independent.json").read_text(encoding="utf-8")); source["runtime"]["phase"]=2
    deqms=dict(source["models"][0]); deqms.update({"id":"deqms-required","engine":"deqms","role":"sensitivity","execution_requirement":"required","count_evidence":"counts.tsv","count_aggregation":"protein_level","count_zero_policy":{"mode":"reject"}}); source["models"].append(deqms)
    for key in ("matrix","observations","features","source_provenance"): source["input"][key]=str((EXAMPLES/source["input"][key]).resolve())
    path=tmp_path/"phase2.json"; path.write_text(json.dumps(source), encoding="utf-8"); config=str(path)
    for command in ("validate","plan","run"):
        args=[command,"--config",config,"--json"]
        if command=="validate":args.append("--no-schema-only")
        if command=="plan":args += ["--output",str(tmp_path/"plan.json")]
        if command=="run":args += ["--output",str(tmp_path/"run")]
        # Amendment A-2026-10-01-09: R06 now implements assay_engines, so the handler's absence is simulated in the child process.
        result=subprocess.run([sys.executable,"-c","import sys; sys.modules['proteomics_pipeline.assay_service']=None; from proteomics_pipeline.cli import main; raise SystemExit(main(sys.argv[1:]))",*args],cwd=ROOT,text=True,capture_output=True,env={**os.environ,"PYTHONPATH":str(ROOT/"src")}, encoding="utf-8")
        assert result.returncode==3,(result.stdout,result.stderr); assert "E_CAPABILITY_NOT_IMPLEMENTED" in result.stdout
def test_unknown_command_is_usage_error(): assert run_cli("not-a-command").returncode==2
def test_json_output_is_utf8_under_non_utf8_console(tmp_path):
    # Amendment A-2026-10-01-12: an ASCII/cp1252 console must not crash the CLI on non-ASCII payload text.
    source=json.loads((EXAMPLES/"example-independent.json").read_text(encoding="utf-8")); source["study"]["title"]="Ω-éß study"
    for key in ("matrix","observations","features","source_provenance"): source["input"][key]=str((EXAMPLES/source["input"][key]).resolve())
    folder=tmp_path/"non ascii é Ω"; folder.mkdir(); path=folder/"config.json"; path.write_text(json.dumps(source,ensure_ascii=False), encoding="utf-8")
    env={k:v for k,v in os.environ.items() if not k.startswith(("LC_","LANG","PYTHONIOENCODING","PYTHONUTF8"))}
    env.update({"PYTHONPATH":str(ROOT/"src"),"LC_ALL":"C","PYTHONUTF8":"0"})
    result=subprocess.run([sys.executable,"-m","proteomics_pipeline","validate","--config",str(path),"--schema-only","--json"],cwd=ROOT,capture_output=True,env=env)
    assert result.returncode==0,result.stderr.decode("utf-8","replace")
    assert json.loads(result.stdout.decode("utf-8"))["config"]["study"]["title"]=="Ω-éß study"
def test_installed_console_and_python_module_work_without_source_path(tmp_path):
    env=os.environ.copy(); env.pop("PYTHONPATH",None); python=Path(sys.executable)
    # Amendment A-2026-10-01-11 (CI run 36982402784): console scripts live in sysconfig's "scripts" directory (<prefix>\Scripts on Windows), not always next to python
    console_path=shutil.which("proteomics",path=os.pathsep.join([sysconfig.get_path("scripts"),str(python.parent)])); assert console_path
    empty=tmp_path/"empty"; empty.mkdir(); env["PATH"]=""; config=str(EXAMPLES/"example-independent.json"); verify=str(empty/"verify-run"); Path(verify).mkdir()
    for executable in ([console_path],[str(python),"-m","proteomics_pipeline"]):
        help_result=subprocess.run(executable+["--help"],cwd=empty,text=True,capture_output=True,env=env, encoding="utf-8")
        version_result=subprocess.run(executable+["--version"],cwd=empty,text=True,capture_output=True,env=env, encoding="utf-8")
        assert help_result.returncode==0 and version_result.returncode==0 and version_result.stdout.strip()=="0.1.0.dev0"
        # Amendment A-2026-10-01-03: a run now leaves an honest FAILED run directory, so each executable gets its own output root.
        tag=Path(executable[0]).name
        commands=[(["doctor","--json"],3),(["validate","--config",config,"--schema-only","--json"],0),(["plan","--config",config,"--output",str(empty/f"plan-{tag}"/"plan.json"),"--json"],3),(["run","--config",config,"--output",str(empty/f"run-{tag}"),"--json"],3),(["verify","--run",verify,"--json"],5)]
        for command,expected in commands:
            result=subprocess.run(executable+command,cwd=empty,text=True,capture_output=True,env=env, encoding="utf-8"); assert result.returncode==expected,(executable,command,result.stdout,result.stderr)
