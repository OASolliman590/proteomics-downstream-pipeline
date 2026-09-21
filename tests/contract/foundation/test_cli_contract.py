import json, os, shutil, subprocess, sys
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[3]; EXAMPLES=ROOT/"configs"/"examples"
def run_cli(*args):
    return subprocess.run([sys.executable,"-m","proteomics_pipeline",*args],cwd=ROOT,text=True,capture_output=True,env={**os.environ,"PYTHONPATH":str(ROOT/"src")})
def test_help_and_version_are_available_without_r():
    h=run_cli("--help"); v=run_cli("--version")
    assert h.returncode==0 and v.returncode==0; assert "doctor" in h.stdout and "validate" in h.stdout; assert v.stdout.strip()=="0.1.0.dev0"
@pytest.mark.parametrize("name",["example-independent.json","example-paired.json","example-effect-threshold.json"])
def test_schema_only_examples_resolve(name):
    result=run_cli("validate","--config",str(EXAMPLES/name),"--schema-only","--json")
    assert result.returncode==0,result.stderr; payload=json.loads(result.stdout); assert payload["valid"] is True; assert payload["config"]["primary_engine"]=="limma"; assert payload["config"]["response_mode"]=="descriptive_only"; assert payload["config"]["score_test"]=="off"
def test_unknown_fields_have_json_pointer_and_typed_error(tmp_path):
    source=json.loads((EXAMPLES/"example-independent.json").read_text()); source["study"]["typo"]=True; path=tmp_path/"bad.json"; path.write_text(json.dumps(source))
    result=run_cli("validate","--config",str(path),"--schema-only","--json"); assert result.returncode==2; error=json.loads(result.stdout)["errors"][0]; assert error["code"]=="E_CONFIG_SCHEMA"; assert error["pointer"]=="/study"
def test_full_validation_and_plan_are_explicitly_unavailable():
    config=str(EXAMPLES/"example-independent.json")
    for command in ("validate","plan","run"):
        args=[command,"--config",config,"--json"]
        if command=="validate":args.append("--no-schema-only")
        if command=="plan":args += ["--output","plan.json"]
        if command=="run":args += ["--output","run"]
        result=run_cli(*args); assert result.returncode==3; assert "E_CAPABILITY_NOT_IMPLEMENTED" in result.stdout
def test_unknown_command_is_usage_error(): assert run_cli("not-a-command").returncode==2
def test_installed_console_and_python_module_work_without_source_path(tmp_path):
    env=os.environ.copy(); env.pop("PYTHONPATH",None); python=Path(sys.executable); console_path=shutil.which("proteomics",path=str(python.parent)); assert console_path
    empty=tmp_path/"empty"; empty.mkdir(); env["PATH"]=""; config=str(EXAMPLES/"example-independent.json"); verify=str(empty/"verify-run"); Path(verify).mkdir()
    for executable in ([console_path],[str(python),"-m","proteomics_pipeline"]):
        help_result=subprocess.run(executable+["--help"],cwd=empty,text=True,capture_output=True,env=env)
        version_result=subprocess.run(executable+["--version"],cwd=empty,text=True,capture_output=True,env=env)
        assert help_result.returncode==0 and version_result.returncode==0 and version_result.stdout.strip()=="0.1.0.dev0"
        commands=[(["doctor","--json"],3),(["validate","--config",config,"--schema-only","--json"],0),(["plan","--config",config,"--output",str(empty/"plan.json"),"--json"],3),(["run","--config",config,"--output",str(empty/"run"),"--json"],3),(["verify","--run",verify,"--json"],5)]
        for command,expected in commands:
            result=subprocess.run(executable+command,cwd=empty,text=True,capture_output=True,env=env); assert result.returncode==expected,(executable,command,result.stdout,result.stderr)
