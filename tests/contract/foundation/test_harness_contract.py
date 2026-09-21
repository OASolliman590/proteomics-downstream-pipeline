import json,re,subprocess,sys
from pathlib import Path
import pytest
from scripts.maintained import check_baseline

def test_intentionally_failing_python_selection_is_nonzero(tmp_path):
    test_file=tmp_path/"test_fail.py"; test_file.write_text("def test_failure(): assert False\n")
    result=subprocess.run([sys.executable,"-m","pytest","-q",str(test_file)],capture_output=True,text=True)
    assert result.returncode!=0

@pytest.mark.parametrize("rel",["pipeline/a.txt","legacy/b.txt","docs/audit/c.txt"])
def test_baseline_mutation_detected_in_each_protected_root(tmp_path,rel):
    root=tmp_path/"root"; (root/"pipeline").mkdir(parents=True); (root/"legacy").mkdir(); (root/"docs/audit").mkdir(parents=True)
    registry={"schema_version":"1.0.0","protected_roots":["pipeline","legacy","docs/audit"],"files":[]}
    path=root/rel; path.write_text("original"); registry["files"].append({"path":rel,"sha256":check_baseline.sha256(path)})
    registry_path=tmp_path/"registry.json"; registry_path.write_text(json.dumps(registry))
    assert check_baseline.check(registry_path,root)[0]
    path.write_text("mutated")
    assert not check_baseline.check(registry_path,root)[0]

def test_regenerated_scratch_registry_is_untrusted(tmp_path):
    root=tmp_path/"root"; (root/"pipeline").mkdir(parents=True); path=root/"pipeline"/"a.txt"; path.write_text("original")
    registry_path=tmp_path/"registry.json"; registry={"schema_version":"1.0.0","protected_roots":["pipeline"],"files":[{"path":"pipeline/a.txt","sha256":check_baseline.sha256(path)}]}; registry_path.write_text(json.dumps(registry))
    assert check_baseline.check(registry_path,root,enforce_trust=True)[0] is False

def test_bootstraps_use_local_gates_and_safe_bash_parameter_expansion():
    root=Path(__file__).parents[3]
    bash=(root/"scripts/maintained/bootstrap.sh").read_text(); powershell=(root/"scripts/maintained/bootstrap.ps1").read_text()
    assert "PYTHON_BIN=${PYTHON_BIN:-python3}" in bash and r"PYTHON_BIN=\${PYTHON_BIN" not in bash
    required=("--help","--version","--schema-only","doctor","unittest","pytest","check_baseline.py","write_solved_versions.py","Rscript","test_r.R")
    for gate in required: assert gate in bash and gate in powershell
    assert bash.index("R CMD INSTALL") < bash.index("doctor")
    assert powershell.index("Invoke-NativeChecked 'R'") < powershell.index("doctor")
    assert 'export R_LIBS_USER="$ROOT/.r-lib"' in bash and '$env:R_LIBS_USER=Join-Path $Root' in powershell
    assert 'LASTEXITCODE' in powershell and 'function Invoke-NativeChecked' in powershell
    assert bash.index("R CMD INSTALL") < bash.rindex("write_solved_versions.py")
    assert powershell.index("Invoke-NativeChecked 'R'") < powershell.rindex("write_solved_versions.py")
    r_harness=(root/"scripts/maintained/test_r.R").read_text(); solved=(root/"scripts/maintained/write_solved_versions.py").read_text()
    assert "length(result)==0L" in r_harness and "status=3L" in r_harness
    for package_name in ("jsonlite","openssl","testthat","proteomicsCore"): assert package_name in solved
    dispatch=(root/"r/proteomicsCore/R/dispatch.R").read_text(); assert 'na="null",null="null"' in dispatch and 'required_fields' in dispatch
    wrapper=(root/"scripts/maintained/run_stage.R").read_text(); assert 'quit(status=3L)' in wrapper and 'result$exit_code' in wrapper

def test_solved_versions_persist_only_privacy_safe_executable_descriptors(tmp_path):
    from scripts.maintained import write_solved_versions
    output=tmp_path/"solved-versions.json"; write_solved_versions.main([str(output)])
    raw=output.read_text(encoding="utf-8"); payload=json.loads(raw)
    assert payload["python"]["path"]=="<redacted>" and payload["python"]["executable"]
    assert str(Path.home()).lower() not in raw.lower()
    assert not re.search(r"(?:[A-Za-z]:[\\/]|/home/|/Users/|\\\\Users\\)",raw)
