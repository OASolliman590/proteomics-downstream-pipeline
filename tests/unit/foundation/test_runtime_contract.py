import hashlib,sys,types,shutil,json,subprocess
from pathlib import Path
import pytest
from proteomics_pipeline import provenance,runtime
def test_canonical_hash_is_key_order_independent(): assert provenance.canonical_json_sha256({"b":2,"a":1})==provenance.canonical_json_sha256({"a":1,"b":2})
def test_plan_hash_oracle_excludes_plan_hash_and_observational_timestamps_recursively():
    value={"plan_hash":"forged","created_at":"a","nested":{"started_at":"b","finished_at":"c","value":1},"value":2}
    oracle={"nested":{"value":1},"value":2}
    assert provenance.canonical_plan_bytes(value)==provenance.canonical_json_bytes(oracle)
    assert provenance.canonical_json_sha256(value)==provenance.canonical_json_sha256(oracle)
    assert provenance.canonical_config_sha256(value)!=provenance.canonical_config_sha256(oracle)
def test_publish_artifact_verifies_hash_and_promotes_atomically(tmp_path):
    temp=tmp_path/"stage"; temp.mkdir(); source=temp/"value.txt"; source.write_text("Ω\n",encoding="utf-8"); digest=hashlib.sha256(source.read_bytes()).hexdigest(); destination=tmp_path/"run"/"value.txt"; runtime.publish_artifact(source,destination,digest); assert destination.read_text(encoding="utf-8")=="Ω\n"
    with pytest.raises(runtime.IntegrityError):runtime.publish_artifact(destination,tmp_path/"run"/"other.txt","0"*64)
def test_capabilities_use_fixed_map_and_absent_handlers_are_unavailable():
    values={item["id"]:item for item in runtime.capabilities()}; assert values["foundation.io_roundtrip"]["implemented"] is True; assert values["foundation.io_roundtrip"]["required_packages"]==["Rscript","jsonlite","openssl","proteomicsCore"]; assert values["limma"]["implemented"] is False; assert values["limma"]["module"]=="proteomics_pipeline.inference_service"
def test_capabilities_lazily_activate_a_real_fixed_module_interface(monkeypatch):
    module=types.ModuleType("proteomics_pipeline.inference_service"); module.capabilities=lambda:[{"id":"limma","dependencies":["design"],"required_r_packages":["futureR"]}]; module.execute=lambda request:{}
    monkeypatch.setitem(sys.modules,"proteomics_pipeline.inference_service",module)
    values={item["id"]:item for item in runtime.capabilities()}; assert values["limma"]["implemented"] is True and values["limma"]["dependencies"]==["design"] and values["limma"]["required_r_packages"]==["futureR"]
    invalid=types.ModuleType("proteomics_pipeline.preprocessing_service"); invalid.execute=lambda request:{}
    monkeypatch.setitem(sys.modules,"proteomics_pipeline.preprocessing_service",invalid)
    values={item["id"]:item for item in runtime.capabilities()}; assert values["preprocessing"]["implemented"] is False
    mismatched=types.ModuleType("proteomics_pipeline.design_service"); mismatched.capabilities=lambda:[{"id":"limma"}]; mismatched.execute=lambda request:{}
    monkeypatch.setitem(sys.modules,"proteomics_pipeline.design_service",mismatched)
    values={item["id"]:item for item in runtime.capabilities()}; assert values["design"]["implemented"] is False

def test_execute_capability_requires_stage_protocol_and_identity(monkeypatch,tmp_path):
    module=types.ModuleType("proteomics_pipeline.inference_service"); module.capabilities=lambda:[{"id":"limma"}]
    module.execute=lambda request: runtime.stage_result(request["run_id"],request["stage_id"],request["capability"],"INAPPLICABLE",plan_hash=request["plan_hash"],exit_code=0,reason_code="E_INAPPLICABLE_TEST")
    monkeypatch.setitem(sys.modules,"proteomics_pipeline.inference_service",module)
    request=_request(tmp_path/"out"); request["capability"]="limma"
    result=runtime.execute_capability("limma",request); assert result["state"]=="INAPPLICABLE" and result["exit_code"]==0
    module.execute=lambda request:{"state":"DELEGATED"}
    with pytest.raises((ValueError,runtime.IntegrityError)): runtime.execute_capability("limma",request)

def test_execute_capability_rejects_argument_request_mismatch_before_handler(monkeypatch,tmp_path):
    called=[]; module=types.ModuleType("proteomics_pipeline.inference_service"); module.capabilities=lambda:[{"id":"limma"}]; module.execute=lambda request: called.append(request)
    monkeypatch.setitem(sys.modules,"proteomics_pipeline.inference_service",module)
    request=_request(tmp_path/"out"); request["capability"]="design"
    with pytest.raises(runtime.IntegrityError): runtime.execute_capability("limma",request)
    assert called==[]

def test_runtime_validators_fail_closed_without_jsonschema(monkeypatch,tmp_path):
    monkeypatch.setitem(sys.modules,"jsonschema",None)
    request=_request(tmp_path/"out"); result=runtime.stage_result("r","s","foundation.io_roundtrip","NOT_RUN",plan_hash=None,exit_code=3,reason_code="E_CAPABILITY_NOT_AVAILABLE")
    for validator,value,args in ((runtime.validate_stage_request,request,()),(runtime.validate_stage_result,result,()),(runtime.validate_run_status,{},(tmp_path,))):
        with pytest.raises(runtime.CapabilityError) as exc: validator(value,*args)
        assert exc.value.code=="E_VALIDATOR_UNAVAILABLE" and exc.value.exit_code==3
def test_stage_state_rejects_nonzero_completed_result():
    record=runtime.stage_result("r","s","foundation.io_roundtrip","COMPLETED",exit_code=7)
    with pytest.raises(ValueError):runtime.validate_stage_result(record)
def test_exclusive_run_lock_allows_one_writer(tmp_path):
    lock=tmp_path/"run.lock"; first=runtime.RunLock(lock); second=runtime.RunLock(lock); first.acquire()
    try:
        with pytest.raises(runtime.CollisionError):second.acquire()
    finally:first.release()

def test_subprocess_uses_literal_argv_and_captures_streams(tmp_path):
    literal=tmp_path/"space Ω & sentinel.txt"
    code="import pathlib,sys; pathlib.Path(sys.argv[1]).write_text('stdout'); print('path-ok'); print('stderr',file=sys.stderr)"
    import os
    result=runtime.run_subprocess([sys.executable,"-c",code,str(literal)],cwd=tmp_path)
    assert result.returncode==0 and "path-ok" in result.stdout and "stderr" in result.stderr
    assert literal.read_text()=="stdout"

def test_real_child_failure_preserves_exit_code():
    result=runtime.run_subprocess([sys.executable,"-c","import sys; sys.exit(7)"])
    assert result.returncode==7 and result.stdout=="" and result.stderr==""

def _request(output_dir):
    return {"schema_version":"1.2.0","run_id":"r","stage_id":"s","capability":"foundation.io_roundtrip","plan_hash":None,"inputs":[],"output_temp_dir":str(output_dir),"config_path":None,"parameters":{},"rng":{"seed":1,"kind":"L'Ecuyer-CMRG","threads":1}}

def test_missing_rscript_is_typed_not_run_without_result_or_session(tmp_path):
    result=runtime.execute_stage(_request(tmp_path/"out"),rscript=str(tmp_path/"missing-Rscript"))
    assert result["state"]=="NOT_RUN" and result["reason_code"]=="E_CAPABILITY_NOT_AVAILABLE" and result["exit_code"]==3
    assert not (tmp_path/"out"/"stage-result.json").exists() and result["session_info_path"] is None

def test_timeout_is_failed_and_never_completed(tmp_path,monkeypatch):
    import subprocess
    def timeout(*_args,**_kwargs): raise subprocess.TimeoutExpired(["Rscript"],0.01,output="partial-out",stderr="partial-err")
    monkeypatch.setattr(runtime,"run_subprocess",timeout)
    result=runtime.execute_stage(_request(tmp_path/"out"),rscript="Rscript",timeout=0.01)
    assert result["state"]=="CANCELLED" and result["exit_code"]==6 and result["reason_code"]=="E_CHILD_TIMEOUT"
    assert (tmp_path/"out"/"stdout.log").read_text()=="partial-out" and (tmp_path/"out"/"stderr.log").read_text()=="partial-err"

def test_nonzero_child_failure_is_preserved(tmp_path,monkeypatch):
    import subprocess
    def child(argv,**_kwargs):
        record=runtime.stage_result("r","s","foundation.io_roundtrip","FAILED",plan_hash=None,exit_code=7,reason_code="E_CHILD_EXIT")
        Path(argv[-1]).write_text(__import__("json").dumps(record))
        return subprocess.CompletedProcess(argv,7,"child-out","child-err")
    monkeypatch.setattr(runtime,"run_subprocess",child)
    result=runtime.execute_stage(_request(tmp_path/"out"),rscript="Rscript")
    assert result["state"]=="FAILED" and result["exit_code"]==7 and result["reason_code"]=="E_CHILD_EXIT"

def test_real_execute_stage_child_failure_without_mock(tmp_path):
    result=runtime.execute_stage(_request(tmp_path/"out"),rscript=sys.executable,wrapper=str(tmp_path/"missing wrapper with Ω & ;.R"))
    assert result["state"]=="FAILED" and result["exit_code"]!=0 and result["reason_code"]=="E_CHILD_EXIT"
    assert (tmp_path/"out"/"stderr.log").is_file() and not (tmp_path/"out"/"stage-result.json").exists()

def test_real_r_execute_stage_integration_when_local_package_exists(tmp_path):
    rscript=shutil.which("Rscript")
    if not rscript: pytest.skip("NOT_RUN: Rscript unavailable")
    check=subprocess.run([rscript,"--vanilla","-e","if(requireNamespace('proteomicsCore', quietly=TRUE)) quit(status=0L) else quit(status=1L)"],shell=False,check=False)
    if check.returncode!=0: pytest.skip("NOT_RUN: local proteomicsCore unavailable")
    source=tmp_path/"input Ω space & ;.tsv"; metadata=tmp_path/"metadata.tsv"
    source.write_text("feature_id\tobs-1\tobs-2\n0001\t1.25\tNA\n",encoding="utf-8")
    metadata.write_text("observation_id\tlabel\nobs-1\tα\nobs-2\tβ\n",encoding="utf-8")
    digest=lambda path: __import__("hashlib").sha256(path.read_bytes()).hexdigest()
    payload="run Ω space & copy NUL sentinel.txt" if __import__("os").name=="nt" else "run Ω space ; touch sentinel.txt"
    request=_request(tmp_path/payload/"stage-temp")
    request["inputs"]=[{"artifact_id":"matrix","path":str(source),"sha256":digest(source)},{"artifact_id":"metadata","path":str(metadata),"sha256":digest(metadata)}]
    request["parameters"]={"matrix_input_id":"matrix","metadata_input_id":"metadata"}
    result=runtime.execute_stage(request,rscript=rscript,run_root=tmp_path/payload,cwd=tmp_path)
    assert result["state"]=="COMPLETED" and result["exit_code"]==0
    promoted=tmp_path/payload/"stages"/"s"; assert promoted.joinpath("sessionInfo.txt").is_file()
    assert promoted.joinpath("matrix.tsv").read_text(encoding="utf-8").splitlines()==["feature_id\tobs-1\tobs-2","0001\t1.25\tNA"]
    assert promoted.joinpath("metadata.tsv").read_text(encoding="utf-8").splitlines()==["observation_id\tlabel","obs-1\tα","obs-2\tβ"]
    sentinel=tmp_path/"outside sentinel Ω & ;.txt"; assert not sentinel.exists()
    assert not (tmp_path/"sentinel.txt").exists()
    failed=subprocess.run([rscript,"--vanilla","-e","stop('R01 intentional failure')"],capture_output=True,text=True,shell=False,check=False)
    assert failed.returncode!=0 and "R01 intentional failure" in failed.stderr

def test_real_r_wrapper_failure_is_failed_with_captured_stderr(tmp_path):
    rscript=shutil.which("Rscript")
    if not rscript: pytest.skip("NOT_RUN: Rscript unavailable")
    wrapper=tmp_path/"failing Ω wrapper.R"; wrapper.write_text("stop('R01 wrapper failure')\n",encoding="utf-8")
    result=runtime.execute_stage(_request(tmp_path/"failed output"),rscript=rscript,wrapper=str(wrapper))
    assert result["state"]=="FAILED" and result["exit_code"]!=0 and result["reason_code"]=="E_CHILD_EXIT"
    assert "R01 wrapper failure" in (tmp_path/"failed output"/"stderr.log").read_text(encoding="utf-8")

def test_missing_result_with_capability_exit_three_is_not_run(tmp_path,monkeypatch):
    import subprocess
    monkeypatch.setattr(runtime,"run_subprocess",lambda argv,**kwargs: subprocess.CompletedProcess(argv,3,"","missing package"))
    result=runtime.execute_stage(_request(tmp_path/"out"),rscript="Rscript")
    assert result["state"]=="NOT_RUN" and result["exit_code"]==3 and result["reason_code"]=="E_CAPABILITY_NOT_AVAILABLE"

def test_zero_process_with_nonzero_result_is_integrity_failure(tmp_path,monkeypatch):
    import json,subprocess
    def child(argv,**_kwargs):
        record=runtime.stage_result("r","s","foundation.io_roundtrip","FAILED",plan_hash=None,exit_code=7,reason_code="E_CHILD_EXIT"); Path(argv[-1]).write_text(json.dumps(record)); return subprocess.CompletedProcess(argv,0,"","")
    monkeypatch.setattr(runtime,"run_subprocess",child)
    with pytest.raises(runtime.IntegrityError): runtime.execute_stage(_request(tmp_path/"out"),rscript="Rscript")

def test_execute_stage_requires_fresh_temp_and_atomically_promotes_complete_directory(tmp_path,monkeypatch):
    import json,subprocess
    run_root=tmp_path/"run"; run_root.mkdir(); request=_request(run_root/"stage-temp")
    def child(argv,**_kwargs):
        request_value=json.loads(Path(argv[-3]).read_text()); out=Path(request_value["output_temp_dir"]); (out/"sessionInfo.txt").write_text("session"); (out/"value.txt").write_text("value")
        record=runtime.stage_result("r","s","foundation.io_roundtrip","COMPLETED",plan_hash=None,exit_code=0); record["session_info_path"]="sessionInfo.txt"; record["outputs"]=[{"artifact_id":"session","relative_path":"sessionInfo.txt","sha256":hashlib.sha256((out/"sessionInfo.txt").read_bytes()).hexdigest(),"result_type":"session_info"},{"artifact_id":"value","relative_path":"value.txt","sha256":hashlib.sha256((out/"value.txt").read_bytes()).hexdigest(),"result_type":"text"}]; Path(argv[-1]).write_text(json.dumps(record)); return subprocess.CompletedProcess(argv,0,"","" )
    monkeypatch.setattr(runtime,"run_subprocess",child)
    result=runtime.execute_stage(request,rscript="Rscript",run_root=run_root)
    destination=run_root/"stages"/"s"; assert result["state"]=="COMPLETED" and destination.is_dir() and not Path(request["output_temp_dir"]).exists(); assert (destination/"sessionInfo.txt").exists()
    with pytest.raises(runtime.CollisionError): runtime.execute_stage(request,rscript="Rscript",run_root=run_root)
    stale=tmp_path/"stale"; stale.mkdir()
    with pytest.raises(runtime.CollisionError): runtime.execute_stage(_request(stale),rscript="Rscript")

def test_incomplete_stage_is_not_promoted_to_completed_destination(tmp_path,monkeypatch):
    import json,subprocess
    run_root=tmp_path/"run"; run_root.mkdir(); request=_request(run_root/"stage-temp")
    def child(argv,**_kwargs):
        out=Path(json.loads(Path(argv[-3]).read_text())["output_temp_dir"]); (out/"diagnostic.txt").write_text("partial")
        record=runtime.stage_result("r","s","foundation.io_roundtrip","FAILED",plan_hash=None,exit_code=7,reason_code="E_CHILD_EXIT"); record["outputs"]=[{"artifact_id":"diagnostic","relative_path":"diagnostic.txt","sha256":hashlib.sha256((out/"diagnostic.txt").read_bytes()).hexdigest(),"result_type":"diagnostic"}]; Path(argv[-1]).write_text(json.dumps(record)); return subprocess.CompletedProcess(argv,7,"", "failed")
    monkeypatch.setattr(runtime,"run_subprocess",child)
    result=runtime.execute_stage(request,rscript="Rscript",run_root=run_root)
    assert result["state"]=="FAILED" and not (run_root/"stages"/"s").exists() and Path(request["output_temp_dir"]).is_dir()

def test_stage_states_distinguish_unavailable_and_scientific_inapplicable():
    unavailable=runtime.stage_result("r","s","foundation.io_roundtrip","NOT_RUN",exit_code=3,reason_code="E_CAPABILITY_NOT_AVAILABLE")
    inapplicable=runtime.stage_result("r","s","design","INAPPLICABLE",exit_code=0,reason_code="E_INAPPLICABLE_NO_ELIGIBLE_GROUPS")
    runtime.validate_stage_result(unavailable); runtime.validate_stage_result(inapplicable)
    assert unavailable["state"]=="NOT_RUN" and inapplicable["state"]=="INAPPLICABLE"
