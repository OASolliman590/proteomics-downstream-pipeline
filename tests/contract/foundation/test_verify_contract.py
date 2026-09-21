import hashlib,json,sys
from pathlib import Path
from proteomics_pipeline import runtime
from proteomics_pipeline.cli import main
from scripts.maintained import check_baseline

def test_default_baseline_registry_uses_committed_lf_bytes():
    registry = check_baseline.REGISTRY
    canonical = registry.read_bytes().replace(b"\r\n", b"\n")
    expected = "289afd374e6d0159eb181526a9724c80833facb195f42d6a427fc4b15a90d9bd"
    assert hashlib.sha256(canonical).hexdigest() == expected
    assert check_baseline.TRUSTED_REGISTRY_SHA256 == expected
    ok, mismatches = check_baseline.check(registry, check_baseline.ROOT, enforce_trust=True)
    assert ok, mismatches

def _status(root, state="FAILED", exit_code=7):
    artifact=root/"value.txt"; artifact.write_text("real")
    digest=hashlib.sha256(artifact.read_bytes()).hexdigest()
    stage=runtime.stage_result("run-1","stage-1","foundation.io_roundtrip",state,exit_code=exit_code,reason_code="E_CHILD_EXIT" if state!="COMPLETED" else None,message="child")
    (root/"stages").mkdir(exist_ok=True); (root/"stages"/"stage-result.json").write_text(json.dumps(stage))
    value={"schema_version":"1.2.0","run_id":"run-1","plan_hash":None,"requested_phase":1,"implemented_capabilities":[],"state":"FAILED","reason_code":"E_CHILD_EXIT","started_at":None,"finished_at":None,"exit_code":7,"stages":[{"stage_id":"stage-1","capability":"foundation.io_roundtrip","required":True,"state":state,"reason_code":stage["reason_code"],"result_path":"stages/stage-result.json"}],"artifacts":[{"artifact_id":"value","relative_path":"value.txt","sha256":digest,"result_type":"text"}],"warnings_file":"warnings.json"}
    (root/"run_status.json").write_text(json.dumps(value))

def test_verify_accepts_failed_child_and_rejects_mutation(tmp_path):
    _status(tmp_path)
    assert main(["verify","--run",str(tmp_path),"--json"])==0
    (tmp_path/"value.txt").write_text("mutated")
    assert main(["verify","--run",str(tmp_path),"--json"])==5

def test_verify_rejects_forged_completed_nonzero_and_empty_completion(tmp_path):
    _status(tmp_path,state="COMPLETED",exit_code=7)
    assert main(["verify","--run",str(tmp_path),"--json"])==5
    root=tmp_path/"empty"; root.mkdir(); value={"schema_version":"1.2.0","run_id":"run-1","plan_hash":None,"requested_phase":1,"implemented_capabilities":[],"state":"COMPLETED","reason_code":None,"started_at":None,"finished_at":"2026-01-01T00:00:00Z","exit_code":0,"stages":[],"artifacts":[],"warnings_file":"warnings.json"}; (root/"run_status.json").write_text(json.dumps(value))
    assert main(["verify","--run",str(root),"--json"])==5

def test_verify_maps_missing_jsonschema_to_typed_exit(monkeypatch,tmp_path,capsys):
    (tmp_path/"run_status.json").write_text("{}")
    monkeypatch.setitem(sys.modules,"jsonschema",None)
    assert main(["verify","--run",str(tmp_path),"--json"])==3
    captured=capsys.readouterr(); payload=json.loads(captured.out); assert payload["errors"][0]["code"]=="E_VALIDATOR_UNAVAILABLE" and "Traceback" not in captured.err

def test_verify_rejects_result_run_id_and_plan_hash_mismatch(tmp_path):
    _status(tmp_path)
    result_path=tmp_path/"stages"/"stage-result.json"
    result=json.loads(result_path.read_text()); result["run_id"]="other-run"; result_path.write_text(json.dumps(result))
    assert main(["verify","--run",str(tmp_path),"--json"])==5
    _status(tmp_path)
    result=json.loads(result_path.read_text()); result["plan_hash"]="0"*64; result_path.write_text(json.dumps(result))
    assert main(["verify","--run",str(tmp_path),"--json"])==5

def test_run_lock_is_exclusive_across_processes(tmp_path):
    lock=runtime.RunLock(tmp_path/"run.lock"); lock.acquire()
    try:
        import subprocess,sys
        code="from proteomics_pipeline.runtime import RunLock; RunLock(r'"+str(tmp_path/"run.lock")+"').acquire()"
        result=subprocess.run([sys.executable,"-c",code],capture_output=True)
        assert result.returncode!=0
    finally: lock.release()

def test_two_subprocess_writers_and_interrupted_owner_leave_no_promotion(tmp_path):
    import subprocess, time
    lock_path=tmp_path/"run.lock"; destination=tmp_path/"run"/"value.txt"; ready=tmp_path/"owner.ready"
    code="""import pathlib,sys,time,hashlib
from proteomics_pipeline.runtime import RunLock,publish_artifact
lock=pathlib.Path(sys.argv[1]); destination=pathlib.Path(sys.argv[2]); ready=pathlib.Path(sys.argv[3]); role=sys.argv[4]
try:
    with RunLock(lock):
        if role=='owner':
            ready.write_text('ready'); time.sleep(30); source=destination.with_suffix('.source'); source.parent.mkdir(parents=True,exist_ok=True); source.write_text('payload'); publish_artifact(source,destination,hashlib.sha256(source.read_bytes()).hexdigest())
except Exception as exc:
    print(type(exc).__name__, flush=True); raise
"""
    owner=subprocess.Popen([sys.executable,"-c",code,str(lock_path),str(destination),str(ready),"owner"],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    try:
        for _ in range(50):
            if ready.exists(): break
            time.sleep(0.02)
        contender=subprocess.run([sys.executable,"-c",code,str(lock_path),str(destination),str(ready),"contender"],capture_output=True,text=True,timeout=5)
        assert contender.returncode!=0 and "CollisionError" in (contender.stdout+contender.stderr)
        owner.terminate(); owner.wait(timeout=5)
        assert not destination.exists() and not destination.with_suffix(".cache").exists()
    finally:
        if owner.poll() is None: owner.kill(); owner.wait(timeout=5)

def _completed(root,optional_state=None):
    stage_dir=root/"stages"; stage_dir.mkdir(parents=True,exist_ok=True); session=stage_dir/"sessionInfo.txt"; session.write_text("real session")
    value=stage_dir/"value.txt"; value.write_text("value")
    stage=runtime.stage_result("run-1","required-stage","foundation.io_roundtrip","COMPLETED",plan_hash=None,exit_code=0); stage["session_info_path"]="sessionInfo.txt"; stage["outputs"]=[]
    for artifact,result_type in ((session,"session_info"),(value,"text")):
        stage["outputs"].append({"artifact_id":artifact.name,"relative_path":artifact.name,"sha256":hashlib.sha256(artifact.read_bytes()).hexdigest(),"result_type":result_type})
    (stage_dir/"stage-result.json").write_text(json.dumps(stage))
    stages=[{"stage_id":"required-stage","capability":"foundation.io_roundtrip","required":True,"state":"COMPLETED","reason_code":None,"result_path":"stages/stage-result.json"}]
    state="COMPLETED"; exit_code=0; reason=None
    if optional_state:
        optional=runtime.stage_result("run-1","optional-stage","design",optional_state,plan_hash=None,exit_code=7,reason_code="E_OPTIONAL_FAILURE"); (stage_dir/"optional-result.json").write_text(json.dumps(optional)); stages.append({"stage_id":"optional-stage","capability":"design","required":False,"state":optional_state,"reason_code":optional["reason_code"],"result_path":"stages/optional-result.json"}); state="PARTIAL"; exit_code=3; reason="E_OPTIONAL_FAILURE"
    status={"schema_version":"1.2.0","run_id":"run-1","plan_hash":None,"requested_phase":1,"implemented_capabilities":[],"state":state,"reason_code":reason,"started_at":"2026-01-01T00:00:00Z","finished_at":"2026-01-01T00:00:01Z","exit_code":exit_code,"stages":stages,"artifacts":[],"warnings_file":"warnings.json"}; (root/"run_status.json").write_text(json.dumps(status))

def test_verify_enforces_aggregate_completed_and_partial_state_table(tmp_path):
    _completed(tmp_path,optional_state="FAILED")
    value=json.loads((tmp_path/"run_status.json").read_text()); value["state"]="COMPLETED"; value["reason_code"]=None; value["exit_code"]=0; (tmp_path/"run_status.json").write_text(json.dumps(value)); assert main(["verify","--run",str(tmp_path),"--json"])==5
    _completed(tmp_path,optional_state="FAILED"); assert main(["verify","--run",str(tmp_path),"--json"])==0
    value=json.loads((tmp_path/"run_status.json").read_text()); value["state"]="FAILED"; value["reason_code"]="E_OPTIONAL_FAILURE"; (tmp_path/"run_status.json").write_text(json.dumps(value)); assert main(["verify","--run",str(tmp_path),"--json"])==5

def test_python_only_completed_stage_does_not_require_r_session_manifest(tmp_path):
    stages=tmp_path/"stages"; stages.mkdir(); value=stages/"value.txt"; value.write_text("python output")
    stage=runtime.stage_result("run-1","python-stage","report_stub","COMPLETED",plan_hash=None,exit_code=0); stage["outputs"]=[{"artifact_id":"value","relative_path":"value.txt","sha256":hashlib.sha256(value.read_bytes()).hexdigest(),"result_type":"text"}]
    (stages/"stage-result.json").write_text(json.dumps(stage))
    status={"schema_version":"1.2.0","run_id":"run-1","plan_hash":None,"requested_phase":1,"implemented_capabilities":["report_stub"],"state":"COMPLETED","reason_code":None,"started_at":"2026-01-01T00:00:00Z","finished_at":"2026-01-01T00:00:01Z","exit_code":0,"stages":[{"stage_id":"python-stage","capability":"report_stub","required":True,"state":"COMPLETED","reason_code":None,"result_path":"stages/stage-result.json"}],"artifacts":[],"warnings_file":"warnings.json"}
    (tmp_path/"run_status.json").write_text(json.dumps(status)); assert main(["verify","--run",str(tmp_path),"--json"])==0
