from __future__ import annotations
import importlib, json, os, subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable
from .errors import CapabilityError, CollisionError, IntegrityError, ProteomicsError
from .paths import atomic_write_bytes
from .provenance import sha256_file
STAGE_STATES={"NOT_RUN","RUNNING","COMPLETED","INAPPLICABLE","FAILED","CANCELLED","NOT_REQUESTED"}
R_BACKED_CAPABILITIES={"foundation.io_roundtrip","preprocessing","design","limma","assay_engines","resources","pathways","response"}
CAPABILITY_MAP={
 "foundation.io_roundtrip":("proteomics_pipeline","io_roundtrip"),
 "intake":("proteomics_pipeline.intake.service","execute"),
 "preprocessing":("proteomics_pipeline.preprocessing_service","execute"),
 "design":("proteomics_pipeline.design_service","execute"),
 "limma":("proteomics_pipeline.inference_service","execute"),
 "assay_engines":("proteomics_pipeline.assay_service","execute"),
 "resources":("proteomics_pipeline.resources","execute"),
 "pathways":("proteomics_pipeline.pathway_service","execute"),
 "response":("proteomics_pipeline.response_service","execute"),
 "report_stub":("proteomics_pipeline.reporting.stub","execute"),
 "report_full":("proteomics_pipeline.reporting.full","execute"),
 "reproduction":("proteomics_pipeline.reproduction","execute"),
 "legacy":("proteomics_pipeline.legacy_service","execute")}
REQUIRED_PACKAGES={"foundation.io_roundtrip":["Rscript","jsonlite","openssl","proteomicsCore"]}
def _discovery(capability,module_name,function_name):
    record={"id":capability,"module":module_name,"function":function_name,
            "implemented":False,"dependencies":[],"required_r_packages":[],
            "required_packages":list(REQUIRED_PACKAGES.get(capability,()))}
    try: module=importlib.import_module(module_name)
    except (ImportError,ModuleNotFoundError): return record
    execute=getattr(module,function_name,None); advertise=getattr(module,"capabilities",None)
    if not callable(execute) or not callable(advertise): return record
    try: advertised=advertise()
    except Exception: return record
    if isinstance(advertised,dict): advertised=[advertised]
    if not isinstance(advertised,(list,tuple)) or not advertised or not all(isinstance(item,dict) for item in advertised): return record
    ids=[item.get("id") for item in advertised]
    if any(item_id not in CAPABILITY_MAP for item_id in ids) or len(set(ids))!=len(ids) or capability not in ids: return record
    item=next(item for item in advertised if item.get("id")==capability)
    dependencies=item.get("dependencies",[])
    r_packages=item.get("required_r_packages",item.get("required_packages",[]))
    if not isinstance(dependencies,list) or not all(isinstance(x,str) and x in CAPABILITY_MAP for x in dependencies): return record
    if not isinstance(r_packages,list) or not all(isinstance(x,str) and x for x in r_packages): return record
    record["implemented"]=True
    record["dependencies"]=list(dict.fromkeys(dependencies))
    record["required_r_packages"]=list(dict.fromkeys(r_packages))
    for package in record["required_r_packages"]:
        if package not in record["required_packages"]: record["required_packages"].append(package)
    return record
def _implemented(capability,module_name,function_name):
    """Compatibility predicate backed by the full validated discovery record."""
    return _discovery(capability,module_name,function_name)["implemented"]
def capabilities(): return [_discovery(i,m,f) for i,(m,f) in CAPABILITY_MAP.items()]
def execute_capability(capability, request):
    """Invoke only a discovered handler from the frozen internal map."""
    if not isinstance(request,dict) or request.get("capability")!=capability:
        raise IntegrityError("stage handler capability argument does not match request")
    validate_stage_request(request)
    module_name,function_name=CAPABILITY_MAP.get(capability,(None,None))
    if module_name is None: raise CapabilityError("E_CAPABILITY_NOT_IMPLEMENTED",f"unknown fixed capability: {capability}")
    record=_discovery(capability,module_name,function_name)
    if not record["implemented"]: raise CapabilityError("E_CAPABILITY_NOT_IMPLEMENTED",f"no validated handler for {capability}")
    module=importlib.import_module(module_name)
    result=getattr(module,function_name)(request)
    if not isinstance(result,dict): raise ProteomicsError("E_CAPABILITY_PROTOCOL",f"{capability} handler must return an object",exit_code=4)
    validate_stage_result(result)
    for key in ("run_id","stage_id","capability","plan_hash"):
        if result.get(key)!=request.get(key): raise IntegrityError(f"stage handler changed {key}")
    if result.get("state") not in STAGE_STATES: raise IntegrityError("stage handler returned an unknown state")
    if result.get("state")=="COMPLETED" and result.get("exit_code")!=0: raise IntegrityError("completed stage handler returned nonzero exit")
    return result
def utc_now(): return datetime.now(timezone.utc).isoformat().replace("+00:00","Z")
def stage_result(run_id,stage_id,capability,state,*,plan_hash=None,exit_code=None,reason_code=None,message=""):
    if state not in STAGE_STATES: raise ValueError(f"unknown stage state: {state}")
    return {"schema_version":"1.2.0","run_id":run_id,"stage_id":stage_id,"capability":capability,"plan_hash":plan_hash,"state":state,"reason_code":reason_code,"message":message,"started_at":None if state=="NOT_REQUESTED" else utc_now(),"finished_at":None if state in {"RUNNING","NOT_REQUESTED"} else utc_now(),"exit_code":exit_code,"outputs":[],"warnings":[],"session_info_path":None}
def validate_stage_result(value):
    state=value.get("state")
    if state not in STAGE_STATES: raise ValueError("invalid stage state")
    if state=="COMPLETED" and value.get("exit_code")!=0: raise ValueError("COMPLETED stage must have exit_code=0")
    if state in {"NOT_RUN","FAILED","INAPPLICABLE","CANCELLED"} and not value.get("reason_code"): raise ValueError("non-completed stage requires reason_code")
    if state=="COMPLETED" and value.get("reason_code") is not None: raise ValueError("COMPLETED stage cannot have reason_code")
    errors=list(_draft_validator("stage-result.schema.json").iter_errors(value))
    if errors: raise ValueError(errors[0].message)
def _schema(name):
    return json.loads((Path(__file__).parent/"schemas"/name).read_text(encoding="utf-8"))
def _draft_validator(schema_name):
    try:
        from jsonschema import Draft202012Validator
    except ImportError as exc:
        raise CapabilityError("E_VALIDATOR_UNAVAILABLE","jsonschema is required for strict protocol validation") from exc
    return Draft202012Validator(_schema(schema_name))
def validate_stage_request(value):
    errors=list(_draft_validator("stage-request.schema.json").iter_errors(value))
    if errors: raise ValueError(errors[0].message)
    return True
def validate_run_status(value, root):
    root=Path(root).resolve()
    errors=list(_draft_validator("run-status.schema.json").iter_errors(value))
    if errors: raise ValueError(errors[0].message)
    stages=value.get("stages",[])
    bad_states={"NOT_RUN","FAILED","CANCELLED"}
    required_bad=[stage for stage in stages if stage.get("required") and stage.get("state") in bad_states]
    optional_bad=[stage for stage in stages if not stage.get("required") and stage.get("state") in bad_states]
    if value.get("state")=="COMPLETED":
        if value.get("exit_code")!=0 or not value.get("finished_at") or not value.get("stages"): raise ValueError("completed run must have successful executed stages")
        if required_bad or optional_bad or any(stage.get("state")=="RUNNING" for stage in stages): raise ValueError("completed run contains a nonexecuted or failed requested stage")
        if any(stage.get("required") and stage.get("state")!="COMPLETED" for stage in stages): raise ValueError("required stage is not completed")
        if not any(stage.get("state")=="COMPLETED" for stage in stages): raise ValueError("completed run requires at least one executed stage")
    elif value.get("state")=="PARTIAL":
        if value.get("exit_code")==0 or required_bad or not optional_bad: raise ValueError("partial run requires optional operational failure and nonzero exit")
    elif value.get("state")=="FAILED":
        if value.get("exit_code")==0 or not required_bad: raise ValueError("failed run requires required stage failure and nonzero exit")
    elif value.get("state")=="NOT_RUN":
        if any(stage.get("state")=="COMPLETED" for stage in stages): raise ValueError("NOT_RUN run cannot contain completed stages")
    if value.get("state") in {"FAILED","PARTIAL","CANCELLED","NOT_RUN"} and not value.get("reason_code"): raise ValueError("non-completed run requires reason_code")
    for stage in value.get("stages",[]):
        result_path=stage.get("result_path")
        if result_path is None:
            if stage.get("state")=="COMPLETED": raise ValueError("completed stage is missing result_path")
            continue
        path=(root/result_path).resolve()
        if root not in path.parents or not path.is_file(): raise ValueError("stage result path escapes or is missing")
        result=json.loads(path.read_text(encoding="utf-8")); validate_stage_result(result)
        if result.get("state")=="COMPLETED" and result.get("capability") in R_BACKED_CAPABILITIES and not result.get("session_info_path"): raise ValueError("completed R-backed stage requires session_info_path")
        if result.get("run_id")!=value.get("run_id") or result.get("plan_hash")!=value.get("plan_hash"): raise ValueError("stage result run_id/plan_hash mismatch")
        if result.get("stage_id")!=stage.get("stage_id") or result.get("capability")!=stage.get("capability") or result.get("state")!=stage.get("state"): raise ValueError("stage result identity/state mismatch")
        if result.get("session_info_path"):
            session_path=(path.parent/result["session_info_path"]).resolve()
            if path.parent not in session_path.parents or not session_path.is_file(): raise ValueError("session_info_path escapes stage output or is missing")
        outputs=result.get("outputs",[])
        if result.get("state")=="COMPLETED" and result.get("capability") in R_BACKED_CAPABILITIES and not any(output.get("relative_path")==result.get("session_info_path") for output in outputs): raise ValueError("session_info_path is not in the completed output manifest")
        for output in outputs:
            output_path=(path.parent/output["relative_path"]).resolve()
            if path.parent not in output_path.parents or not output_path.is_file() or sha256_file(output_path)!=output["sha256"]: raise ValueError("stage output is outside its stage or hash-mismatched")
    for artifact in value.get("artifacts",[]):
        path=(root/artifact["relative_path"]).resolve()
        if root not in path.parents or not path.is_file() or sha256_file(path)!=artifact["sha256"]: raise ValueError("run artifact is outside root, missing, or hash-mismatched")
    return True
class RunLock:
    def __init__(self,path): self.path=Path(path); self._fd=None
    def acquire(self):
        self.path.parent.mkdir(parents=True,exist_ok=True)
        try:
            self._fd=os.open(self.path,os.O_CREAT|os.O_EXCL|os.O_WRONLY); os.write(self._fd,f"pid={os.getpid()}\n".encode())
        except FileExistsError as exc: raise CollisionError() from exc
    def release(self):
        if self._fd is not None:
            os.close(self._fd); self._fd=None
            try:self.path.unlink()
            except FileNotFoundError:pass
    def __enter__(self): self.acquire(); return self
    def __exit__(self,*_): self.release()
def publish_artifact(source,destination,expected_sha256):
    source=Path(source); destination=Path(destination)
    if sha256_file(source)!=expected_sha256: raise IntegrityError(f"source hash does not match expected SHA-256: {source}")
    if destination.exists(): raise CollisionError(f"destination already exists: {destination}")
    destination.parent.mkdir(parents=True,exist_ok=True); atomic_write_bytes(destination,source.read_bytes())
    if sha256_file(destination)!=expected_sha256: raise IntegrityError(f"promoted artifact hash mismatch: {destination}")
def run_subprocess(argv:Iterable[str],*,cwd=None,timeout=None):
    return subprocess.run(list(argv),cwd=cwd,shell=False,check=False,capture_output=True,text=True,encoding="utf-8",errors="replace",timeout=timeout)
def execute_stage(request, *, rscript="Rscript", wrapper=None, cwd=None, timeout=None, run_root=None, promoted_stage_dir=None):
    validate_stage_request(request)
    output_dir=Path(request["output_temp_dir"]).resolve()
    if output_dir.exists(): raise CollisionError(f"stage temp already exists: {output_dir}")
    if run_root is not None:
        root=Path(run_root).resolve(); destination=(Path(promoted_stage_dir).resolve() if promoted_stage_dir else root/"stages"/request["stage_id"])
        lock_path=root/".run.lock"
    else:
        destination=Path(promoted_stage_dir).resolve() if promoted_stage_dir else None; lock_path=output_dir.parent/".stage.lock"
    lock_path.parent.mkdir(parents=True,exist_ok=True)
    with RunLock(lock_path):
        if output_dir.exists(): raise CollisionError(f"stage temp already exists: {output_dir}")
        if destination is not None and destination.exists(): raise CollisionError(f"promoted stage already exists: {destination}")
        output_dir.mkdir(parents=True,exist_ok=False)
        request_path=output_dir/"stage-request.json"; result_path=output_dir/"stage-result.json"
        atomic_write_bytes(request_path, json.dumps(request, ensure_ascii=False, sort_keys=True, separators=(",",":")).encode("utf-8"))
        wrapper_path=Path(wrapper) if wrapper else Path(__file__).resolve().parents[2]/"scripts"/"maintained"/"run_stage.R"
        argv=[rscript,"--vanilla",str(wrapper_path),"--request",str(request_path),"--result",str(result_path)]
        try:
            process=run_subprocess(argv,cwd=cwd,timeout=timeout)
        except FileNotFoundError as exc:
            (output_dir/"stdout.log").write_text("",encoding="utf-8"); (output_dir/"stderr.log").write_text(str(exc),encoding="utf-8")
            return stage_result(request["run_id"],request["stage_id"],request["capability"],"NOT_RUN",plan_hash=request["plan_hash"],exit_code=3,reason_code="E_CAPABILITY_NOT_AVAILABLE",message=f"Rscript executable is unavailable: {rscript}")
        except subprocess.TimeoutExpired as exc:
            stdout=exc.stdout.decode("utf-8","replace") if isinstance(exc.stdout,bytes) else (exc.stdout or "")
            stderr=exc.stderr.decode("utf-8","replace") if isinstance(exc.stderr,bytes) else (exc.stderr or "")
            (output_dir/"stdout.log").write_text(stdout,encoding="utf-8"); (output_dir/"stderr.log").write_text(stderr,encoding="utf-8")
            return stage_result(request["run_id"],request["stage_id"],request["capability"],"CANCELLED",plan_hash=request["plan_hash"],exit_code=6,reason_code="E_CHILD_TIMEOUT",message="stage subprocess exceeded its timeout")
        (output_dir/"stdout.log").write_text(process.stdout or "",encoding="utf-8")
        (output_dir/"stderr.log").write_text(process.stderr or "",encoding="utf-8")
        if not result_path.is_file():
            state="NOT_RUN" if process.returncode==3 else "FAILED"; reason="E_CAPABILITY_NOT_AVAILABLE" if state=="NOT_RUN" else "E_CHILD_EXIT"
            return stage_result(request["run_id"],request["stage_id"],request["capability"],state,plan_hash=request["plan_hash"],exit_code=process.returncode,reason_code=reason,message=(process.stderr or "stage did not produce a result"))
        result=json.loads(result_path.read_text(encoding="utf-8")); validate_stage_result(result)
        if process.returncode != result.get("exit_code"): raise IntegrityError("child exit code does not match stage result")
        if result.get("run_id") != request["run_id"] or result.get("stage_id") != request["stage_id"] or result.get("capability") != request["capability"] or result.get("plan_hash") != request["plan_hash"]:
            raise IntegrityError("stage result identity does not match request")
        if result.get("session_info_path"):
            session_path=(output_dir/result["session_info_path"]).resolve()
            if output_dir not in session_path.parents or not session_path.is_file(): raise IntegrityError("session_info_path is outside output_temp_dir or missing")
        for output in result.get("outputs",[]):
            path=(output_dir/output["relative_path"]).resolve()
            if output_dir not in path.parents or not path.is_file() or sha256_file(path)!=output["sha256"]: raise IntegrityError(f"stage output hash/path verification failed: {output.get('relative_path')}")
        if result.get("state")=="COMPLETED" and result.get("capability") in R_BACKED_CAPABILITIES and not any(output.get("relative_path")==result.get("session_info_path") for output in result.get("outputs",[])): raise IntegrityError("session_info_path is not in the completed output manifest")
        if destination is not None and result.get("state")=="COMPLETED":
            destination.parent.mkdir(parents=True,exist_ok=True)
            if destination.exists(): raise CollisionError(f"promoted stage already exists: {destination}")
            os.replace(output_dir,destination)
        return result
