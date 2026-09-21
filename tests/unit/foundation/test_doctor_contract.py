import hashlib,json,shutil,subprocess,sys,types
import pytest
from proteomics_pipeline.doctor import inspect,exit_code

def test_doctor_reports_actual_python_identity():
    report=inspect()
    assert report["python"]["path"]==sys.executable
    assert report["python"]["version"]==sys.version.split()[0]
    assert {item["name"] for item in report["packages"]}=={"jsonlite","openssl","proteomicsCore"}

def test_doctor_resource_hash_and_absolute_path_are_truthful(tmp_path):
    resource=tmp_path/"resource.dat"; resource.write_bytes(b"stable-resource")
    digest=hashlib.sha256(resource.read_bytes()).hexdigest()
    source=json.loads((__import__("pathlib").Path(__file__).parents[3]/"configs/examples/example-independent.json").read_text())
    source["runtime"]["phase"]=2
    source["resources"]=[{"id":"r","kind":"mapping","path":str(resource),"sha256":digest,"version":"1","source":"synthetic","source_taxonomy_id":1,"target_taxonomy_id":1,"terms":"term"}]
    config_path=tmp_path/"analysis.json"; config_path.write_text(json.dumps(source))
    report=inspect(config_path=config_path); assert report["resources"][0]["available"] is True; assert exit_code(report)==3
    resource.write_bytes(b"mutated"); report=inspect(config_path=config_path); assert report["resources"][0]["status"]=="NOT_AVAILABLE" and exit_code(report)==3
    resource.unlink(); report=inspect(config_path=config_path); assert report["resources"][0]["status"]=="NOT_AVAILABLE" and exit_code(report)==3

def test_doctor_config_reports_required_limma_unavailable():
    root=__import__("pathlib").Path(__file__).parents[3]; config_path=root/"configs/examples/example-independent.json"
    report=inspect(config_path=config_path)
    required={item["id"] for item in report["requested_capabilities"] if item["required"]}
    assert "limma" in required and any(item["id"]=="limma" and item["status"]=="NOT_IMPLEMENTED" for item in report["requested_capabilities"])
    assert exit_code(report)==3

def test_doctor_activates_only_a_valid_discovered_future_module(monkeypatch):
    root=__import__("pathlib").Path(__file__).parents[3]; config_path=root/"configs/examples/example-independent.json"
    design=types.ModuleType("proteomics_pipeline.design_service"); design.capabilities=lambda:[{"id":"design"}]; design.execute=lambda request:{}
    limma=types.ModuleType("proteomics_pipeline.inference_service"); limma.capabilities=lambda:[{"id":"limma","dependencies":["design"]}]; limma.execute=lambda request:{}
    monkeypatch.setitem(sys.modules,"proteomics_pipeline.design_service",design); monkeypatch.setitem(sys.modules,"proteomics_pipeline.inference_service",limma)
    report=inspect(config_path=config_path); statuses={item["id"]:item for item in report["capabilities"]}
    assert statuses["design"]["status"]=="AVAILABLE" and statuses["limma"]["status"]=="AVAILABLE"

    limma.capabilities=lambda:[{"id":"limma","dependencies":["design"],"required_r_packages":["futureR"]}]
    report=inspect(config_path=config_path); statuses={item["id"]:item for item in report["capabilities"]}
    assert statuses["limma"]["status"]=="NOT_AVAILABLE" and statuses["limma"]["missing_r_packages"]==["futureR"]

def test_doctor_matches_direct_real_runtime_inventory_when_available():
    rscript=shutil.which("Rscript")
    if not rscript: pytest.skip("NOT_RUN: Rscript unavailable")
    for package in ("jsonlite","openssl","proteomicsCore"):
        probe=subprocess.run([rscript,"--vanilla","-e",f"if(requireNamespace('{package}', quietly=TRUE)) quit(status=0L) else quit(status=1L)"],capture_output=True,shell=False,check=False)
        if probe.returncode!=0: pytest.skip("NOT_RUN: required R package unavailable")
    report=inspect(); packages={item["name"]:item for item in report["packages"]}; caps={item["id"]:item for item in report["capabilities"]}
    assert all(packages[name]["available"] for name in ("jsonlite","openssl","proteomicsCore")); assert caps["foundation.io_roundtrip"]["status"]=="AVAILABLE"
