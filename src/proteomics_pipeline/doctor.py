from __future__ import annotations
import shutil, subprocess, sys, hashlib
from importlib import metadata
from pathlib import Path
from .runtime import CAPABILITY_MAP, capabilities
from .config import load_config
def _r_inventory():
    executable=shutil.which("Rscript")
    if executable is None:return {"name":"Rscript","available":False,"status":"NOT_AVAILABLE","path":None,"version":None}
    try: result=subprocess.run([executable,"--version"],capture_output=True,text=True,shell=False,check=False); version=(result.stdout or result.stderr).strip()
    except OSError as exc:return {"name":"Rscript","available":False,"status":"NOT_AVAILABLE","path":executable,"version":None,"error":str(exc)}
    return {"name":"Rscript","available":result.returncode==0,"status":"AVAILABLE" if result.returncode==0 else "NOT_AVAILABLE","path":executable,"version":version}
def _requested(config):
    merged={}
    def add(capability,required,source):
        item=merged.setdefault(capability,{"id":capability,"required":False,"sources":[]})
        item["required"] = bool(item["required"] or required)
        item["sources"].append(source)
    add("foundation.io_roundtrip",True,"foundation")
    add("intake",True,"config.inputs")
    if "preprocessing" in config:
        add("preprocessing",True,"config.preprocessing")
    if "design" in config or config.get("additional_designs"):
        add("design",True,"config.design")
    for model in config.get("models",[]):
        engine=model.get("engine")
        capability={"limma":"limma","deqms":"assay_engines","proda":"assay_engines"}.get(engine,engine)
        if capability in CAPABILITY_MAP:
            add(capability,model.get("execution_requirement")=="required",f"model:{model.get('id')}")
    if config.get("resources"):
        add("resources",True,"config.resources")
    if config.get("pathways",{}).get("enabled"):
        add("pathways",True,"config.pathways")
    if config.get("response",{}).get("enabled"):
        add("response",True,"config.response")
    if config.get("report"):
        add("report_stub",False,"config.report")
    return list(merged.values())
def inspect(*,config_path: str|Path|None=None):
    r=_r_inventory(); python={"name":"python","available":True,"status":"AVAILABLE","path":sys.executable,"version":sys.version.split()[0]}
    discovered=capabilities()
    package_names=[]
    for item in discovered:
        for name in item.get("required_r_packages",[]):
            if name not in package_names: package_names.append(name)
    for name in ("jsonlite","openssl","proteomicsCore"):
        if name not in package_names: package_names.append(name)
    packages=[{"name":name,"available":False,"status":"NOT_AVAILABLE","version":None,"reason":"Rscript is unavailable"} for name in package_names]
    if r["available"]:
        packages=[]
        for package_name in package_names:
            result=subprocess.run([r["path"],"--vanilla","-e",f"if(requireNamespace('{package_name}', quietly=TRUE)) quit(status=0L) else quit(status=1L)"],capture_output=True,shell=False,check=False)
            version_result=subprocess.run([r["path"],"--vanilla","-e",f"if(requireNamespace('{package_name}', quietly=TRUE)) cat(as.character(packageVersion('{package_name}')))"],capture_output=True,text=True,shell=False,check=False)
            packages.append({"name":package_name,"available":result.returncode==0,"status":"AVAILABLE" if result.returncode==0 else "NOT_AVAILABLE","version":version_result.stdout.strip() or None})
    package_by_name={item["name"]:item for item in packages}
    discovered_by_id={item["id"]:item for item in discovered}
    availability={}
    def is_available(capability,stack=()):
        if capability in availability: return availability[capability]
        item=discovered_by_id.get(capability)
        if not item or not item["implemented"] or capability in stack:
            availability[capability]=False; return False
        deps_ok=all(is_available(dep,stack+(capability,)) for dep in item.get("dependencies",[]))
        packages_ok=all(r["available"] and package_by_name.get(pkg,{}).get("available",False) for pkg in item.get("required_r_packages",[]))
        availability[capability]=deps_ok and packages_ok
        return availability[capability]
    caps=[]
    for item in discovered:
        entry=dict(item); entry["available"]=is_available(item["id"])
        entry["status"]="AVAILABLE" if entry["available"] else ("NOT_IMPLEMENTED" if not item["implemented"] else "NOT_AVAILABLE")
        missing_deps=[dep for dep in item.get("dependencies",[]) if not availability.get(dep,False)]
        missing_packages=[pkg for pkg in item.get("required_r_packages",[]) if not r["available"] or not package_by_name.get(pkg,{}).get("available",False)]
        if missing_deps: entry["missing_dependencies"]=missing_deps
        if missing_packages: entry["missing_r_packages"]=missing_packages
        caps.append(entry)
    resources=[]; requested=[]
    if config_path:
        try:
            config=load_config(config_path)
            requested=_requested(config)
            config_root=Path(config_path).resolve().parent
            for item in config.get("resources",[]):
                path=(config_root/item["path"]).resolve()
                observed=None
                if path.is_file():
                    digest=hashlib.sha256(path.read_bytes()).hexdigest(); observed=digest
                expected=item.get("sha256")
                resources.append({"id":item["id"],"path":str(path),"available":path.is_file() and (expected is None or observed==expected),"status":"AVAILABLE" if path.is_file() and (expected is None or observed==expected) else "NOT_AVAILABLE","sha256":observed,"expected_sha256":expected})
        except Exception as exc:
            resources.append({"available":False,"status":"NOT_AVAILABLE","reason":str(exc)})
    cap_by_id={item["id"]:item for item in caps}
    requested_report=[dict(item,available=bool(cap_by_id.get(item["id"],{}).get("available",False)),status=cap_by_id.get(item["id"],{}).get("status","NOT_IMPLEMENTED"),required_packages=cap_by_id.get(item["id"],{}).get("required_packages",[])) for item in requested]
    return {"python":python,"r":r,"packages":packages,"resources":resources,"capabilities":caps,"requested_capabilities":requested_report,"config":str(config_path) if config_path else None}
def exit_code(report):
    foundation=next(item for item in report["capabilities"] if item["id"]=="foundation.io_roundtrip")
    if any(item.get("required") and not item.get("available") for item in report.get("requested_capabilities",[])): return 3
    return 0 if foundation["available"] and all(item.get("available",False) for item in report.get("resources",[])) else 3
