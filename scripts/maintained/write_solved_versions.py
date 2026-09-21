#!/usr/bin/env python3
from __future__ import annotations
import json,platform,shutil,subprocess,sys
from importlib import metadata
from datetime import datetime,timezone
from pathlib import Path
def _safe_executable(path):
 if not path:return None
 name=str(path).replace("\\","/").rsplit("/",1)[-1]
 return {"executable":name,"path":"<redacted>"}
def main(argv=None):
 output=Path(argv[0]) if argv else Path("docs/validation/002-foundation/solved-versions.json"); rscript=shutil.which("Rscript"); r=None; r_packages={name:None for name in ("jsonlite","openssl","testthat","proteomicsCore")}
 if rscript:
  process=subprocess.run([rscript,"--version"],capture_output=True,text=True,check=False,shell=False); packages={}
  for package_name in ("jsonlite","openssl","testthat","proteomicsCore"):
   version_process=subprocess.run([rscript,"--vanilla","-e",f"if(requireNamespace('{package_name}', quietly=TRUE)) cat(as.character(packageVersion('{package_name}')))"],capture_output=True,text=True,check=False,shell=False)
   packages[package_name]=version_process.stdout.strip() or None
  r_packages=packages
  r={**_safe_executable(rscript),"version":(process.stdout or process.stderr).strip(),"exit_code":process.returncode,"packages":packages}
 dependencies={}
 for name in ("setuptools","jsonschema","PyYAML","pytest"):
  try: dependencies[name]=metadata.version(name)
  except metadata.PackageNotFoundError: dependencies[name]=None
 value={"generated_at":datetime.now(timezone.utc).isoformat(),"python":{**_safe_executable(sys.executable),"version":platform.python_version()},"dependencies":dependencies,"r_packages":r_packages,"rscript":r,"status":"AVAILABLE" if r else "NOT_RUN"}; output.parent.mkdir(parents=True,exist_ok=True); output.write_text(json.dumps(value,indent=2,sort_keys=True)+"\n",encoding="utf-8"); print(json.dumps(value,sort_keys=True)); return 0 if r else 3
if __name__=="__main__":raise SystemExit(main())
