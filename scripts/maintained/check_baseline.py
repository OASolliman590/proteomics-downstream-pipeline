#!/usr/bin/env python3
from __future__ import annotations
import argparse,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; REGISTRY=ROOT/"docs"/"validation"/"baseline-hashes.json"
TRUSTED_REGISTRY_SHA256="a634351e0f2692d65bc6c6786e0b20b79594bde18c48485424b16255a9513a80"
def sha256(path):
 import hashlib
 digest=hashlib.sha256()
 with Path(path).open("rb") as handle:
  for chunk in iter(lambda:handle.read(1024*1024),b""):digest.update(chunk)
 return digest.hexdigest()
def check(registry_path=REGISTRY,root=ROOT,enforce_trust=False):
 registry_path=Path(registry_path); trusted=(registry_path.resolve()==REGISTRY.resolve() and sha256(registry_path)==TRUSTED_REGISTRY_SHA256)
 data=json.loads(registry_path.read_text(encoding="utf-8")); mismatches=[]
 if enforce_trust and not trusted: mismatches.append({"path":str(registry_path),"reason":"untrusted_registry"})
 for item in data.get("files",[]):
  path=Path(root)/item["path"]
  if not path.is_file():mismatches.append({"path":item["path"],"reason":"missing"})
  else:
   observed=sha256(path)
   if observed!=item["sha256"]:mismatches.append({"path":item["path"],"reason":"hash_mismatch","expected":item["sha256"],"observed":observed})
 return not mismatches,mismatches
def main(argv=None):
 parser=argparse.ArgumentParser(); parser.add_argument("--registry",type=Path,default=REGISTRY); args=parser.parse_args(argv); ok,mismatches=check(args.registry,enforce_trust=True); print(json.dumps({"valid":ok,"registry":str(args.registry),"mismatches":mismatches},sort_keys=True)); return 0 if ok else 5
if __name__=="__main__":raise SystemExit(main())
