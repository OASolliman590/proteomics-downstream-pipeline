from __future__ import annotations
import argparse,json
from pathlib import Path
from . import __version__
from .config import load_config
from .doctor import exit_code as doctor_exit_code, inspect as inspect_doctor
from .errors import ProteomicsError
from .provenance import canonical_config_sha256, sha256_file
from .runtime import validate_run_status
def _parser():
    parser=argparse.ArgumentParser(prog="proteomics",description="Maintained downstream proteomics workflow")
    parser.add_argument("--version",action="version",version=__version__)
    commands=parser.add_subparsers(dest="command",metavar="COMMAND")
    doctor=commands.add_parser("doctor",help="inspect the actual local runtime"); doctor.add_argument("--json",action="store_true"); doctor.add_argument("--config")
    validate=commands.add_parser("validate",help="validate an analysis configuration"); validate.add_argument("--config",required=True); validate.add_argument("--json",action="store_true"); validate.add_argument("--schema-only",action="store_true"); validate.add_argument("--no-schema-only",action="store_false",dest="schema_only"); validate.set_defaults(schema_only=False)
    plan=commands.add_parser("plan",help="freeze a validated analysis plan"); plan.add_argument("--config",required=True); plan.add_argument("--output",required=True); plan.add_argument("--json",action="store_true")
    run=commands.add_parser("run",help="execute an analysis run"); run.add_argument("--config",required=True); run.add_argument("--output",required=True); run.add_argument("--json",action="store_true")
    verify=commands.add_parser("verify",help="verify foundation artifacts"); verify.add_argument("--run",required=True); verify.add_argument("--json",action="store_true")
    return parser
def _emit(payload,*,json_mode=True): print(json.dumps(payload,ensure_ascii=False,sort_keys=True,separators=(",",":") if json_mode else None))
def _error(error,json_mode): _emit({"valid":False,"errors":[error.as_dict()]},json_mode=json_mode); return error.exit_code
def _validate(args):
    try: config=load_config(args.config)
    except ProteomicsError as error:return _error(error,args.json)
    if not args.schema_only:return _unsupported(args,"validate")
    _emit({"valid":True,"schema_only":True,"config":config,"config_hash":canonical_config_sha256(config)},json_mode=args.json); return 0
def _unsupported(args,command): return _error(ProteomicsError("E_CAPABILITY_NOT_IMPLEMENTED",f"{command} is unavailable until a validated command-level handler is implemented",exit_code=3),args.json)
def _verify(args):
    root=Path(args.run).resolve(); status=root/"run_status.json"
    if not status.is_file():return _error(ProteomicsError("E_INTEGRITY","run_status.json is missing",exit_code=5),args.json)
    try:value=json.loads(status.read_text(encoding="utf-8")); validate_run_status(value,root)
    except ValueError as exc:return _error(ProteomicsError("E_INTEGRITY",str(exc),exit_code=5),args.json)
    except ProteomicsError as exc:return _error(exc,args.json)
    except (OSError,json.JSONDecodeError) as exc:return _error(ProteomicsError("E_INTEGRITY",str(exc),exit_code=5),args.json)
    failures=[]
    for artifact in value.get("artifacts",[]):
        path=(root/artifact["relative_path"]).resolve()
        if root not in path.parents and path!=root: failures.append({"code":"E_PATH_ESCAPE","path":artifact["relative_path"]})
        elif not path.is_file() or sha256_file(path)!=artifact["sha256"]: failures.append({"code":"E_ARTIFACT_HASH","path":artifact["relative_path"]})
    if failures:_emit({"valid":False,"scope":"foundation","errors":failures},json_mode=args.json); return 5
    _emit({"valid":True,"scope":"foundation","run_status":value},json_mode=args.json); return 0
def main(argv=None):
    parser=_parser()
    try:args=parser.parse_args(argv)
    except SystemExit as exc:return int(exc.code)
    if args.command is None: parser.print_help(); return 0
    if args.command=="doctor":
        if args.config:
            try: load_config(args.config)
            except ProteomicsError as error: return _error(error,args.json)
        report=inspect_doctor(config_path=args.config); _emit(report,json_mode=args.json); return doctor_exit_code(report)
    if args.command=="validate":return _validate(args)
    if args.command in {"plan","run"}:return _unsupported(args,args.command)
    if args.command=="verify":return _verify(args)
    parser.error(f"unknown command: {args.command}"); return 2
