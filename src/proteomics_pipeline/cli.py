from __future__ import annotations
import argparse,json,sys
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
    # Amendment A-2026-10-01-02: R07 activates the explicit resource preparation command (contract: resources prepare).
    resources=commands.add_parser("resources",help="prepare local hashed resource snapshots (explicit action)"); actions=resources.add_subparsers(dest="action",metavar="ACTION")
    prepare=actions.add_parser("prepare",help="build a local snapshot from local source files"); prepare.add_argument("--manifest",required=True); prepare.add_argument("--output",required=True); prepare.add_argument("--json",action="store_true")
    # Amendment A-2026-10-01-06: R10b activates report and compare (contract: report --run, compare --left --right --output).
    report=commands.add_parser("report",help="render the full offline report of an existing run"); report.add_argument("--run",required=True); report.add_argument("--output"); report.add_argument("--json",action="store_true")
    compare=commands.add_parser("compare",help="compare two verified runs by stable keys"); compare.add_argument("--left",required=True); compare.add_argument("--right",required=True); compare.add_argument("--output",required=True); compare.add_argument("--json",action="store_true")
    # Amendment A-2026-10-01-07: R11 activates resume (contract: resume --run).
    resume=commands.add_parser("resume",help="resume or refresh a run, reusing only verified unchanged stages"); resume.add_argument("--run",required=True); resume.add_argument("--json",action="store_true")
    return parser
# Amendment A-2026-10-01-12: CLI output is always UTF-8, independent of the console/locale encoding (cp1252 on Windows, ASCII under LC_ALL=C).
def _emit(payload,*,json_mode=True):
    text=json.dumps(payload,ensure_ascii=False,sort_keys=True,separators=(",",":") if json_mode else None)+"\n"
    buffer=getattr(sys.stdout,"buffer",None)
    if buffer is None: sys.stdout.write(text); return
    sys.stdout.flush(); buffer.write(text.encode("utf-8")); buffer.flush()
def _error(error,json_mode): _emit({"valid":False,"errors":[error.as_dict()]},json_mode=json_mode); return error.exit_code
def _workflow(args,command):
    """Amendment A-2026-10-01-02: delegate full validate/plan/run to the integration seam."""
    from . import workflow
    try:
        if command=="validate": payload,code=workflow.validate_full(args.config)
        elif command=="plan": payload,code=workflow.plan_command(args.config,args.output)
        else: payload,code=workflow.run_command(args.config,args.output)
    except ProteomicsError as error:return _error(error,args.json)
    _emit(payload,json_mode=args.json); return code
def _validate(args):
    try: config=load_config(args.config)
    except ProteomicsError as error:return _error(error,args.json)
    if not args.schema_only:return _workflow(args,"validate")
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
    from .workflow import verify_run
    try: plan_check=verify_run(root)
    except ProteomicsError as exc:return _error(exc,args.json)
    _emit({"valid":True,"scope":"integrity_and_plan_consistency","scientific_validation":"not_claimed","plan":plan_check,"run_status":value},json_mode=args.json); return 0
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
    if args.command in {"plan","run"}:return _workflow(args,args.command)
    if args.command=="verify":return _verify(args)
    if args.command=="resume":
        from . import workflow
        try: payload,code=workflow.resume_command(args.run)
        except ProteomicsError as error:return _error(error,args.json)
        _emit(payload,json_mode=args.json); return code
    if args.command in {"report","compare"}:
        from . import workflow
        try: payload,code=(workflow.report_command(args.run,args.output) if args.command=="report" else workflow.compare_command(args.left,args.right,args.output))
        except ProteomicsError as error:return _error(error,args.json)
        _emit(payload,json_mode=args.json); return code
    if args.command=="resources":
        if getattr(args,"action",None)!="prepare": parser.error("resources requires the prepare action"); return 2
        from .resources import prepare
        try: payload=prepare(args.manifest,args.output)
        except ProteomicsError as error:return _error(error,args.json)
        _emit(payload,json_mode=args.json); return 0
    parser.error(f"unknown command: {args.command}"); return 2
