"""Command-level DAG integration for validate/plan/run (Maintainer amendment A-2026-10-01-02).

This module is the integration seam that the R01 CLI delegates to once
real stage handlers exist.  It executes only capabilities discovered
through the fixed R01 map, freezes the AnalysisPlan before any model fit,
promotes stage outputs atomically, keeps failed outputs in a diagnostic
area, and derives run state strictly from executed stage results.
"""
from __future__ import annotations

import json
import os
import secrets
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

from . import __version__
from .config import _read as read_raw_config, load_config
from .errors import CapabilityError, CollisionError, IntegrityError, ProteomicsError
from .provenance import canonical_json_bytes, content_sha256, sha256_bytes, sha256_file
from .runtime import RunLock, capabilities as discovered_capabilities, execute_capability, stage_result, utc_now, validate_run_status

PLANNING_CAPABILITIES = ("intake", "preprocessing", "design")
REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
R_PACKAGES = ("jsonlite", "openssl", "proteomicsCore", "limma", "statmod", "impute", "vegan", "permute", "svglite")


class PlanRejected(ProteomicsError):
    pass


# --------------------------------------------------------------------------- configuration
def precheck(config_path: str | Path) -> dict:
    """Scientific guards that must win over generic schema messages."""
    raw = read_raw_config(Path(config_path))
    from .intake.hierarchy import scale_transition
    source = raw.get("source_scale"); transform = (raw.get("preprocessing") or {}).get("transform")
    if source and transform:
        scale_transition(source, transform)
    runtime = raw.get("runtime") or {}
    if runtime.get("phase") == 1 and raw.get("score_test", "off") != "off":
        raise PlanRejected("E_PHASE_CAPABILITY", "Phase 1 has no score-testing capability; score_test must be 'off'", "/score_test")
    if (raw.get("post_de") or {}).get("enabled"):   # A-2026-10-01-14: typed post-DE declaration refusals win over schema messages
        from . import post_de
        post_de.precheck(raw)
    return raw


def load(config_path: str | Path) -> dict:
    precheck(config_path)
    return load_config(config_path)


def requested_capabilities(config: dict) -> list[dict]:
    """CapabilityPlan entries with explicit requiredness (independent of installation)."""
    analysis = config["runtime"]["scope"] == "analysis"
    items = [{"capability": "intake", "required": True, "source": "input"},
             {"capability": "preprocessing", "required": True, "source": "preprocessing"}]
    if analysis:
        items.append({"capability": "design", "required": True, "source": "design"})
        for model in config["models"]:
            capability = {"limma": "limma", "deqms": "assay_engines", "proda": "assay_engines"}[model["engine"]]
            items.append({"capability": capability, "required": model["execution_requirement"] == "required", "source": f"model:{model['id']}"})
        if config.get("resources"):
            items.append({"capability": "resources", "required": True, "source": "resources"})
        if config.get("pathways", {}).get("enabled"):
            items.append({"capability": "pathways", "required": True, "source": "pathways"})
        if config.get("response", {}).get("enabled"):
            items.append({"capability": "response", "required": True, "source": "response"})
    multivariate = config.get("multivariate") or {}
    if multivariate.get("enabled"):
        items.append({"capability": "permanova", "required": multivariate.get("execution_requirement", "optional") == "required", "source": "multivariate"})
    post_de_block = config.get("post_de") or {}
    if post_de_block.get("enabled"):   # A-2026-10-01-14 (ADR 0009): declared post-DE modules, then the R14f eligibility report
        from . import post_de
        for module in post_de.requested_modules(config):
            items.append({"capability": post_de.CAPABILITIES[module], "required": post_de.required(config, module), "source": f"post_de.{module}"})
        if post_de.module_impl("eligibility") is not None:   # the R14f eligibility/dependency report, once that packet exists
            items.append({"capability": post_de.ELIGIBILITY_CAPABILITY, "required": False, "source": "post_de"})
    items.append({"capability": "report_stub", "required": False, "source": "report"})
    items.append({"capability": "report_full", "required": False, "source": "report"})
    merged: dict[str, dict] = {}
    for item in items:
        entry = merged.setdefault(item["capability"], {"capability": item["capability"], "required": False, "sources": []})
        entry["required"] = entry["required"] or item["required"]
        entry["sources"].append(item["source"])
    implemented = {c["id"]: c["implemented"] for c in discovered_capabilities()}
    for entry in merged.values():
        entry["implemented"] = bool(implemented.get(entry["capability"], False))
    return list(merged.values())


def environment_inventory() -> dict:
    from .preprocessing_service import rscript_executable
    rscript = rscript_executable()
    info = {"python": {"version": __import__("sys").version.split()[0]}, "proteomics_pipeline": __version__, "rscript": None, "r_packages": {}}
    if not shutil.which(rscript):
        info["rscript"] = {"available": False}
        return info
    code = ("p <- commandArgs(TRUE); v <- lapply(p, function(x) if (requireNamespace(x, quietly=TRUE)) as.character(utils::packageVersion(x)) else NULL); names(v) <- p; "
            "cat(jsonlite::toJSON(list(r=R.version.string, packages=v), auto_unbox=TRUE, null='null'))")
    from .runtime import run_r_code
    result = run_r_code(code, R_PACKAGES, rscript=rscript)
    if result.returncode != 0:
        info["rscript"] = {"available": False, "error": result.stderr.strip()[-500:]}
        return info
    payload = json.loads(result.stdout.strip().splitlines()[-1])
    info["rscript"] = {"available": True, "version": payload["r"]}
    info["r_packages"] = payload["packages"]
    return info


def code_manifest(package_root: Path | None = None, repository_root: Path | None = None) -> dict:
    """Code identity: LF-normalised content hashes (D-42), so a CRLF checkout (for example DESCRIPTION or the HTML template
    under `* text=auto` on Windows) has the same identity as an LF checkout."""
    files = []
    package_root = Path(package_root) if package_root is not None else Path(__file__).resolve().parent
    repository_root = Path(repository_root) if repository_root is not None else REPOSITORY_ROOT
    for path in sorted(package_root.rglob("*")):
        if path.is_file() and "__pycache__" not in path.parts and path.suffix in (".py", ".json", ".html", ".typed"):
            files.append((f"src/proteomics_pipeline/{path.relative_to(package_root).as_posix()}", content_sha256(path)))
    r_root = repository_root / "r" / "proteomicsCore"
    if r_root.is_dir():
        for path in sorted((r_root / "R").glob("*.R")) + [r_root / "DESCRIPTION", r_root / "NAMESPACE", repository_root / "scripts" / "maintained" / "run_stage.R"]:
            if path.is_file():
                files.append((path.relative_to(repository_root).as_posix(), content_sha256(path)))
    lines = "".join(f"{digest}  {name}\n" for name, digest in sorted(files))
    return {"sha256": sha256_bytes(lines.encode("utf-8")), "n_files": len(files), "r_source_present": r_root.is_dir()}


# --------------------------------------------------------------------------- stage execution
class Ledger:
    """Records executed stages; run state is derived only from these records."""

    def __init__(self, root: Path, run_id: str, *, reuse: bool = False):
        self.root, self.run_id, self.reuse = root, run_id, reuse
        self.stages: list[dict] = []
        self.warnings: list[dict] = []
        self.reused: list[str] = []
        self.invalidated: list[str] = []
        self._code = self._env = None

    def _fingerprint(self, request: dict) -> str:
        from . import cache
        if self._code is None:
            self._code = code_manifest()["sha256"]
            self._env = sha256_bytes(canonical_json_bytes(environment_inventory()))
        return cache.fingerprint(request, code_sha256=self._code, environment_sha256=self._env)

    def run_stage(self, capability: str, stage_id: str, required: bool, build: Callable[[Path], dict], destination: Path) -> dict:
        from . import cache
        temp = self.root / f".stage-{stage_id}-{secrets.token_hex(4)}"
        request = build(temp)
        fp = self._fingerprint(request)
        if destination.exists():
            if self.reuse and cache.reusable(destination, fp):
                result = json.loads((destination / "stage-result.json").read_text(encoding="utf-8"))
                self.reused.append(stage_id)
                self.warnings.extend(result.get("warnings", []))
                self.stages.append({"stage_id": stage_id, "capability": capability, "required": required, "state": "COMPLETED", "reason_code": None,
                                    "result_path": (destination / "stage-result.json").relative_to(self.root).as_posix(), "exit_code": 0, "message": result.get("message", "")})
                return result
            if not self.reuse:
                raise CollisionError(f"stage destination already exists: {destination}")
            self.invalidated.append(stage_id)
            cache.retire(destination, self.root, f"stale-{stage_id}")
        try:
            result = execute_capability(capability, request)
        except CapabilityError as error:
            result = stage_result(self.run_id, stage_id, capability, "NOT_RUN", plan_hash=request["plan_hash"], exit_code=3, reason_code=error.code, message=error.message)
        except ProteomicsError as error:
            result = stage_result(self.run_id, stage_id, capability, "FAILED", plan_hash=request["plan_hash"], exit_code=error.exit_code, reason_code=error.code, message=error.message)
        result_file = temp / "stage-result.json"
        if temp.is_dir() and not (temp / "stage-request.json").is_file():
            (temp / "stage-request.json").write_bytes(canonical_json_bytes(request))
        if temp.is_dir() and not result_file.is_file():
            temp.joinpath("stage-result.json").write_bytes(json.dumps(result, ensure_ascii=False, indent=2).encode("utf-8"))
        if result["state"] == "COMPLETED":
            if destination.exists():
                raise CollisionError(f"stage destination already exists: {destination}")
            destination.parent.mkdir(parents=True, exist_ok=True)
            os.replace(temp, destination)
            cache.record(destination, fp)
            result_path = (destination / "stage-result.json").relative_to(self.root).as_posix()
        else:
            failed = self.root / "logs" / f"{stage_id}-failed"
            result_path = None
            if temp.is_dir():
                failed.parent.mkdir(parents=True, exist_ok=True)
                os.replace(temp, failed)
                if not (failed / "stage-result.json").is_file():
                    (failed / "orchestrator-stage-result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
                else:
                    result_path = (failed / "stage-result.json").relative_to(self.root).as_posix()
        for warning in result.get("warnings", []):
            self.warnings.append(warning)
        self.stages.append({"stage_id": stage_id, "capability": capability, "required": required, "state": result["state"],
                            "reason_code": result.get("reason_code"), "result_path": result_path, "exit_code": result.get("exit_code"),
                            "message": result.get("message", "")})
        return result

    def virtual(self, stage_id: str, capability: str, required: bool, state: str, reason_code: str | None, message: str, exit_code: int | None = None):
        self.stages.append({"stage_id": stage_id, "capability": capability, "required": required, "state": state, "reason_code": reason_code,
                            "result_path": None, "exit_code": exit_code, "message": message})


def _python_request(run_id: str, stage_id: str, capability: str, temp: Path, *, config_path: str | None, inputs: list[dict], parameters: dict, plan_hash=None, seed=0, threads=1) -> dict:
    return {"schema_version": "1.2.0", "run_id": run_id, "stage_id": stage_id, "capability": capability, "plan_hash": plan_hash, "inputs": inputs,
            "output_temp_dir": str(temp.resolve()), "config_path": config_path, "parameters": parameters,
            "rng": {"seed": int(seed), "kind": "L'Ecuyer-CMRG", "threads": int(threads)}}


def _required_unavailable(requested: list[dict]) -> None:
    for item in requested:
        if item["required"] and not item["implemented"]:
            raise CapabilityError("E_CAPABILITY_NOT_IMPLEMENTED", f"required capability {item['capability']!r} (requested by {', '.join(item['sources'])}) has no validated handler in this release")


def _stage_failure(result: dict) -> ProteomicsError:
    code = result.get("reason_code") or "E_STAGE_FAILED"
    return ProteomicsError(code, f"stage {result['stage_id']} {result['state']}: {result.get('message', '')}".strip(), exit_code=result.get("exit_code") or 4)


def plan_into(config_path: str | Path, config: dict, *, plan_root: Path, base: Path, ledger: Ledger, requested: list[dict]) -> dict:
    """Run intake, preprocessing and design, then freeze the AnalysisPlan (no fit)."""
    from . import design_service, preprocessing_service
    from .intake import service as intake
    config_path = Path(config_path).resolve()
    run_id = ledger.run_id
    seed, threads = config["runtime"]["seed"], config["runtime"]["threads"]
    if config.get("resources"):
        from . import resources as resource_service
        resource_service.verify_resources(config, config_path.parent)   # every snapshot verified before analysis (SM14)
    if any(m["engine"] in ("deqms", "proda") for m in config["models"]):
        from . import assay_service, design_service as _ds
        # adapter eligibility depends only on declarations, so it is decided before any data stage (SM10)
        for row in _ds.eligibility_table(config, mask_state="known", prior_imputation=config["input"]["prior_imputation"]):
            if row["engine"] != "limma" and not row["eligible"] and row["execution_requirement"] == "required":
                raise PlanRejected(row["reason_code"], f"required model {row['model_id']!r} ({row['engine']}) is scientifically ineligible: {row['reason_code']}", f"/models[{row['model_id']}]")
        assay_service.plan_checks(config, config_path.parent)
    sources = intake.source_inputs(config, config_path.parent)
    result = ledger.run_stage("intake", "intake", True, lambda temp: _python_request(run_id, "intake", "intake", temp, config_path=str(config_path), inputs=sources,
                              parameters={"source_scale": config["source_scale"]}, seed=seed, threads=threads), base / "inputs")
    if result["state"] != "COMPLETED":
        raise _stage_failure(result)
    result = ledger.run_stage("preprocessing", "preprocessing", True, lambda temp: preprocessing_service.build_request(config, config_path=config_path, bundle_dir=base / "inputs",
                              run_id=run_id, output_temp_dir=temp), base / "preprocessing")
    if result["state"] != "COMPLETED":
        raise _stage_failure(result)
    stage_dirs = {"intake": base / "inputs", "preprocessing": base / "preprocessing"}
    analysis = config["runtime"]["scope"] == "analysis"
    if analysis:
        manifest = json.loads((base / "inputs" / "manifest.json").read_text(encoding="utf-8"))
        scale = manifest["scale"]
        if not scale["inference_eligible"]:
            raise PlanRejected(scale["inference_reason_code"], f"abundance inference is not eligible on scale {scale['output_scale']!r}; QC-only scope remains available", "/source_scale")
        eligibility = design_service.eligibility_table(config, mask_state=manifest["original_observed_mask"]["state"], prior_imputation=config["input"]["prior_imputation"])
        for row in eligibility:
            if not row["eligible"] and row["execution_requirement"] == "required":
                raise PlanRejected(row["reason_code"], f"required model {row['model_id']!r} ({row['engine']}) is scientifically ineligible: {row['reason_code']}", f"/models[{row['model_id']}]")
        if config.get("omnibus_tests"):
            for index, test in enumerate(config["omnibus_tests"]):
                if test.get("reduced_design_id") is not None:
                    raise PlanRejected("E_OMNIBUS_REDUCED_UNSUPPORTED", "reduced-design omnibus tests are not exposed in this release; declare coefficient_names only", f"/omnibus_tests/{index}/reduced_design_id")
        try:
            result = ledger.run_stage("design", "design", True, lambda temp: design_service.build_request(config, config_path=config_path, preprocessing_dir=base / "preprocessing",
                                      run_id=run_id, output_temp_dir=temp), base / "designs")
        except ProteomicsError as error:
            raise PlanRejected(error.code, error.message, error.pointer) from error
        if result["state"] != "COMPLETED":
            raise _stage_failure(result)
        design_service.interpret(base / "designs", config)
        stage_dirs["design"] = base / "designs"
        if (config.get("pathways") or {}).get("enabled"):
            from . import pathway_service
            pathway_service.plan_checks(config)
        if (config.get("response") or {}).get("enabled"):
            from . import response_service
            design_request = json.loads((base / "designs" / "stage-request.json").read_text(encoding="utf-8"))
            response_service.plan_checks(config, design_request["parameters"]["contrasts"])
        import csv
        with (base / "designs" / "estimability.tsv").open(encoding="utf-8", newline="") as handle:
            estimability = list(csv.DictReader(handle, delimiter="\t"))
        primary = next(m for m in config["models"] if m["role"] == "primary")
        if not any(r["model_id"] == primary["id"] and r["eligibility"] == "eligible" for r in estimability):
            raise PlanRejected("E_NO_ELIGIBLE_FEATURES", "the primary model has no eligible protein features; zero eligible tests is not a zero-discovery analysis", "/models")
    permanova_refusals: list[tuple[str, str]] = []
    if (config.get("multivariate") or {}).get("enabled"):
        from . import permanova_service
        observations, feature_ids = permanova_service.read_observations(base / "preprocessing")
        permanova_refusals = permanova_service.plan_checks(config, observations, feature_ids)
        if permanova_refusals and multivariate_required(config):
            code, message = permanova_refusals[0]
            raise PlanRejected(code, f"required PERMANOVA is scientifically ineligible: {message}", "/multivariate")
    post_de_decisions = post_de_plan(config, config_path, base, analysis)
    environment = environment_inventory()
    code = code_manifest()
    capability_plan = []
    for item in requested:
        entry = dict(item, scientific_eligibility="eligible", reason_code=None)
        if item["capability"] == "permanova" and permanova_refusals:
            entry.update(scientific_eligibility="inapplicable", reason_code=permanova_refusals[0][0], reason=permanova_refusals[0][1])
        decision = post_de_decisions.get(item["capability"])
        if decision is not None:   # SM41: the planner records each post-DE module's eligibility before any computation
            entry["post_de"] = {k: v for k, v in decision.items() if k not in ("resolved",)}
            if decision["state"] == "INAPPLICABLE":
                entry.update(scientific_eligibility="inapplicable", reason_code=decision["reason_code"], reason=decision.get("reason"))
        capability_plan.append(entry)
    plan = design_service.build_plan(config, plan_root=plan_root, stage_dirs=stage_dirs, sources=sources, environment=environment,
                                     code_sha256=code["sha256"], capability_plan=capability_plan)
    return plan


# --------------------------------------------------------------------------- commands
def _new_run_id() -> str:
    return f"run-{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}-{secrets.token_hex(3)}"


def validate_full(config_path: str | Path) -> tuple[dict, int]:
    config = load(config_path)
    requested = requested_capabilities(config)
    _required_unavailable(requested)
    with tempfile.TemporaryDirectory(prefix="proteomics-validate-") as scratch:
        root = Path(scratch)
        ledger = Ledger(root, _new_run_id())
        plan = plan_into(config_path, config, plan_root=root, base=root, ledger=ledger, requested=requested)
    return ({"valid": True, "scope": "full", "plan_hash_preview": plan["plan_hash"], "stages": [{k: s[k] for k in ("stage_id", "state", "reason_code")} for s in ledger.stages],
             "note": "Full validation executed intake, preprocessing and design checks without fitting; the preview plan was discarded."}, 0)


def plan_command(config_path: str | Path, output: str | Path) -> tuple[dict, int]:
    config = load(config_path)
    requested = requested_capabilities(config)
    _required_unavailable(requested)
    output = Path(output).resolve()
    artifacts = output.parent / "plan-artifacts"
    if output.exists() or artifacts.exists():
        raise CollisionError(f"plan output or plan-artifacts already exists next to {output.name}")
    artifacts.mkdir(parents=True)
    ledger = Ledger(output.parent, _new_run_id())
    try:
        plan = plan_into(config_path, config, plan_root=output.parent, base=artifacts, ledger=ledger, requested=requested)
    except BaseException:
        # a refused plan leaves nothing behind: plan-artifacts was created by this call
        for leftover in output.parent.glob(".stage-*"):
            shutil.rmtree(leftover, ignore_errors=True)
        shutil.rmtree(artifacts, ignore_errors=True)
        shutil.rmtree(output.parent / "logs", ignore_errors=True) if not any((output.parent / "logs").glob("*")) else None
        raise
    from .design_service import write_plan
    write_plan(plan, output)
    return ({"plan": str(output.name), "plan_hash": plan["plan_hash"], "scope": plan["scope"], "artifacts": len(plan["artifacts"]),
             "fits_performed": 0}, 0)


def run_command(config_path: str | Path, output: str | Path) -> tuple[dict, int]:
    config = load(config_path)
    requested = requested_capabilities(config)
    _required_unavailable(requested)
    root = Path(output).resolve()
    if root.exists() and any(root.iterdir()):
        raise CollisionError(f"run output already exists and overwrite=false: {root.name}")
    root.mkdir(parents=True, exist_ok=True)
    run_id = _new_run_id()
    started = utc_now()
    ledger = Ledger(root, run_id)
    with RunLock(root / ".run.lock"):
        (root / "config.resolved.json").write_bytes(json.dumps(config, indent=2, ensure_ascii=False, sort_keys=True).encode("utf-8"))
        (root / "provenance").mkdir(exist_ok=True)
        (root / "provenance" / "config_source.json").write_text(json.dumps({"config_path": str(Path(config_path).resolve())}), encoding="utf-8")
        plan = None
        failure: ProteomicsError | None = None
        try:
            plan = plan_into(config_path, config, plan_root=root, base=root, ledger=ledger, requested=requested)
            from .design_service import write_plan
            write_plan(plan, root / "plan.json")
        except ProteomicsError as error:
            failure = error
            if not ledger.stages or ledger.stages[-1]["state"] == "COMPLETED":
                ledger.virtual("plan", "design", True, "FAILED", error.code, error.message, error.exit_code)
        plan_hash = plan["plan_hash"] if plan else None
        if failure is None:
            _run_analysis(config, plan, root, ledger, Path(config_path).resolve().parent)
        for item in requested:
            if not item["implemented"] and not item["required"]:
                ledger.virtual(item["capability"], item["capability"], False, "NOT_RUN", "E_CAPABILITY_NOT_IMPLEMENTED", "optional requested capability has no handler in this release", 3)
        _write_provenance(root, config, plan)
        status = _finalize(root, config, ledger, run_id, plan_hash, started, report=True)
    payload = {"run_id": run_id, "state": status["state"], "exit_code": status["exit_code"], "plan_hash": plan_hash,
               "stages": [{k: s[k] for k in ("stage_id", "state", "reason_code")} for s in ledger.stages], "report": "report/index.html" if (root / "report" / "index.html").is_file() else None}
    if failure is not None:
        payload["error"] = failure.as_dict()
    return payload, status["exit_code"]


def _run_analysis(config: dict, plan: dict, root: Path, ledger: Ledger, config_dir: Path) -> None:
    from . import inference_service
    analysis = config["runtime"]["scope"] == "analysis"
    eligibility = {row["model_id"]: row for row in plan["engine_eligibility"]}
    if analysis:
        for row in plan["engine_eligibility"]:
            if not row["eligible"]:
                ledger.virtual(f"model.{row['model_id']}", "limma" if row["engine"] == "limma" else "assay_engines", False, "INAPPLICABLE", row["reason_code"],
                               f"optional model {row['model_id']} is scientifically ineligible", 0)
        limma_models = [m for m in config["models"] if m["engine"] == "limma" and eligibility[m["id"]]["eligible"]]
        if limma_models:
            result = ledger.run_stage("limma", "limma", True, lambda temp: inference_service.build_request(plan, plan_path=root / "plan.json", config=config, run_id=ledger.run_id,
                                      output_temp_dir=temp), root / "dea")
            if result["state"] == "COMPLETED":
                for status in inference_service.model_status(root / "dea"):
                    if status["state"] != "COMPLETED":
                        ledger.virtual(f"model.{status['model_id']}", "limma", bool(status.get("required")), status["state"], status.get("reason_code"), status.get("message", ""), 4)
        assay = [m for m in config["models"] if m["engine"] in ("deqms", "proda") and eligibility[m["id"]]["eligible"]]
        if assay:
            _run_assay_engines(config, plan, root, ledger, config_dir, assay)
        if config.get("resources"):
            _run_resources(config, plan, root, ledger, config_dir)
        if (config.get("response") or {}).get("enabled"):
            _run_response(config, plan, root, ledger, config_dir)
        if (config.get("pathways") or {}).get("enabled"):
            _run_pathways(config, plan, root, ledger)
    multivariate = config.get("multivariate") or {}
    if multivariate.get("enabled"):
        _run_permanova(config, plan, root, ledger)
    if analysis and (config.get("post_de") or {}).get("enabled"):
        _run_post_de(config, plan, root, ledger, config_dir)


def _run_assay_engines(config: dict, plan: dict, root: Path, ledger: Ledger, config_dir: Path, models: list[dict]) -> None:
    from . import assay_service
    required = any(m["execution_requirement"] == "required" for m in models)
    if not (root / "dea" / "stage-result.json").is_file():
        ledger.virtual("assay_engines", "assay_engines", required, "NOT_RUN", "E_PREREQUISITE_FAILED", "alternative engines are compared with the completed primary model", 4)
        return
    result = ledger.run_stage("assay_engines", "assay_engines", required, lambda temp: assay_service.build_request(plan, plan_path=root / "plan.json", config=config, config_dir=config_dir,
                              run_id=ledger.run_id, output_temp_dir=temp, dea_dir=root / "dea"), root / "assay_engines")
    if result["state"] == "COMPLETED":
        for status in assay_service.model_status(root / "assay_engines"):
            if status["state"] != "COMPLETED":
                exit_code = 3 if status["state"] == "NOT_RUN" else (0 if status["state"] == "INAPPLICABLE" else 4)
                ledger.virtual(f"model.{status['model_id']}", "assay_engines", bool(status.get("required")), status["state"], status.get("reason_code"), status.get("message", ""), exit_code)


def _run_resources(config: dict, plan: dict, root: Path, ledger: Ledger, config_dir: Path) -> None:
    from . import resources as resource_service
    ledger.run_stage("resources", "resources", True, lambda temp: resource_service.build_request(plan, plan_path=root / "plan.json", config=config, config_dir=config_dir,
                     run_id=ledger.run_id, output_temp_dir=temp), root / "resources")


def _run_pathways(config: dict, plan: dict, root: Path, ledger: Ledger) -> None:
    try:
        from . import pathway_service
    except ImportError:
        return
    if not (root / "resources" / "stage-result.json").is_file() or not (root / "dea" / "stage-result.json").is_file():
        ledger.virtual("pathways", "pathways", True, "NOT_RUN", "E_PREREQUISITE_FAILED", "pathways need completed resources and limma stages", 4)
        return
    ledger.run_stage("pathways", "pathways", True, lambda temp: pathway_service.build_request(plan, plan_path=root / "plan.json", config=config, run_id=ledger.run_id,
                     output_temp_dir=temp, dea_dir=root / "dea", resources_dir=root / "resources"), root / "pathways")


def _run_response(config: dict, plan: dict, root: Path, ledger: Ledger, config_dir: Path) -> None:
    from . import response_service
    if not (root / "dea" / "stage-result.json").is_file():
        ledger.virtual("response", "response", True, "NOT_RUN", "E_PREREQUISITE_FAILED", "response needs the completed limma stage", 4)
        return
    ledger.run_stage("response", "response", True, lambda temp: response_service.build_request(plan, plan_path=root / "plan.json", config=config, run_id=ledger.run_id,
                     output_temp_dir=temp, dea_dir=root / "dea", config_dir=config_dir), root / "response")


def _run_permanova(config: dict, plan: dict, root: Path, ledger: Ledger) -> None:
    required = multivariate_required(config)
    entry = next((c for c in plan["capability_plan"] if c["capability"] == "permanova"), None)
    if entry and entry.get("scientific_eligibility") == "inapplicable":
        ledger.virtual("permanova", "permanova", required, "INAPPLICABLE", entry["reason_code"], entry.get("reason", ""), 0)
        return
    try:
        from . import permanova_service
    except ImportError:
        return
    if not any(c["id"] == "permanova" and c["implemented"] for c in discovered_capabilities()):
        return
    dea_ready = (root / "dea" / "stage-result.json").is_file()
    ledger.run_stage("permanova", "permanova", required, lambda temp: permanova_service.build_request(plan, plan_path=root / "plan.json", config=config, run_id=ledger.run_id,
                     output_temp_dir=temp, dea_dir=(root / "dea") if dea_ready else None), root / "permanova")


def post_de_plan(config: dict, config_path: Path, base: Path, analysis: bool) -> dict[str, dict]:
    """Plan-time eligibility of every declared post-DE module (A-2026-10-01-14, SM41); required ineligible modules reject the plan."""
    if not (config.get("post_de") or {}).get("enabled"):
        return {}
    from . import post_de
    if not analysis:
        raise PlanRejected("E_PHASE_CAPABILITY", "post-differential analysis needs runtime.scope = analysis", "/post_de")
    import csv
    with (base / "preprocessing" / "primary" / "observations.tsv").open(encoding="utf-8", newline="") as handle:
        observations = list(csv.DictReader(handle, delimiter="\t"))
    with (base / "preprocessing" / "primary" / "features.tsv").open(encoding="utf-8", newline="") as handle:
        features = list(csv.DictReader(handle, delimiter="\t"))
    design_request = json.loads((base / "designs" / "stage-request.json").read_text(encoding="utf-8"))
    context = {"observations": observations, "features": features, "design_request": design_request, "config_dir": Path(config_path).parent,
               "preprocessing_dir": base / "preprocessing", "design_dir": base / "designs"}
    decisions = post_de.plan_checks(config, context)
    out = {}
    for module, decision in decisions.items():
        if decision["state"] == "INAPPLICABLE" and post_de.required(config, module):
            raise PlanRejected(decision["reason_code"], f"required post-DE module {module!r} is scientifically ineligible: {decision.get('reason')}", f"/post_de/{module}")
        out[post_de.CAPABILITIES[module]] = decision
    return out


def _run_post_de(config: dict, plan: dict, root: Path, ledger: Ledger, config_dir: Path) -> None:
    """Run the declared post-DE modules in dispatch order (R14a-R14e), then the R14f eligibility/dependency stage."""
    from . import post_de
    implemented = {c["id"] for c in discovered_capabilities() if c["implemented"]}
    entries = {c["capability"]: c for c in plan["capability_plan"]}
    for module in post_de.requested_modules(config):
        capability = post_de.CAPABILITIES[module]
        required = post_de.required(config, module)
        entry = entries.get(capability, {})
        if entry.get("scientific_eligibility") == "inapplicable":
            ledger.virtual(capability, capability, required, "INAPPLICABLE", entry["reason_code"], entry.get("reason", ""), 0)
            continue
        if capability not in implemented:
            continue   # recorded NOT_RUN (E_CAPABILITY_NOT_IMPLEMENTED) by the caller's optional-capability pass; required ones never reach here
        impl = post_de.module_impl(module)
        missing = [p for p in getattr(impl, "PREREQUISITES", ("dea",)) if not (root / p / "stage-result.json").is_file()]
        if missing:
            ledger.virtual(capability, capability, required, "NOT_RUN", "E_PREREQUISITE_FAILED", f"post-DE {module} needs completed {', '.join(missing)}", 4)
            continue
        decision = entry.get("post_de", {})
        try:
            ledger.run_stage(capability, capability, required, lambda temp, impl=impl, decision=decision: impl.build_request(
                plan, plan_path=root / "plan.json", config=config, run_id=ledger.run_id, output_temp_dir=temp, root=root, decision=decision, config_dir=config_dir),
                root / "post_de" / module)
        except ProteomicsError as error:   # a request that cannot be built is a typed stage failure, never a crash
            ledger.virtual(capability, capability, required, "FAILED", error.code, error.message, error.exit_code)
    if post_de.ELIGIBILITY_CAPABILITY in implemented:
        from .post_de import eligibility
        snapshot = [dict(s) for s in ledger.stages]
        ledger.run_stage(post_de.ELIGIBILITY_CAPABILITY, post_de.ELIGIBILITY_CAPABILITY, False, lambda temp: _python_request(
            ledger.run_id, post_de.ELIGIBILITY_CAPABILITY, post_de.ELIGIBILITY_CAPABILITY, temp, config_path=None, inputs=[],
            parameters={"run_root": str(root), "stages": snapshot}, plan_hash=plan["plan_hash"]), root / "post_de" / "eligibility")


def multivariate_required(config: dict) -> bool:
    return (config.get("multivariate") or {}).get("execution_requirement", "optional") == "required"


def _write_provenance(root: Path, config: dict, plan: dict | None) -> None:
    provenance = root / "provenance"
    provenance.mkdir(exist_ok=True)
    (provenance / "environment.json").write_text(json.dumps(plan["environment"] if plan else environment_inventory(), indent=2, sort_keys=True), encoding="utf-8")
    (provenance / "code_manifest.json").write_text(json.dumps(code_manifest(), indent=2, sort_keys=True), encoding="utf-8")
    sources = plan["sources"] if plan else []
    (provenance / "input_hashes.json").write_text(json.dumps(sources, indent=2, sort_keys=True), encoding="utf-8")


def _derive_state(stages: list[dict]) -> tuple[str, int, str | None]:
    bad = {"NOT_RUN", "FAILED", "CANCELLED"}
    required_bad = [s for s in stages if s["required"] and s["state"] in bad]
    optional_bad = [s for s in stages if not s["required"] and s["state"] in bad]
    if required_bad:
        first = required_bad[0]
        return "FAILED", int(first.get("exit_code") or 4) or 4, first["reason_code"] or "E_STAGE_FAILED"
    if optional_bad:
        first = optional_bad[0]
        code = int(first.get("exit_code") or 3)
        return "PARTIAL", code if code in (3, 4) else 3, first["reason_code"] or "E_OPTIONAL_STAGE_FAILED"
    return "COMPLETED", 0, None


def _finalize(root: Path, config: dict, ledger: Ledger, run_id: str, plan_hash: str | None, started: str, *, report: bool) -> dict:
    (root / "warnings.json").write_bytes(json.dumps(ledger.warnings, indent=2, ensure_ascii=False).encode("utf-8"))
    state, exit_code, reason = _derive_state(ledger.stages)
    implemented = [c["id"] for c in discovered_capabilities() if c["implemented"]]
    if report and any(c["id"] == "report_stub" and c["implemented"] for c in discovered_capabilities()):
        from .reporting import stub
        snapshot = {"run_id": run_id, "plan_hash": plan_hash, "state": state, "exit_code": exit_code, "reason_code": reason,
                    "stages": [dict(s) for s in ledger.stages], "requested_phase": config["runtime"]["phase"], "implemented_capabilities": implemented}
        ledger.run_stage("report_stub", "report", False, lambda temp: _python_request(run_id, "report", "report_stub", temp, config_path=None, inputs=[],
                         parameters={"run_root": str(root), "run_snapshot": snapshot}, plan_hash=plan_hash), root / "report")
        if any(c["id"] == "report_full" and c["implemented"] for c in discovered_capabilities()):
            ledger.run_stage("report_full", "report_full", False, lambda temp: _python_request(run_id, "report_full", "report_full", temp, config_path=None, inputs=[],
                             parameters={"run_root": str(root), "run_snapshot": snapshot, "figure_formats": list(config["report"].get("figure_formats", []))}, plan_hash=plan_hash), root / "report-full")
        state, exit_code, reason = _derive_state(ledger.stages)
        if (root / "report" / "methods.md").is_file():
            shutil.copyfile(root / "report" / "methods.md", root / "methods.md")
    artifacts = []
    for relative, kind in (("config.resolved.json", "ResolvedConfig"), ("plan.json", "AnalysisPlan"), ("warnings.json", "Warnings"), ("methods.md", "MethodsSummary"),
                           ("provenance/environment.json", "Provenance"), ("provenance/code_manifest.json", "Provenance"), ("provenance/input_hashes.json", "Provenance")):
        path = root / relative
        if path.is_file():
            artifacts.append({"artifact_id": relative.replace("/", "_").replace(".", "_"), "relative_path": relative, "sha256": sha256_file(path), "result_type": kind})
    status = {"schema_version": "1.2.0", "run_id": run_id, "plan_hash": plan_hash, "requested_phase": config["runtime"]["phase"],
              "implemented_capabilities": implemented, "state": state, "reason_code": reason, "started_at": started, "finished_at": utc_now(),
              "exit_code": exit_code, "stages": [{k: s[k] for k in ("stage_id", "capability", "required", "state", "reason_code", "result_path")} for s in ledger.stages],
              "artifacts": artifacts, "warnings_file": "warnings.json"}
    validate_run_status(status, root)
    (root / "run_status.json").write_bytes(json.dumps(status, indent=2, ensure_ascii=False).encode("utf-8"))
    return status


def verify_run(root: str | Path) -> dict:
    """Plan-level consistency on top of the R01 manifest verification."""
    from .design_service import verify_plan
    root = Path(root)
    status = json.loads((root / "run_status.json").read_text(encoding="utf-8"))
    if (root / "plan.json").is_file():
        plan = verify_plan(root / "plan.json", expected_hash=status["plan_hash"])
        return {"plan_verified": True, "plan_hash": plan["plan_hash"], "planned_artifacts": len(plan["artifacts"])}
    if status["plan_hash"] is not None:
        raise IntegrityError("run_status names a plan hash but plan.json is missing", code="E_PLAN_REQUIRED")
    return {"plan_verified": False, "reason": "no plan was frozen for this run"}


def report_command(run: str | Path, output: str | Path | None = None) -> tuple[dict, int]:
    """Re-render the full report of an existing verified run without modifying it (R10b)."""
    from .reporting import full
    root = Path(run).resolve()
    verify_run(root)
    status = json.loads((root / "run_status.json").read_text(encoding="utf-8"))
    validate_run_status(status, root)
    config = json.loads((root / "config.resolved.json").read_text(encoding="utf-8"))
    out = Path(output).resolve() if output else root / "reports" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    if out.exists():
        raise CollisionError(f"report output already exists: {out.name}")
    out.mkdir(parents=True)
    snapshot = {k: status[k] for k in ("run_id", "plan_hash", "state", "exit_code", "reason_code", "requested_phase", "implemented_capabilities")}
    snapshot["stages"] = [dict(s, message="") for s in status["stages"]]
    files = full.build_report(root, snapshot, out, list(config["report"].get("figure_formats", [])))
    return {"report": str(out / "index.html"), "files": len(files), "run_state": status["state"]}, 0


def compare_command(left: str | Path, right: str | Path, output: str | Path) -> tuple[dict, int]:
    from .reporting.comparison import compare
    summary = compare(left, right, output)
    return {"comparison": str(Path(output).resolve() / "index.html"), "same_plan": summary["same_plan"], "config_differences": len(summary["config_differences"])}, 0


def resume_command(run: str | Path) -> tuple[dict, int]:
    """Resume or refresh a run in place, reusing only verified stages with identical fingerprints (R11, V103)."""
    from . import cache
    from .design_service import write_plan
    root = Path(run).resolve()
    status_path = root / "run_status.json"
    source = root / "provenance" / "config_source.json"
    if not source.is_file():
        raise ProteomicsError("E_RESUME_UNSUPPORTED", "this run has no recorded configuration source; start a new run", exit_code=2)
    config_path = Path(json.loads(source.read_text(encoding="utf-8"))["config_path"])
    config = load(config_path)
    requested = requested_capabilities(config)
    _required_unavailable(requested)
    run_id = json.loads(status_path.read_text(encoding="utf-8"))["run_id"] if status_path.is_file() else _new_run_id()
    started = utc_now()
    ledger = Ledger(root, run_id, reuse=True)
    with RunLock(root / ".run.lock"):
        abandoned = cache.abandon_temporaries(root)
        plan = None; failure = None
        try:
            plan = plan_into(config_path, config, plan_root=root, base=root, ledger=ledger, requested=requested)
            existing = root / "plan.json"
            if existing.is_file() and json.loads(existing.read_text(encoding="utf-8")).get("plan_hash") != plan["plan_hash"]:
                cache.retire(existing, root, "stale-plan"); ledger.invalidated.append("plan")
            if not existing.is_file():
                write_plan(plan, existing)
        except ProteomicsError as error:
            failure = error
            if not ledger.stages or ledger.stages[-1]["state"] == "COMPLETED":
                ledger.virtual("plan", "design", True, "FAILED", error.code, error.message, error.exit_code)
        plan_hash = plan["plan_hash"] if plan else None
        if failure is None:
            _run_analysis(config, plan, root, ledger, config_path.parent)
        for item in requested:
            if not item["implemented"] and not item["required"]:
                ledger.virtual(item["capability"], item["capability"], False, "NOT_RUN", "E_CAPABILITY_NOT_IMPLEMENTED", "optional requested capability has no handler in this release", 3)
        for name in ("provenance/environment.json", "provenance/code_manifest.json", "provenance/input_hashes.json"):
            if (root / name).exists():
                (root / name).unlink()
        _write_provenance(root, config, plan)
        (root / "provenance" / "config_source.json").write_text(json.dumps({"config_path": str(config_path)}), encoding="utf-8")
        status = _finalize(root, config, ledger, run_id, plan_hash, started, report=True)
    payload = {"run_id": run_id, "state": status["state"], "exit_code": status["exit_code"], "plan_hash": plan_hash, "reused": ledger.reused,
               "invalidated": ledger.invalidated, "abandoned_temporaries": abandoned}
    if failure is not None:
        payload["error"] = failure.as_dict()
    return payload, status["exit_code"]
