# CLI and artifact contract v1.2.0

**Implementation target, not current functionality.** Today the repository has the recovery scripts, not these entry points. R01 creates the command/parser/runtime foundation; no statistical stage is claimed until its owning packet is accepted.

## Entry points and Phase 1 commands

`proteomics` and `python -m proteomics_pipeline` call the same main function, exit codes and services. A global `--help`/`--version` needs neither R nor input data.

| Command | Contract | R01 behavior / later owner |
|---|---|---|
| `proteomics doctor --json [--config analysis.json]` | Read-only actual binary/package/resource inventory. Without config requires the foundation R I/O capability; with config inspects all requested requirements. Missing required capability exits 3. | Real R01. |
| `proteomics validate --config analysis.json --json [--schema-only]` | Apply the four declared defaults and structural checks. Default full validation also checks data/identity/design and scientific eligibility without fitting. Errors contain code, JSON pointer and message. | Schema-only is real R01. Full validation returns NOT_RUN/E_CAPABILITY_NOT_IMPLEMENTED until R02/R04 handlers exist; it never reports a complete scientific validation from schema alone. |
| `proteomics plan --config analysis.json --output plan.json` | Perform read-only intake/transform/design calculations, save validated design/contrast/inclusion artifacts and freeze the complete plan before any model fit. Refuse collision or unresolved required dependencies. | R01 parser/guards only; no false AnalysisPlan until R02–R04 are accepted. |
| `proteomics run --config analysis.json --output runs/example` | Use the same validation/planning services, freeze first, then run requested accepted capabilities. R10a integrates the Phase 1 DAG and report; required failure stops dependants and produces an honest partial report. | Unsupported stages exit 3/NOT_RUN in R01; no placeholder backend. |
| `proteomics verify --run runs/example --json` | Verify schema, content hashes, lineage, planned coverage and state consistency; identify integrity scope versus scientific validation. Does not fit or claim calibrated methods. | R01 verifies real foundation manifests; scientific consistency grows through typed artifacts without modifying CLI ownership. |

Configurations may be JSON or safe YAML parsed into the same JSON-compatible object. YAML aliases/tags that execute code are disallowed. No shell interpolation. Defaults do not depend on a username, drive or D1/T1 label. Standard output under --json is exactly one JSON object; human diagnostics are stderr. Every command has --help. Unknown command/option exits 2.

Later commands already represented in the twelve-packet scope are activated only with their real handlers: `resources prepare --manifest <file> --output <directory>` (R07), `report --run <directory>` (R10b), `compare --left <run> --right <run> --output <directory>` (R10b), `resume --run <directory>` (R11), `import-legacy --archive-root <local path> --output <directory>` (R12, Maintainer-only). QC-only uses runtime.scope=qc_only with `run`; template generation uses documented example copies, not an invented slash-command system. Unimplemented optional commands are not advertised as supported.

## Exit codes and states

0: requested execution completed, including zero discoveries or an explicitly optional scientifically inapplicable analysis. 2: input/config/design/eligibility rejection. 3: required runtime/resource/capability unavailable. 4: numerical/backend failure. 5: integrity/collision/cache error. 6: interruption/cancellation.

Stage states are NOT_RUN, RUNNING, COMPLETED, INAPPLICABLE, FAILED, CANCELLED and NOT_REQUESTED. A completed stage requires successful execution plus validated output schema/hashes. NOT_RUN includes an absent implementation or prerequisite; it is never PASS. INAPPLICABLE requires a scientific reason and an implemented eligibility check. NOT_REQUESTED means the frozen plan did not ask for the capability. Required scientific ineligibility rejects the plan (exit 2); optional ineligibility may be documented without running that method.

Run states are NOT_RUN, RUNNING, COMPLETED, PARTIAL, FAILED, CANCELLED. Model requiredness comes only from the frozen `model.execution_requirement`; the primary model is required. Scientific eligibility is evaluated without consulting installation state. A required stage failure/unavailability makes the run FAILED (exit 4/3); required dependants are NOT_RUN. A requested optional, scientifically eligible operational failure remains requested and produces PARTIAL with nonzero code 3/4 while preserving other real outputs. It cannot be removed from the plan or relabeled INAPPLICABLE. A failed report renderer cannot upgrade a failed model. A true zero-discovery result with all requested required work executed is COMPLETED/0. Validation-evidence PASS/FAIL/NOT_RUN is a separate vocabulary.

## Phase 1 layout

```text
<run>/
  config.resolved.json   plan.json   run_status.json   warnings.json   methods.md
  inputs/               # canonical values, ordered metadata and masks, manifest
  preprocessing/        # factors, transform lineage, feature/group eligibility
  designs/              # design/contrast matrices, term map, estimability, covariance
  qc/                   # sample_n.tsv, missingness.tsv, pca_scores.tsv, etc.
  dea/                  # zero_null.tsv OR treat.tsv, families.tsv, diagnostics
  report/index.html     # thin R10a report; links include complete dea tables
  provenance/           # input/resource/code hashes, actual environment/R session
  stages/<stage_id>/    # stage-result.json and immutable output manifests
  logs/                 # stdout, stderr and structured errors
```

`plan --output plan.json` stores supporting artifacts next to that plan under `plan-artifacts/`; plan.json contains only relative references to these immutable artifacts. `run` uses the same planner directly within its fresh run root. It does not rerun a logarithm on a transformed planning artifact. Plan hash is SHA-256 of canonical JSON (UTF-8, sorted object keys, compact separators, no NaN) excluding plan_hash and observational timestamps; referenced artifact and actual environment/code hashes are included. Repeated sample/row order is normalized by stable IDs before semantic comparisons, never by discarding duplicates.

`warnings.json` is an array of {code, stage_id, severity, message, artifact_ref}; severity is info, warning or error. methods.md is generated from resolved config and executed stage metadata. The R10a version is a small factual methods/status summary; full publication methods/limitations belong to FR-098/R10b. All files required by a completed stage have hashes in its manifest. A missing optional numeric table is represented by a typed stage reason, not fabricated numbers. A successful zero-row result retains the correct headers.

Phase 2 adds pathways/ and response/ with separate descriptive, equivalence, conjunction and independent-score tables. Phase 1 emits neither score tests nor an empty score-P-value table. Phase 3 adds full figures/report source data, lock/reproduction and calibration evidence. Every scientific result row follows [data-model.md](../data-model.md), including hypothesis_type, family_id, engine, estimable and n_obs_by_required_group.

<a id="runtime-extension-seam"></a>

## Frozen runtime extension seam — disjoint packet ownership

R01 owns the CLI, runtime and R dispatch entry points permanently. It defines a fixed internal capability-to-module/function map in `runtime.py`, not a user plugin registry: intake→intake.service, preprocessing→preprocessing_service, design→design_service, limma→inference_service, assay_engines→assay_service, resources→resources, pathways→pathway_service, response→response_service, report_stub→reporting.stub, report_full→reporting.full, reproduction→reproduction, legacy→legacy_service. All paths are under proteomics_pipeline. Each later module exports `capabilities()` and `execute(request)` using the stage schemas; no module import installs/fetches anything or runs analysis. Nonexistent modules are unavailable, never placeholder successes. `capabilities()` returns only actually implemented capability IDs, dependencies and required R packages. R01 inspects the finite map; arbitrary module paths from config are rejected.

The fixed R dispatch map uses the following internal handlers. Except for the exported foundation path function, every handler has signature `function(request) -> result`, where `request` is the already schema-validated decoded stage-request object and `result` is a decoded stage-result object. `dispatch_stage(request_path, result_path)` alone performs JSON file I/O, selects only a registered existing handler and writes the result atomically.

| Capability | Fixed R handler | Owning packet |
|---|---|---|
| `foundation.io_roundtrip` | exported `io_roundtrip(request_path, result_path)` | R01 |
| `preprocessing` | internal `preprocess_stage(request)` | R03 |
| `design` | internal `design_stage(request)` | R04 |
| `limma` | internal `limma_stage(request)` | R05 |
| `assay_engines` | internal `assay_engine_stage(request)` | R06 |
| `resources` | internal `mapping_stage(request)` | R07 |
| `pathways` | internal `pathway_stage(request)` | R08 |
| `response` | internal `response_stage(request)` | R09 |

Absent handlers are operationally unavailable. There is no separate registration metadata file: R01's hardcoded `dispatch.R` map stores each capability and the fixed function name above, and handler existence in the installed package namespace is the registration test. Later packets add only their owned handler file with that exact function; they do not edit `dispatch.R`, `DESCRIPTION` or `NAMESPACE`.

The R package exports only an implemented `dispatch_stage(request_path, result_path)` and `io_roundtrip(request_path, result_path)` at foundation. Later R functions are internal package functions, selected through a fixed dispatch map by actual existence; no nonexistent functions are exported in NAMESPACE. DESCRIPTION at R01 separates foundation Imports from future optional Suggests, never fabricates package versions. Shared R utilities expose typed file/result helpers; later files do not edit DESCRIPTION/NAMESPACE just to expose internal methods. Any unavoidable dependency/API correction is a Maintainer spec amendment before dispatch, not permission to reopen another packet's files.

R03 owns its detection-only Fisher/BH endpoint in `detection.R`; that narrow function is not the reusable cross-packet adjustment API. R05 owns the shared exact-fit and general family-adjustment APIs, and R06/R08/R09 consume them read-only. Phase 1 integration reruns V029 after R05 to confirm detection-family separation remains intact. R10a owns the stub/assembler; R10b adds a full renderer and reads the same ReportData, without changing the stub. R11 adds reproduction/cache services behind the established interface rather than reopening runtime/provenance. R12 documentation uses distinct release-guide paths. The [packet index](../packet-index.md) is the sole ownership list.

## Python ↔ R stage protocol

Request fields: schema_version=1.2.0; run_id; stage_id; capability; plan_hash (null only for foundation/validation work); inputs (array of artifact_id, path, sha256); output_temp_dir; config_path; parameters (a JSON object validated for that capability, never R source); rng (seed, kind, threads). Paths with spaces/Unicode are passed as argv list entries to `Rscript --vanilla scripts/maintained/run_stage.R --request <file> --result <file>` using shell=false. The wrapper calls the installed package, not `source()` on a user-supplied path.

Result fields: schema_version, run_id, stage_id, capability, plan_hash, state, reason_code, message, started_at, finished_at, exit_code, outputs (artifact_id, relative_path, sha256, result_type), warnings, session_info_path. R captures warnings and sessionInfo in the same process. The orchestrator validates matching identifiers/hash, child exit and every manifest before atomic promotion. A missing result after a child crash produces an orchestrator-generated FAILED record, clearly sourced to the orchestrator; no fake R session is created. Only output_temp_dir is writable by that stage.

An exclusive run lock prevents two writers. Stage output lives in a same-filesystem temporary directory; completed promotion uses an atomic rename after validation. Interrupted temporaries are never valid caches. Lock recovery needs an explicit recorded stale-owner check, not blind deletion. Cancellation terminates the owned child, saves logs and records CANCELLED. No background daemon/relay or global configuration mutation is required.

A required statistical stage with eligible numerical failures is FAILED and its family is incomplete. Preserve its partial raw/diagnostic outputs beneath the failed stage's diagnostic area, with hashes and explicit incomplete status; do not promote them as completed inference inputs. The report may show the failure and link diagnostic tables labeled incomplete, but downstream statistics never consume them. Prespecified feature exclusions/nonestimability are different: a successful fit of the remaining eligible universe can complete and still exports excluded rows with reasons. This is the SM12/SM25 status behavior tested by V048/V091/V099.
