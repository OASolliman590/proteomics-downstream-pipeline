# R01 implementation brief — the only first packet

**Status:** 1.2.0-frozen, ELIGIBLE BUT NOT AUTHORIZED FOR DISPATCH. Baseline inspected: `3aefdd95c46a1c56dae89e79ddf441b6cc978148`. No maintained package or engine was implemented in the specification pass. Begin only after the user explicitly says GO R01 using the recorded frozen contract identity.

## Deliverable and stopping boundary

Create a genuine installable Python CLI and minimal real R I/O package with strict schema/default handling, state/provenance, atomic publication, doctor and tests. This packet owns FR-001–FR-010 / T001–T010 / V001–V010 only. Full intake/QC/design/limma is R02–R05; do not implement any of it to make a demo succeed. A real I/O roundtrip is required; a fake limma result, dummy P table or mocked R run is forbidden.

The [packet index](packet-index.md) and [exact R01 ownership](packet-ownership.json) are mandatory. The [CLI/stage seam](contracts/cli-and-artifacts.md), [analysis schema](contracts/analysis.schema.json), [request schema](contracts/stage-request.schema.json), [result schema](contracts/stage-result.schema.json), [run-state schema](contracts/run-status.schema.json), [semantic rules](contracts/semantic-validation.md) and [data model](data-model.md) define fields; do not invent replacements.

## First files that must become real

| Files (future R01 targets) | Required first behavior |
|---|---|
| pyproject.toml; src/proteomics_pipeline/__init__.py, __main__.py, cli.py | setuptools-based install, package name proteomics-downstream-pipeline, development version 0.1.0.dev0; console entry `proteomics = proteomics_pipeline.cli:main`; python -m uses the same main; help/version need no R/data. This version denotes a foundation development build, not a released Phase 1 product. |
| config.py; schemas/analysis.schema.json; errors.py | Safe JSON/YAML loading (reject non-JSON NaN/Infinity constants and executable YAML tags), only the four documented defaults, strict Draft 2020-12 structural checks and data-independent cross-reference validation; JSON pointer error records. Copy the candidate contract schema exactly into package data, record its hash. |
| runtime.py, provenance.py, paths.py; three protocol schema copies | Typed states, fixed lazy capability map, actual subprocess management, exclusive run lock, verified atomic publication and hash-based manifest checks. No arbitrary user-specified import/source path. |
| doctor.py | Real Python/R/package discovery. Missing binary/package/resource is NOT_AVAILABLE with exit 3 when required. No hidden install/download. |
| r/proteomicsCore/DESCRIPTION, NAMESPACE, R/io.R, R/dispatch.R | Real installable package with only implemented dispatch_stage/io_roundtrip exports; foundation JSON/TSV I/O, warnings/session capture, safe literal argument handling. No statistics adapters. |
| scripts/maintained/run_stage.R, test_r.R, check_baseline.py, write_solved_versions.py, bootstrap.sh, bootstrap.ps1 | Actually runnable wrapper, test harness, byte-preservation check, solved-version recording and project-local explicit setup. |
| configs/examples/example-*.json and configs/examples/fixtures/ | Exact copies of the three public candidate examples and tiny fixture bytes, preserving relative paths and input hashes. These are examples, not fitted expected values. |
| tests/contract/foundation/, tests/unit/foundation/, tests/fixtures/foundation/, R tests/testthat/test-foundation.R | Real positive/negative tests for V001–V010. Preserve tests/test_runner.py unchanged. |
| docs/development.md and docs/validation/002-foundation/ | Actual exercised commands, results, limitations and receipt; baseline-hashes.json records a verified inventory, never rewrites recovered provenance. |

Use Python standard library plus jsonschema/PyYAML for the initial configuration/runtime layer; use pytest for Python tests and jsonlite/testthat for real R I/O/tests. Declare known later analysis dependencies as optional package extras/Suggests without importing them at startup; required packages are enforced per actual capability. Record actual solved versions, rather than pretending this document is a tested lockfile. R11 owns requirements.lock/renv.lock. Do not invent a License field/file to satisfy packaging tools. If missing owner/license metadata prevents actual installation, retain the real error and mark that gate NOT_RUN/unresolved; continue independent foundation work but do not close R01. No public package upload is authorized.

## Exact function and I/O seam

`cli.main(argv=None) -> int`; `config.load_config(path) -> dict` returns only resolved/schema-validated configuration, not a scientific plan. `runtime.capabilities() -> list[dict]` inspects the fixed module map. Each later owned service exports `capabilities()` and `execute(request: dict) -> dict` per stage schemas. `provenance.canonical_json_sha256(value: dict) -> str` implements the plan-hash canonicalization in the shared contract and can hash a complete synthetic plan envelope for pre-planner contract tests; it does not claim a real AnalysisPlan exists. R01 may inspect absent modules, but absence is E_CAPABILITY_NOT_IMPLEMENTED/NOT_RUN, never a placeholder success. The exact later R handler names/signatures are frozen in [cli-and-artifacts.md](contracts/cli-and-artifacts.md#frozen-runtime-extension-seam--disjoint-packet-ownership); R01 registers those names but does not implement nonexistent handlers.

`io_roundtrip` is a foundation-only test capability. Its parameters object has exactly `matrix_input_id` and `metadata_input_id`, naming two declared input artifacts. It reads their hashed paths, checks unique matrix/metadata keys, writes matrix.tsv and metadata.tsv beneath output_temp_dir with preserved values/NA/UTF-8 fields, captures a real sessionInfo.txt and returns a schema-valid stage result. It computes no scale transform, model, P value, imputation or biological claim. This narrow test I/O does not count as completion of R02's intake scope.

R result_path is supplied outside its promoted output manifest but inside the orchestrator's temporary stage directory. Output relative paths must resolve within output_temp_dir. Validate matching run_id/stage_id/capability/plan_hash, actual process exit and every output hash before promotion. A child crash cannot provide invented R session metadata. Root run_status.json uses the supplied schema; every referenced result path and artifact must remain within the run root after resolution. Empty or nonexecuted runs cannot verify as scientifically completed.

`validate --schema-only` genuinely works in R01. Default full validate, plan and run recognize their arguments but fail honestly until the real downstream handlers exist; do not generate a fictitious AnalysisPlan or fit. `verify` verifies actual foundation manifests/state and labels its scope; it does not declare calibrated science. This is an explicit unavailable capability, not a replacement backend.

## Commands to create and actually exercise

The following are future R01 gate commands; none is claimed to work at the recovered baseline:

```bash
python -m venv .venv
# Activate .venv using the documented platform-specific command.
python -m pip install -e ".[test]"
proteomics --help
python -m proteomics_pipeline --version
proteomics validate --config configs/examples/example-independent.json --schema-only --json
proteomics doctor --json
python -m unittest discover -s tests -v
python -m pytest tests/contract/foundation tests/unit/foundation -q
Rscript --vanilla scripts/maintained/test_r.R --suite foundation
python scripts/maintained/check_baseline.py
```

Bootstrap uses an existing documented R binary and creates only project-local Python/R libraries. It invokes a real R CMD INSTALL into that local library and the commands above. The test harness accepts --suite foundation (and unit as the same foundation selection), runs actual testthat files, and returns nonzero for no selected executable tests. Later scientific suites discover only the approved named tests; an absent method suite is NOT_RUN, not PASS. Installation may require explicitly permitted network access; analysis never does.

Expected unsupported checks during R01: full validate/plan/run on a statistical config return exit 3/E_CAPABILITY_NOT_IMPLEMENTED with no fit artifacts. With R deliberately unavailable, doctor/I/O gates report unavailable/NOT_RUN; --help/version and schema-only validation still work. The existing historical Python tests must remain intact and continue to be evaluated separately.

## Fixtures and oracles

Use only the eight-feature/twelve-observation independent fixture and eight-feature/eight-observation paired fixture under contracts/fixtures, plus a tiny path/Unicode/NA roundtrip. For argv safety use a directory containing a space, Unicode and shell metacharacter; independently check that no outside sentinel is created. For failures use actual subprocess exit/error, not a mocked R backend. For atomicity use two real writer processes and terminate one before promotion. For source preservation mutate only scratch copies and show the original registry detects it.

Each V001–V010 case in [the foundation spec](../002-foundation/spec.md) contains its fixture, oracle, exact assertion and negative case. Those cases are the completion criterion. A proposed test/evidence filename is not evidence; populate actual artifacts only after executing it. R01's packet receipt includes changed paths, exact commands/results, actual software versions/hashes and remaining external blockers. Maintainer independently verifies before advancing to R02. No commit/push or spec modification is authorized to the packet implementer.
