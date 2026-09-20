# CLI, execution and artifact contract v1.0.0

These interfaces are implementation targets. The baseline currently exposes only the recovery scripts. A command becomes supported only when its behavior and acceptance tests exist.

## Commands

| Command | Behavior | Failure contract |
|---|---|---|
| `proteomics doctor --json` | Read-only inspection of binaries, packages, platform, schema and resource availability | Nonzero when a requested required capability is missing; versions are actual inspected values |
| `proteomics init-config --profile canonical-lfq --output analysis.yml` | Write a documented valid template using synthetic/local-relative paths; no data discovery | Refuse overwrite unless an explicit scoped overwrite option is supported and tested |
| `proteomics validate --config analysis.yml [--json]` | Parse source headers/schema/metadata, scale, grain, design and method eligibility; no model fits or source writes | Structured errors with path, field and reason; JSON safe for automation |
| `proteomics plan --config analysis.yml --output plan.json` | Validate fully and freeze ordered IDs, filters, models, contrast/family definitions and provenance | Fail aliased/nonestimable scientific plan or missing required snapshot; no silently chosen substitute |
| `proteomics run --config analysis.yml --output runs/run-01` | Validate, freeze plan, execute all eligible required stages and write offline report | Refuse existing unrelated run; stop failed required stage, preserve logs and honest partial report |
| `proteomics qc --config analysis.yml --output runs/qc-01` | Intake/preprocess/QC only, labeled limited scope | Never claims differential/pathway execution |
| `proteomics resume --run runs/run-01` | Continue original frozen plan using only verified compatible artifacts | Changed config/source/code/resources require a new run or explicit rebuilt plan; cannot mix stale and new science |
| `proteomics report --run runs/run-01` | Render solely from validated scientific artifacts and statuses | Missing inputs stay missing; no recomputation of P/q/model claims inside templates |
| `proteomics verify --run runs/run-01 [--json]` | Validate hashes, schema, lineage, required stages, semantic consistency and no false completion | Failure and NOT_RUN distinguish integrity versus scientific-validation scope |
| `proteomics compare --left runs/a --right runs/b --output comparison` | Compare effects/universes/config/resources and classify intentional/unknown differences | Join by stable scientific keys, not row order or contrast nickname |
| `proteomics resources prepare --manifest resources.yml --output cache` | Explicit external-resource acquisition/preparation with version/hash/terms | No opaque latest download, no credentials printed, no production auto-fetch |
| `proteomics import-legacy --archive-root path --output canonical` | Local conversion of audited archive into canonical bundle/config with inherited unknown metadata | No archive instructions executed; no original files modified or uploaded |

All commands support `--help` and emit a final concise status. Programmatic results use JSON; human diagnostics go to stderr or a structured log. CLI and Python API call the same core services. Defaults never depend on a particular analyst's home directory, drive letter or sample-prefix spelling.

## Exit codes

0 = requested scope completed (possibly scientific warnings, never hidden required failure). 2 = input/config/schema/eligibility error. 3 = required environment/resource unavailable. 4 = statistical/numerical stage failure. 5 = artifact integrity/cache conflict. 6 = interrupted/cancelled. A run may successfully describe an inapplicable optional analysis, but a requested required backend that is not installed exits 3 rather than pretending inapplicability. A result with no significant proteins exits 0.

## Run layout

```text
run-01/
  run.json config.resolved.yml plan.json
  inputs/ canonical.tsv observations.tsv features.tsv masks/ input-manifest.json
  preprocessing/ transformations.tsv normalization-factors.tsv eligibility.tsv
  qc/ tables/ figures/ figure-data/
  designs/ <design-id>.tsv contrasts.tsv feature-estimability.tsv
  models/ <model-id>/ fit.rds fit-metadata.json diagnostics/
  results/ differential.tsv families.tsv detection.tsv
  pathways/ results.tsv mapping.tsv universes/ memberships/ diagnostics/
  response/ feature-response.tsv equivalence.tsv scores.tsv score-tests.json
  report/ index.html assets/ methods.md limitations.md report-data.json
  provenance/ manifest.json environment.json package-versions.tsv resources.json
  logs/ <stage>.jsonl <stage>.stdout.log <stage>.stderr.log
  stages/ <stage-id>/ status.json
```

Optional outputs have a status/eligibility artifact even when no numeric table is applicable. Do not create a table with fabricated columns/values merely so a report can read it; an empty valid schema is appropriate for a successful zero-row result. Each model/contrast/resource subdirectory uses validated safe IDs to prevent path traversal.

## Stage protocol

The orchestrator supplies a JSON stage request containing schema_version, run_id, stage_id, plan_hash, canonical input paths/hashes, explicit output temporary directory, model/contrast definitions and deterministic RNG settings. R validates the request, performs actual work, captures warnings/sessionInfo in the same process and writes stage-result.json plus typed outputs. Output status is derived from the call and validations, not from a prose report. The orchestrator checks process exit, output completeness and hashes before promoting the stage atomically.

Timeout and cancellation preserve partial logs and mark the stage incomplete. Killing a relay or R child is not completion. Do not start hidden background daemons; bounded processes are owned by the task. Resource limits are explicit and benchmarked. No interpolation of user-provided formulas or filenames into shell command strings.

## Reproducibility

Exact content hashes identify input/resource/code artifacts. Semantic output comparisons align sorted keys and compare numeric fields within prespecified backend tolerance; unchanged meanings cannot be masked by sorting away duplicate keys or ignoring changed n/universes. Save seed, RNG kind, threads and package versions for stochastic methods. Production network access is unnecessary when required resources are available; missing cache fails cleanly. Installation belongs to bootstrap/doctor, not run.

## Release protocol

The maintainer exports reviewed source/spec/tests/docs and provenance to a fresh `<release-root>/releases/<version-or-commit>/` folder. Original `<private-evidence-root>/`, original ZIP and original `package_manifest.json` remain unchanged. Local scientific run results may live in a separate explicitly local results folder; they are not implicitly part of a GitHub source release. Build output checks for included credentials/private data and records a file manifest. No license is selected without owner input.
