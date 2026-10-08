# Using the maintained pipeline (synthetic examples)

Everything below runs offline on the public synthetic examples. Real data
must stay outside the repository.

## Setup (once)

```bash
~/.local/bin/uv venv --python 3.13 .venv
~/.local/bin/uv pip install --python .venv/bin/python -e '.[test]'
export R_LIBS_USER="$PWD/.r-lib" LANG=en_US.UTF-8
Rscript --vanilla scripts/maintained/install_r_dependencies.R
R CMD INSTALL --library="$R_LIBS_USER" r/proteomicsCore
.venv/bin/proteomics doctor --json
```

Phase 2 engines additionally need Bioconductor `DEqMS`, `proDA` and `fgsea` in the
same library (`BiocManager::install(c("DEqMS","proDA","fgsea"), lib = Sys.getenv("R_LIBS_USER"))`).

For the pinned, reproducible environment (R11) restore the locks instead and qualify them:

```bash
.venv/bin/python -m pip install -r requirements.lock      # or: uv pip install -r requirements.lock
Rscript --vanilla -e 'lib <- Sys.getenv("R_LIBS_USER"); .libPaths(c(lib, .Library)); install.packages("renv", lib = lib); renv::restore(lockfile = "renv.lock", library = lib, prompt = FALSE)'
.venv/bin/python -c "from proteomics_pipeline import reproduction as r; print(r.check_r_lock('renv.lock')['qualified'], r.check_python_lock('requirements.lock'))"
```

## Workflow

```bash
.venv/bin/proteomics validate --config configs/examples/example-independent.json --json
.venv/bin/proteomics plan --config configs/examples/example-independent.json --output runs/plan-only/plan.json --json
.venv/bin/proteomics run --config configs/examples/example-independent.json --output runs/example --json
.venv/bin/proteomics verify --run runs/example --json
.venv/bin/proteomics report --run runs/example --output runs/example-report --json
.venv/bin/proteomics run --config tests/fixtures/reports/demo-phase2.json --output runs/demo-phase2 --json
.venv/bin/proteomics compare --left runs/example --right runs/demo-phase2 --output runs/comparison --json
```

On Windows (PowerShell) the same workflow is:

```powershell
$env:R_LIBS_USER = "$PWD\.r-lib"
.venv\Scripts\proteomics.exe validate --config configs/examples/example-independent.json --json
.venv\Scripts\proteomics.exe plan --config configs/examples/example-independent.json --output runs/plan-only/plan.json --json
.venv\Scripts\proteomics.exe run --config configs/examples/example-independent.json --output runs/example --json
.venv\Scripts\proteomics.exe verify --run runs/example --json
.venv\Scripts\proteomics.exe report --run runs/example --output runs/example-report --json
.venv\Scripts\proteomics.exe run --config tests/fixtures/reports/demo-phase2.json --output runs/demo-phase2 --json
.venv\Scripts\proteomics.exe compare --left runs/example --right runs/demo-phase2 --output runs/comparison --json
```

`python -m proteomics_pipeline <command> ...` works on every platform. Forward-slash paths such as `runs/example` work unchanged. All files are read and written as UTF-8 on every platform, and R child processes always run with UTF-8 input/output. The plan hash is platform-independent: the same configuration and inputs give the same hash on Linux, macOS and Windows (line endings are normalised; the solved environment is recorded in the plan but outside the hash, see D-42).

- `validate` runs intake, preprocessing and design checks without fitting.
- `plan` freezes the AnalysisPlan next to `plan-artifacts/`. The plan hash covers configuration, LF-normalised input and planning-artifact content, design, families and code identity; the solved environment and R session information are recorded beside it and protected by `integrity_sha256`.
- `run` re-plans in the run directory, fits, and writes `report/index.html` (thin report) and `report-full/index.html` (full report with figures).
- `verify` checks every manifest hash and the plan hash; it does not claim scientific validation.
- `report` re-renders the full report of an existing verified run into a new directory.
- `compare` writes a keyed comparison of two verified runs.
- `resume --run <dir>` re-uses every stage whose fingerprint (request, code manifest, environment) and outputs still verify and re-executes the first changed stage and everything after it.
- `resources prepare --manifest <file> --output <dir>` builds a local hashed snapshot of mapping/gene-set/orthology tables (explicit action; analysis never downloads).

Exit codes: 0 completed, 2 configuration/design/eligibility rejection, 3 missing capability or package, 4 numerical/backend failure, 5 integrity or collision error.

## Implemented capabilities

intake (R02), preprocessing/QC (R03), design and plan (R04), limma (R05), thin report (R10a), DEqMS/proDA (R06), resources and gene matrix (R07), pathways CAMERA/ROAST/fgsea/ORA (R08), response/equivalence/scores (R09), full report and comparison (R10b), PERMANOVA/PERMDISP (R13). Locks, offline reproduction, resume, the golden reference matrix, release calibration and the benchmark are R11 (see `docs/validation/012-validation`); the private legacy regression (R12) is Maintainer-only and not run here.

## Validation tooling (R11)

```bash
.venv/bin/python scripts/maintained/run_calibration.py            # release calibration, ~5 min; writes docs/validation/012-validation/calibration
.venv/bin/python scripts/maintained/benchmark.py --scale 0.1      # smoke benchmark; --scale 1.0 is the 20000 x 100 release workload
.venv/bin/python -c "from proteomics_pipeline import reproduction as r; r.write_acceptance_ledger()"   # docs/validation/acceptance.json
```
