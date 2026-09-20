# Proteomics pipeline: recovered baseline and downstream implementation Spec Kit

The detailed maturation mandate is in [Spec Kit START_HERE](specs/001-downstream-proteomics/START_HERE.md): twelve implementation packets, 120 requirement/task/acceptance links, scientific/data/CLI contracts, and a corrected [packet index](specs/001-downstream-proteomics/packet-index.md). The kit is a draft awaiting independent specification audit and maintainer freeze; no implementation packet is authorized yet. These documents specify the successor; the recovered functionality below is the current baseline.

Recovered and audited on 2026-09-12 from `dapagliflozin_two_method_annotated_package_2026-07-21.zip`.

**Status: recovered research code plus a corrected, not-yet-frozen implementation Spec Kit. Full maintained implementation, R execution and scientific validation remain outstanding.** The recovered code is a study-specific downstream analysis of processed log2 protein abundance. It does not perform raw mass-spectrometry identification or quantification.

Start with the [appraisal](docs/audit/APPRAISAL.md), [pipeline reconstruction](docs/PIPELINE.md), [scientific audit](docs/audit/SCIENTIFIC_AUDIT.md), and [technical audit](docs/audit/TECHNICAL_AUDIT.md). The [release checklist](docs/RELEASE_CHECKLIST.md) separates completed recovery from remaining work.

## Recovered analysis

Method A reconstructs a 3,714-protein, 20-sample matrix from four workbook sheets and standardizes the supplier's four limma tables. Method B refits the matrix using limma, first in a joint five-group model and then in separate PreDM and DM models. It adds QC, sensitivity fits, reversal summaries, rank-based pathway enrichment and figures. The archive describes four independent rats per group; the underlying provider confirmation is not included.

The four main disease/direct-treatment contrasts have no proteins at FDR < 0.05. All six axis-specific contrasts also have none at FDR < 0.10. An important exception to the archive's broad summary is **Negr1 (A0A8I6AJV3), joint-model PreDM+dapa versus CTRL**, with log2FC -1.07835 and FDR 0.0169815. This is a sensitivity contrast comparing treated disease with healthy controls; it does not demonstrate treatment rescue. Selected reversal-score tests and pathway findings retain the limitations in the audit.

## Folder boundary

Only this `repository/` folder is intended as the future Git working tree. Its sibling `private_evidence/` holds the original ZIP, byte-preserved extracted project and prepared local input. Data and generated outputs are excluded from the repository. Audit summaries include study-level results; review their disclosure together with authorship and licensing before public release.

```text
repository/
  pipeline/       six unchanged historical R stages and portable Method A intake
  legacy/         original finalizer and CSV export builder, retained for inspection
  scripts/        preparation, validation, execution and independent table audit
  tests/          synthetic input and failure-path checks
  docs/           reconstructed workflow, appraisal and audit evidence
private_evidence/ original ZIP, original/ export, prepared_run/
```

The R source hashes match the archive. Only Method A's analyst-specific default input path was removed; `--input` is now required. See [source_provenance.json](docs/audit/source_provenance.json). Scientific logic has not been silently corrected. The portable runner excludes the legacy finalizer because it hardcodes study claims and PASS statuses.

## Local use

Python 3.12.14 was used for input preparation, Method A intake and the unit tests. R 4.4.2 is the historical interpreter recorded in the archive. `environment.yml` lists the recovered R dependencies; it is an **untested recovery specification, not a lockfile**. Bioconda-based setup is intended for Linux/macOS. On native Windows, use an installed Rscript with the listed packages. The independent table-audit requirements pin the NumPy/pandas versions actually used.

From this folder, the packaged study inputs can be checked with:

```bash
python scripts/run_pipeline.py preflight --run-dir ../private_evidence/prepared_run
python -m unittest discover -s tests -v
```

Preflight exits nonzero when Rscript or required packages are unavailable. This happened during recovery: input validation passed, but no R runtime was found. A successful preflight alone does not establish numerical reproduction.

After provisioning R and reviewing the documented scientific limitations:

```bash
python scripts/run_pipeline.py run --run-dir ../private_evidence/prepared_run
```

Use `--rscript "path/to/Rscript"` when it is not on PATH. The runner executes six R stages in the archived order, writes each stage's log and same-process `sessionInfo()`, and stops on failure. Enrichment may retrieve msigdbr resources from the network; historical gene-set snapshots were not included. Results may therefore differ. Every execution uses a separate run directory, and existing execution results are never overwritten by the runner.

To prepare another fresh run from the archived Method B export:

```bash
python scripts/run_pipeline.py prepare --export ../private_evidence/original/02_METHOD_B_ABUNDANCE_REFIT_FULL_PIPELINE/01_FULL_VERIFIED_EXPORT --run-dir runs/reproduction-01
```

Preparation converts quoted CSV into the TSV files expected by R, preserves missingness and Method B metadata, validates sample/rat IDs, group counts, treatment labels and annotation alignment, and records input hashes. It does not apply log transformation, imputation or normalization. The validator intentionally enforces this historical study's five groups and four rats per group; it is not a generic experiment-design interface.

For Method A, supply the workbook and all output destinations explicitly:

```bash
python pipeline/20_external_limma_intake.py --input ../private_evidence/original/01_METHOD_A_SOURCE_PRECOMPUTED_LIMMA_INTAKE/01_RAW_SOURCE_WORKBOOK/limma_results.xlsx --outdir runs/method-a/data --results runs/method-a/results --report runs/method-a/report
```

Method A metadata uses PDM/DM_Treated/PDM_Treated and unresolved biological fields. It is not directly interchangeable with Method B metadata.

For independent arithmetic checks of the archived results, install `requirements-audit.txt` and run:

```bash
python scripts/check_data.py --original ../private_evidence/original --output runs-data-checks.json
```

This checks archived numbers; it does not rerun limma or fgsea. The JSON contains integrity exceptions as findings, so inspect its content rather than treating process success as a scientific PASS.

## Publishing

This repository is a sanitized source/specification export. It contains no private evidence directory or original study files. Review audit evidence, study-level disclosures, authorship, ownership and redistribution terms before changing repository visibility or publishing a release. No license has been invented. See [release checklist](docs/RELEASE_CHECKLIST.md). The current CI workflow checks recovered Python code and synthetic inputs only; it does not validate the future maintained R package, gene-set resources or biological findings.
