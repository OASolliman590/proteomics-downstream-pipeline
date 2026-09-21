# Proteomics downstream pipeline

This repository contains a **recovered study baseline** and a **frozen v1.2.0 Spec Kit for a maintained successor**. It is **not yet a validated general platform**, and **no implementation packet starts without explicit authorization**.

## What exists today

The recovered analysis is inspectable in pipeline/ and legacy/, with findings in [docs/PIPELINE.md](docs/PIPELINE.md) and the immutable docs/audit/ records. The runnable surface is scripts/run_pipeline.py, six historical R stages, pipeline/20_external_limma_intake.py and tests/test_runner.py. The recovery did not rerun R. environment.yml is not a lockfile, gene-set snapshots are not pinned, and no repository license has been selected. Private study files are not included and must not be uploaded or committed.

The existing inspection/test commands are:

```bash
python scripts/run_pipeline.py --help
python -m unittest discover -s tests -v
python scripts/check_spec_kit.py
```

They do not constitute scientific validation of a maintained successor. The historical runner/importer usage remains documented in [docs/PIPELINE.md](docs/PIPELINE.md); its code and audit are unchanged.

## Maintained successor — specified, not shipped

Start with [START_HERE](specs/001-downstream-proteomics/START_HERE.md), [PHASES](specs/001-downstream-proteomics/PHASES.md) and the [packet index](specs/001-downstream-proteomics/packet-index.md). The kit is **v1.2.0-frozen**; the [scientific methods contract](specs/001-downstream-proteomics/contracts/scientific-methods.md) is v1.1.0. The only first implementation packet eligible for explicit authorization is [R01](specs/001-downstream-proteomics/IMPLEMENTATION_BRIEF_R01.md).

| Scope | Status / intended milestone |
|---|---|
| Recovered study scripts + audit | Current baseline-recovered; not a general-platform validation claim. |
| R01–R05 + thin R10a report | Future v0.1-limma-core. |
| R06–R09 qualified methods/resources/response | Future v0.2-qualified-methods. |
| Full R10b + R11 + R12 | Future v1.0-defensible, only after its actual gates. |

The target is declared protein-abundance intake/QC → frozen design → qualified limma or eligible adapters → design-valid enrichment → descriptive response and restricted independent inference → honest offline reporting. It excludes raw-MS search/quantification, PTM localization, single-cell, classifier training, web UI/databases, network deployment and causal drug-mechanism claims.

Paths such as src/proteomics_pipeline/, r/proteomicsCore/, configs/examples/ and renv.lock are implementation targets, not existing maintained software at this baseline. The three [configuration examples](specs/001-downstream-proteomics/contracts/example-independent.json) use [tiny synthetic fixtures](specs/001-downstream-proteomics/contracts/fixtures/README.md), not study evidence. No phase is completed merely because its specifications exist.

See [progress](specs/001-downstream-proteomics/progress.md), [hardening report](docs/SPEC_HARDENING_REPORT.md) and [release checklist](docs/RELEASE_CHECKLIST.md) for review, verification and unresolved release boundaries.
