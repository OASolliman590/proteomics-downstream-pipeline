# Target quickstart and implementation verification

This describes the maintained workflow to be built. The commands are not yet implemented at specification delivery. Maintainer must execute the finished versions and replace proposed examples with tested commands in the maintained user guide.

## New analysis

1. Install the project-local pinned Python and R environments using the implemented bootstrap for the platform; run `proteomics doctor --json`.
2. Generate a canonical LFQ/TMT/processed-protein configuration. Supply abundance, observation metadata, feature annotations and provenance. Declare scale, upstream normalization/imputation, biological-unit hierarchy, primary contrasts/hypothesis/families and optional eligible methods.
3. If pathways are requested, explicitly prepare pinned resource/mapping snapshots and verify their hashes. Production execution will remain offline.
4. Run `proteomics validate --config analysis.yml`, inspect the design/applicability report and fix input errors without altering source measurements.
5. Run `proteomics plan --config analysis.yml --output plan.json`; inspect the frozen estimands, biological n, filters, model method and testing families.
6. Run `proteomics run --config analysis.yml --output runs/run-01`; open the generated offline report.
7. Run `proteomics verify --run runs/run-01`. For a declared sensitivity configuration, use a different run folder and compare with `proteomics compare`.

The three JSON examples in contracts are schema examples with intentionally local fixture paths, not shipped study data. Foundation/intake slices create matching synthetic fixtures and YAML runnable examples. They demonstrate independent, paired and effect-threshold designs. The descriptive response epsilon and denominator floor are explicitly example policies, not universally optimal biological thresholds.

## Archived project

Only Maintainer accesses `<private-archive-root>/private_evidence/original`. It runs the maintained import-legacy/regression command locally after implementation. Packet Implementer, Independent Reviewer and any documentation worker receive only synthetic/redistributable fixtures and non-private contract facts. New outputs go to a separate run directory. No original file or report is modified.

The importer records that sample/treatment/tissue/preprocessing provenance is partly inherited or unknown. It does not transform the historical assertions into externally confirmed biological facts. A legacy comparison uses matched settings/universes when possible, then documents scientifically intentional differences in the maintained default.

## Developer acceptance

Follow the exact gates in plan.md and each completed slice. Validate configuration structure and semantics, exercise the real R bridge, compare numerical outputs against independent reference calls, run scientific calibration, render actual reports and verify final archive hashes. Missing R or resource access must be resolved or explicitly recorded as a remaining prerequisite, never hidden behind a successful Python-only test suite.
