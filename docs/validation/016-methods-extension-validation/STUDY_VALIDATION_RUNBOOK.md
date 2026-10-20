# Study validation runbook (R15d; slots V189–V194, deferred)

Status: the real-study validation is **deferred by the Maintainer (2026-10-13)**. The framework (manifest, validate/run/verify,
public-data approval, private status records) is implemented and tested on synthetic studies (V185–V188). The five design
slots and the reference re-run (V189–V194) are NOT_RUN in this repository. Run them locally, on your own machine, with the
steps below. Nothing in this file contains a study name, a file name, a sample name or a path.

## Before any study

1. Work from a clean checkout of the branch, with the project R library: `export R_LIBS_USER="$PWD/.r-lib"` and
   `export LANG=en_US.UTF-8`. The Python environment is `.venv/bin/python`.
2. Keep every private input and every private output **outside the repository**. Choose a run root such as
   `$HOME/private-runs/<study>` (not under the repository). Only the non-identifying status record may be committed.
3. Write one study manifest per study next to its inputs, outside the repository. The schema is version `1.0.0`:
   `study_id` (a local identifier, no slashes or spaces), `origin` (`private` or `public`), `slot`, `design`, `input`
   (`format`, and `files` with `name`, `sha256` and `bytes`), `organism`, `columns` (`subject`, `group`, `batch`, `phenotype`),
   `licence`, `config` (the analysis configuration, relative to the manifest) and, for a public study, `approval_ref`.
   Compute each file's SHA-256 with `shasum -a 256 <file>`; a changed file is refused (E_MANIFEST_HASH).
4. Choose the slot and the design that match the data: slot `dda_label_free` with design `two_group`; slot
   `paired_repeated` with `paired` or `repeated`; slot `tmt_batched` with `two_group` or `multi_group` and a batch column;
   slot `dia` with input format `diann` or `spectronaut`; slot `multigroup_phenotype_nonhuman` with `multi_group` and a
   phenotype column; slot `reference_welch` with `multi_group` (the three-group serum re-run, Welch engine).

## Slot by slot

| Slot | Case | Manifest | Command outline | Expected evidence |
|---|---|---|---|---|
| 1, label-free DDA (MaxQuant-style protein groups) | V189 | `slot: dda_label_free`, `input.format: maxquant_proteingroups` | `multistudy.run_study(manifest, run_root)` | mapping recorded; reverse hits and contaminants removed and counted; planted contrasts recovered |
| 2, paired or repeated measures | V190 | `slot: paired_repeated`, `design: paired` (or `repeated`) | the same call | subject structure preserved; the paired variant or subject-blocked analysis ran; repeated observations never treated as independent |
| 3, TMT with plex or batch | V191 | `slot: tmt_batched`, `columns.batch` set | the same call | group-by-batch table in the plan; a confounded batch recorded as an adaptation and its term not fitted |
| 4, DIA (DIA-NN or Spectronaut style) | V192 | `slot: dia`, `input.format: diann` or `spectronaut` | the same call | the column mapping and its SHA-256 recorded; an unmapped column refused on its own |
| 5, multi-group with a continuous phenotype, non-human organism | V193 | `slot: multigroup_phenotype_nonhuman`, `columns.phenotype` set, `organism` declared | the same call | complete-case association; identifiers mapped only through the declared organism resource |
| Reference re-run, three-group serum (Welch t) | V194 | `slot: reference_welch`, `design: multi_group`, a config with `methods_extension.welch` enabled and the paper-compatible rule | the same call | PASS or FAIL with check counts; the `not_multiplicity_controlled` label beside BH q |

For each slot, run from the repository root:

```
.venv/bin/python -c "import json; from proteomics_pipeline.multistudy.runner import run_study, commit_status; r = run_study('<manifest path>', '<run root outside the repository>'); print(json.dumps(r['checks'], indent=1))"
```

Read the four checks (`validate`, `run`, `verify`, `acceptance`). Every one must be PASS. A FAIL is a finding, not a reason to change the
gate; record it and report it.

## Committing the result (only the status record)

```
.venv/bin/python -c "from pathlib import Path; from proteomics_pipeline.multistudy.runner import run_study, commit_status; r = run_study('<manifest>', '<run root>'); commit_status(r, Path('docs/validation/016-methods-extension-validation/R15d/status/<slot_or_label>.json'))"
.venv/bin/python scripts/maintained/verify_private_status.py
```

`commit_status` refuses (E_PRIVATE_IDENTIFIER) any record that would carry the study identifier or a file name, and
`verify_private_status.py` fails if one is present. Do not edit a status record by hand.

## Public datasets (only with an approval)

A public dataset is downloaded only after an approval record names its accession, its URL and the SHA-256 of the file.
Write the approval record outside the repository, then:

```
.venv/bin/python scripts/maintained/fetch_public_study.py --approval <approval.json> --accession <ACC> --dest <directory outside the repository>
```

A missing approval, a changed URL or a file whose SHA-256 differs is refused before anything is written. Analysis never downloads.

## Reference re-run (V194) notes

The re-run is Maintainer-only and private. Use the three-group serum configuration with the Welch engine declared under
`methods_extension.welch` (Phase 5), the paper-compatible raw-P rule enabled (its rows are labelled `not_multiplicity_controlled`),
and BH reported beside it. Commit only the PASS or FAIL status record.
