# R15d receipt — Multi-study validation framework (FR-185 to FR-194; SM46; ADR 0011)

Packet R15d (FR-185 to FR-194 / T185 to T194 / V185 to V194) was implemented on 2026-10-13 by the **Claude route (Haiku builder, Opus review)**. The **framework** (FR-185 to FR-188) is implemented and tested on synthetic studies: V185 to V188 are PASS. The real-study slots (V189 to V194) are **NOT_RUN**: the Maintainer deferred real-study validation on 2026-10-13. They are run locally with the study runbook. Independent audit: **NOT_RUN**.

## Identities

- Contracts: scientific methods v1.4.0 (SM46); `specs/016-methods-extension-validation/contracts/methods-extension.md` (study manifest, error codes, outputs).
- Working-source manifest (9 files: the framework, the two maintained scripts, the integration test, the fixture factory and the runbook): `evidence/working-source-manifest.txt`.
- `verified_commit`: **null** (local commits on `claude/methods-016`, not pushed).
- Environment: as R15c (macOS, Python 3.13.15, R 4.6.1, `.r-lib`, `LANG=en_US.UTF-8`).
- Commit: `4b3fb9f` (framework, tests, scripts, runbook, amendment A-2026-10-01-48).

## What was built

- `src/proteomics_pipeline/multistudy/manifest.py`: the study manifest (schema 1.0.0). Every input file is checked against its declared SHA-256 and byte count before anything is validated or run (`E_MANIFEST_HASH`). Unknown slots, designs and formats are refused (`E_MANIFEST_DESIGN`). A study identifier with a path separator or a space is refused (`E_PRIVATE_IDENTIFIER`). A public manifest needs a matching approval record (`E_PUBLIC_DATA_UNAPPROVED`). The slot adapters: the column mapping of each format (`map_columns`), recorded with its SHA-256; an unmapped column is refused on its own (`E_MAPPING_UNMAPPED`); the MaxQuant reverse-hit and contaminant filters with counts (`maxquant_filters`); the TMT group-by-batch confounding check, where a confounded batch becomes an adaptation and its term is not fitted (`tmt_batch_confounding`).
- `src/proteomics_pipeline/multistudy/runner.py`: `run_study` runs validate, run and verify for one study in a run root the caller chooses (never in the repository), each reporting its own state, then the design acceptance. `status_record` keeps only the slot, design type, state and check counts. `commit_status` refuses (`E_PRIVATE_IDENTIFIER`) any record that would carry the study identifier or an input file name.
- `src/proteomics_pipeline/multistudy/acceptance.py`: the design acceptance (tested endpoints in `dea/zero_null.tsv`; the subject structure for paired and repeated designs).
- `scripts/maintained/fetch_public_study.py`: downloads a public file only against a matching approval (accession, URL and SHA-256). A missing approval refuses before any network call; a changed file is refused and removed; the destination must lie outside the repository. The network layer is an injectable opener; the tests patch it, and no real dataset is downloaded.
- `scripts/maintained/verify_private_status.py`: checks committed status records for extra keys and for file names or paths.
- `docs/validation/016-methods-extension-validation/STUDY_VALIDATION_RUNBOOK.md`: how to run each slot locally, what to commit (the status record only), and how to fetch a public dataset with an approval.

## Gates executed

| Command | Exit | Result | Log |
|---|---|---|---|
| `python -m pytest tests/integration/test_multistudy.py -q -rs -p no:cacheprovider` | 0 | 15 passed, 0 skipped | `evidence/python-multistudy-tests.log` |
| `python -m pytest tests/regression tests/unit tests/scientific/reference/test_text_encoding.py -q -p no:cacheprovider` | 0 | 91 passed | `evidence/gates-regression-unit-encoding.log` |
| `python scripts/check_spec_kit.py` | 0 | PASS | `evidence/spec-kit-check.log` |

Full-suite results are recorded in the section "Full suites" below.

## Acceptance cases

| Case | Status | Positive | Negative |
|---|---|---|---|
| V185 | PASS | A valid manifest validates and hashes all four input files. | A changed input file fails `E_MANIFEST_HASH` before any run; an unknown design and an unknown format are each refused (`E_MANIFEST_DESIGN`); an identifier naming a path or a person is refused (`E_PRIVATE_IDENTIFIER`). |
| V186 | PASS | A synthetic three-group study runs through validate (PASS), run (COMPLETED), verify (PASS) and the acceptance (PASS); each step reports its own state (`run_study`). | An output changed after the run fails `proteomics verify` with `valid: false` and exit 5 (the existing code for a changed stage output is `E_INTEGRITY`; see the deviation below). |
| V187 | PASS | A public file is downloaded against a matching approval (the fake opener is called with the approved URL; the observed SHA-256 equals the approved one; the record holds the approval, the URL and the hash); the file is written outside the repository. | No approval: the network layer is never reached and no file is written. A changed file is refused (`E_MANIFEST_HASH`), a destination inside the repository is refused, and a public manifest without an approval record fails validation (`E_PUBLIC_DATA_UNAPPROVED`). |
| V188 | PASS | A private study commits only the status record (keys: schema, slot, design, state, checks, check counts); no study identifier, file name or path is in it; `verify_private_status.py` passes it. | A status record that would carry the study identifier is refused (`E_PRIVATE_IDENTIFIER`), and a committed record with a file name fails the verifier (exit 5). |

Slot adapters (FR-189, FR-191, FR-192; tested on synthetic inputs; they are not cases of V189 to V194): the DIA-NN columns map through a recorded mapping with its hash, and a missing `genes` column is refused alone (`E_MAPPING_UNMAPPED`); the MaxQuant filters remove one reverse hit and one contaminant and count them; a TMT batch without one group level is confounded, its adaptation is recorded and its batch term is not fitted, while a balanced layout fits the term.

## Real-study slots (NOT_RUN)

V189 (slot 1, label-free DDA), V190 (slot 2, paired or repeated), V191 (slot 3, TMT with batches), V192 (slot 4, DIA), V193 (slot 5, multi-group with a continuous phenotype, non-human organism) and V194 (the three-group serum re-run with Welch t) are **NOT_RUN** with the reason "real-study validation deferred by the Maintainer (2026-10-13); run locally with the study manifest". No real study was used, read or downloaded. The runbook (`STUDY_VALIDATION_RUNBOOK.md`) states how to run each slot, what to commit (the status record only), and how a public dataset is fetched after an approval.

## Deviations and doubts (for the reviewer)

1. **Changed stage outputs report `E_INTEGRITY`.** The existing `proteomics verify` reports a changed stage output as `E_INTEGRITY` ("stage output is outside its stage or hash-mismatched"), not `E_ARTIFACT_HASH` as the spec's error table names. The framework does not rename the existing code; the V186 negative accepts the existing code. The Maintainer should decide whether the spec name is to be added.
2. **Public-data destination.** A destination inside the repository is refused with `E_PUBLIC_DATA_UNAPPROVED` (the spec names no separate code for it); no new code was introduced.
3. **No status record is committed yet.** `R15d/status/` does not exist: no real study has run, and the synthetic tests write their records to temporary folders only.
4. **Approval records.** The approval format (accession, URL, file name, SHA-256, approver, date) is the implementer's; the Maintainer confirms it before any public download.
5. **Acceptance is simple.** The design checks are the tested-endpoint count and the subject structure. A real slot needs its own acceptance rules, to be specified with the Maintainer when the real study runs.
6. **Amendment A-2026-10-01-48** records the runbook (outside the R15d folder) and the traceability, ledger and task updates. Pending Maintainer confirmation.

## Full suites

See the table at the bottom of this file.
