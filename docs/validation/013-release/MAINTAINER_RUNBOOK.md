# Maintainer runbook — R12 private gates and release decisions

For Omar only. Every step runs on your own machine. Private inputs, private logs and filled-in configuration files stay **outside** the repository and are never uploaded, pasted into a chat or committed. Public summaries need your disclosure review first (spec 013, Boundary).

Shell setup for every step (macOS/Linux; for PowerShell see `docs/user-guide/usage.md`):

```bash
cd /path/to/proteomics-downstream-pipeline
export LANG=en_US.UTF-8 R_LIBS_USER="$PWD/.r-lib"
PRIVATE=/path/outside/the/repository/r12-private     # choose; must not be inside the checkout
mkdir -p "$PRIVATE"
```

## 1. V111 — verify the original archive (Maintainer-only)

```bash
.venv/bin/python scripts/maintained/legacy_regression.py verify-archive --archive /path/to/dapagliflozin_two_method_annotated_package_2026-07-21.zip > "$PRIVATE/v111-archive.json"; echo "exit=$?"
```

- **Pass:** exit 0 and `"verified": true`, `"file_count": 165`, `"archive_sha256_matches": true`. The archive and all 165 member hashes equal the public recovered inventory `docs/audit/archive_inventory.json`. Nothing is extracted or modified.
- **NOT_RUN:** exit 3 (`E_LEGACY_UNAVAILABLE`), meaning the archive is not on this machine. Record NOT_RUN and do not reconstruct inputs.
- **FAIL:** exit 5, meaning a hash differs. Do not continue; the count of problems is printed, and the member names are in the private JSON.

## 2. V114 — full maintained study run (Maintainer-only)

1. Write a study configuration outside the repository (start from `configs/examples/example-independent.json`; the post-DE block is in `specs/015-post-de-analysis/contracts/post-de.md`). Point `input` at your local study files. Keep unknown metadata unknown.
2. Run, verify and keep the run outside the repository:

```bash
.venv/bin/proteomics validate --config "$PRIVATE/study.json" --json > "$PRIVATE/v114-validate.json"
.venv/bin/proteomics run --config "$PRIVATE/study.json" --output "$PRIVATE/v114-run" --json > "$PRIVATE/v114-run.json"; echo "exit=$?"
.venv/bin/proteomics verify --run "$PRIVATE/v114-run" --json > "$PRIVATE/v114-verify.json"; echo "exit=$?"
.venv/bin/python scripts/maintained/check_baseline.py                    # the recovered originals are unchanged
```

- **Pass:** `run` exits 0 (or 3 with only optional stages NOT_RUN for absent optional packages, each named in the run status), `verify` exits 0, `check_baseline.py` prints `"valid": true`, and the run directory is new. Inspect `report-full/index.html` yourself. Unlike the synthetic tests, this is a real-data acceptance.
- **NOT_RUN:** study inputs or R unavailable. Synthetic data cannot substitute.

## 3. V112 — legacy core numerical comparison (Maintainer-local)

1. Copy `configs/local-legacy.template.json` to `$PRIVATE/local-legacy.json`. Set `archive`, `maintained_run` (`$PRIVATE/v114-run`), `output_dir` (a new directory under `$PRIVATE`), and one `tables` entry per historical table. Historical tables must first be keyed like the maintained table (`feature_id`, `contrast_id`; use the intake crosswalk for group names). Map the columns in `fields`.
2. Run:

```bash
.venv/bin/python scripts/maintained/legacy_regression.py run --config "$PRIVATE/local-legacy.json"; echo "exit=$?"
```

- The result `legacy_regression.json` (private) gives per field `MATCH` (frozen tolerances: effects abs ≤ 1e-10 or rel ≤ 1e-8; statistics/P abs ≤ 1e-8 or rel ≤ 1e-6), `DIFFERENT` (preserved, never reconciled) and `NOT_REPRODUCED` (rows on one side only). A maintained table byte-identical to a historical one is refused (`E_LEGACY_COPIED`).
- V112 passes only if a compatible historical R/resource environment was actually provisioned and the reference fits were rerun. Otherwise record NOT_RUN/NOT_REPRODUCED with the reason.

## 4. V113 — scientific-change reconciliation (Maintainer-approved)

Write the change ledger yourself, or approve one, from the public audit (`docs/audit/SCIENTIFIC_AUDIT.md`, `docs/audit/APPRAISAL.md`) and your approved local comparison summaries. It must explain the Negr1 auxiliary exception, the 2/70 versus 3/71 and 6/70 versus 7/71 arithmetic, and the intentional model, pathway and response differences, without production constants. It must not call the Negr1 treated-versus-control result a significance rescue, or principal null results discoveries. Place it in `docs/validation/legacy-comparison.md` only after disclosure review. The Claude route did not write it, because it would put study-derived numbers into the public repository.

## 5. Release decisions (operator decisions, AGENTS.md)

- **License and ownership:** choose them yourself. Record the decision in `docs/validation/013-release/operator-decisions.json` (`{"license": "...", "ownership": "...", "disclosure_review": "<date and scope>"}`) and add the license file in the same commit. Until then `build_release.py check` reports the license as unresolved, and a public export is refused.
- **Authorization:** for each release, write an authorization record outside the repository (see `docs/user-guide/release.md`).

## 6. V117–V119 on the final tree

```bash
.venv/bin/python scripts/maintained/build_release.py check; echo "exit=$?"            # V117: exit 0, no findings
.venv/bin/python -m pytest tests -q -p no:cacheprovider -rs                             # V119: full Python suite
Rscript --vanilla scripts/maintained/test_r.R --suite foundation                        # V119: every testthat file
.venv/bin/python scripts/check_spec_kit.py                                              # V119: spec kit
.venv/bin/python -c "from proteomics_pipeline import reproduction as r; print(r.check_r_lock('renv.lock')['qualified'])"
.venv/bin/python scripts/maintained/build_release.py export --version X.Y.Z --destination /path/outside/release-root --authorization /path/outside/authorization.json   # V118
```

V119 also needs the CI workflow (`.github/workflows/maintained.yml`) to pass on the pushed release commit, on Ubuntu and Windows. A green result for an earlier commit does not count.

## 7. Recording results

Update `specs/001-downstream-proteomics/traceability.json` for V111–V120 with the date, the commands, the exit codes and a non-private evidence summary. Regenerate the ledger with `.venv/bin/python -c "import proteomics_pipeline.reproduction as r; r.write_acceptance_ledger()"` (never hand-edit it). Run `tests/regression/test_traceability_reconciliation.py`. Private outputs stay in `$PRIVATE`.
