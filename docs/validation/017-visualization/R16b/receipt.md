# R16b Claude-route receipt: Prism `.pzfx` writer (ADR 0012)

Packet R16b (FR-199–FR-201 / T199–T201 / V199–V201) was implemented on 2026-10-08 by the Claude route (Haiku builder, Opus review), at the operator's direction. This is not work by the routes named in AGENTS.md. The receipt means ready for review, not accepted. Independent audit: **NOT_RUN** (the Opus review follows this receipt).

## Identities

- Contract: `specs/017-visualization/contracts/figures.md` ("Prism `.pzfx` table rules", "Dependencies"). Decisions: D-12 (round-trip precision), D-68 (three outputs), D-74.
- Working-source manifest: `evidence/working-source-manifest.txt`, 14 files, SHA-256 `bf41d8163410369931fa79f6ead19996c3545f9c60fce92636c64b2177ec62ef`. It covers the writer, its test and fixtures, the R read-back test, and the three shared-interface files of A-2026-10-01-30 (`r/proteomicsCore/DESCRIPTION`, `renv.lock`, `docs/methods/reason-codes.md`). The ledger, traceability, tasks and this receipt are excluded because they record the result.
- Shared-interface amendment made first: **A-2026-10-01-30** (`packet-ownership.json`): pzfx in the R Suggests and `renv.lock`, the reason-code catalogue, traceability, task ticks T199–T201, and the ledger.
- Optional read-back oracle: CRAN **pzfx 0.3.1** (with **xml2 1.6.0**) installed in the project-local `.r-lib` (`evidence/r-install.log` for the package install). `reproduction.check_r_lock('renv.lock')` reports qualified = True with no mismatches.
- `verified_commit`: **null**. The R16b commit is local on `claude/visualization` and is not pushed.
- Environment: macOS, Python 3.13.15 (`.venv`), R 4.6.1, `LANG=en_US.UTF-8`.

## Implementation

- `src/proteomics_pipeline/figures/pzfx_writer.py` writes the XML with `xml.etree.ElementTree` (standard library). The pzfx R package is not used to write.
- Table types: `Column` = Prism `OneWay`, `XY` = `XY`, `Grouped` = `TwoWay`. Values are written with `repr()`, which reads back as the same float. A missing value is an empty cell.
- The writer's self-check parses the bytes it wrote and compares every value and title with the source before the file is returned (`E_PZFX_VALUE_MISMATCH`).
- `tests/fixtures/figures/pzfx/make_fixtures.py` generates the synthetic source tables and the committed `.pzfx` fixtures. The source values come from fixed trigonometric formulas with planted offsets. No study data are used.

## Gates executed

| Command | Exit | Result | Log |
|---|---|---|---|
| `.venv/bin/python -m pytest tests/scientific/test_pzfx_writer.py -q -rs -p no:cacheprovider` | 0 | 11 passed | `evidence/python-tests-pzfx.log` |
| `Rscript --vanilla -e "testthat::test_file('r/proteomicsCore/tests/testthat/test-pzfx-readback.R', …)"` | 0 | 30 expectations, 0 failed, 0 skipped | `evidence/r-test-pzfx-readback.log` |
| `R CMD INSTALL --library="$R_LIBS_USER" r/proteomicsCore` | 0 | installed (DESCRIPTION changed) | `evidence/r-install.log` |
| `.venv/bin/python scripts/check_spec_kit.py` | 0 | PASS | `evidence/spec-kit-check.log` |
| `.venv/bin/python -m pytest tests/regression tests/unit -q -p no:cacheprovider` | 0 | 85 passed | `evidence/regression-unit.log` |
| Full Python suite, `.venv/bin/python -m pytest tests -q -p no:cacheprovider` (background, on the final tree after the review fix) | 0 | 475 passed in 25 min 36 s | `evidence/full-python-suite.log` |

## Acceptance cases

| Case | Status | Positive oracle | Negative case |
|---|---|---|---|
| V199 | PASS | Column tables read back from the XML with `xml.etree` equal the source TSV values exactly (three proteins, two groups, six replicates). Round-trip precision (17-digit and tiny values). R: `pzfx::read_pzfx` returns the same values (as numeric). Prism opening: manual, NOT_RUN. | one changed value (by 1e-9) fails the comparison; values rounded to three decimals fail. |
| V200 | PASS | volcano (200 rows) and ROC (30 rows) XY tables: every X and Y equals its source, one point per row, table type XY. R read-back of both tables. Prism opening: manual, NOT_RUN. | values rounded to three decimals fail for the volcano. |
| V201 | PASS | grouped panel (four proteins, two groups, six replicates): rows are replicates, 8 data sets (protein by group), each value equals the source. R read-back. Prism opening: manual, NOT_RUN. | a missing protein data set fails the comparison. |

## Review fix to R16a (commit d2b8414)

The coordinator's review of R16a (16bb03e) found that the exact P label rounded across the 0.05 threshold (`0.04999` gave `p = 0.05`) and dropped trailing zeros (`0.01` gave `p = 0.01`). The fix is the separate commit **d2b8414** with amendment A-2026-10-01-31 (contract note in `figures.md`). The new V198 cases are negative cases: `evidence/fix-fail-before-pass-after.log` shows 14 of them failing on the old implementation and all 37 V198 tests passing on the fix. The R16a receipt, manifest, traceability rows and gate logs were updated in that commit.

**Open point for the Maintainer.** The V198 example values in `specs/017-visualization/spec.md` (`p = 0.04`, `p = 0.0004`) do not follow the corrected rule, which gives `p = 0.040` and `p = 0.00040`. The spec text is Maintainer-frozen and was not edited here.

## Notes and open points

- **Prism-openability is a manual Maintainer check**, recorded as **NOT_RUN** with the reason "manual Maintainer check in GraphPad Prism" for V199, V200 and V201 (`gate-results.json`, `manual_checks`).
- **`TwoWay` for grouped tables is an assumption to verify in Prism.** The pzfx package writes only `OneWay` and `XY` tables, so no reference file exists for the grouped code. The grouped layout (rows = replicates, one data set per protein and group) follows FR-201 as specified, not a Prism template.
- **Precision of the R read-back.** `pzfx::read_pzfx` parses numbers with R's `as.numeric`, so the R test compares the parsed values with `as.numeric` of the same strings. The exact 17-digit equality is verified in Python through `xml.etree`, which is the oracle the spec names.
- **Read-back precision caveat for the test code.** `as.character()` on a double keeps 15 significant digits, so the R test does not call it on the read values. Only `as.numeric()` is used.
- No R code was added to the package; only the testthat file and the DESCRIPTION Suggests entry changed.
- The `renv.lock` entries were generated from the installed DESCRIPTION files with R's `read.dcf`, and they follow the fields of the existing entries.
- No private data and no study output was used. All fixtures are synthetic.

## Acceptance

Every packet gate above exited 0, and every case, including its negative case, passed. R16b is recorded as ready for review (PASS by the Claude route; independent audit NOT_RUN). The Prism-openability checks are NOT_RUN as recorded.
