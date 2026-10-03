# R02 Maintainer-route receipt — Canonical protein input and biological identity

Packet R02 (FR-011–FR-020 / T011–T020 / V011–V020) was implemented and verified on 2026-10-01 (all gates re-executed 2026-10-02 after the Phase 2/3 amendments) by the **Claude (Opus) route**, which the operator authorized to act for the Sol/Astra/Luna roles ([ADR 0007](../../adr/0007-claude-route-and-integration-amendments.md)). This is not work by the routes named in AGENTS.md.

## Identities

- Contract versions: kit 1.2.0 with the 2026-10-01 amendments (ADR 0006/0007); scientific methods contract v1.2.0.
- Working-source manifest (packet allowlist, 32 files): `ec8a84f27729be9b6845b16dfe33b28bd6609214e212623a9ffca833f2472f4d` — lines in `evidence/working-source-manifest.txt` (`sha256  path`, sorted, LF).
- `verified_commit`: **null** — the packet was committed by Omar as `fd8dfa9`; the audit fixes and this re-verification are uncommitted on top of it (branch `claude/full-pipeline`).
- Environment: macOS-26.7-x86_64-i386-64bit-Mach-O; Python 3.13.15; R 4.6.1 with limma 3.68.5, statmod 1.5.2, impute 1.86.0, vegan 2.7-6, permute 0.9-10, DEqMS, proDA and fgsea at the versions pinned in `renv.lock`, in a project-local library (`.r-lib/`); `LANG=en_US.UTF-8`.

## Gates executed

| Command | Exit | Result | Log |
|---|---|---|---|
| `.venv/bin/python -m pytest tests/contract/intake tests/scientific/test_technical_units.py -q -rs -p no:cacheprovider` | 0 | 47 passed in 3.23s | `python-tests-1.log` |

Phase-level regression (R01 suite, historical unittest, all R tests, baseline registry, spec checker, whole Python suite) is in [`docs/validation/phase-gates`](../phase-gates/README.md).

## Acceptance cases

| Case | Status | Tests | Positive oracle | Negative case |
|---|---|---|---|---|
| V011 | PASS | test_canonical_intake.py::test_v011_* | wide and shuffled long give byte-identical matrix/masks; NA cells P07/U2, P08/T1 | duplicate long pair (row-located), duplicate header, missing long key; no publication |
| V012 | PASS | test_canonical_intake.py::test_v012_* | log2 0/−1 preserved; linear [1,2,4]→[0,1,2]; explicit zero → NA and counted | log2→log2 E_SCALE_SECOND_LOG before any transform (log2 spy never called); linear nonpositive E_LINEAR_NONPOSITIVE |
| V013 | PASS | test_technical_units.py::test_v013_* | 4 injections → 1 biological unit; paired example keeps 4 subjects, not 8 | unaggregated injections E_TECHNICAL_REPLICATION_UNMODELED |
| V014 | PASS | test_technical_units.py::test_v014_* | mean_linear = log2(10), mean_log2 = 3; distinct visit kept; coverage lineage | aggregation across visits/subjects E_TECHNICAL_AGGREGATION_INVALID |
| V015 | PASS | test_features_and_edges.py::test_v015_* | three opaque rows sharing accession/symbol; multi-gene members; unknown flags kept | duplicate feature_id with different annotation; invalid flag |
| V016 | PASS | test_vendor_mappings.py::test_v016_* | DIA-NN, MaxQuant, FragPipe, Spectronaut header-set profiles match hand-authored rows | peptide grain E_UNSUPPORTED_SCOPE; unknown header E_MAPPING_HEADER; version mismatch; undeclared zero encoding |
| V017 | PASS | test_legacy_workbook.py::test_v017_* | recovered parser wrapped unchanged (hash checked); sheet identities and Method A/B crosswalk | D1/D2 or T1/T2 swap E_CONTRAST_IDENTITY; mixed Method names E_CROSSWALK_METHOD/E_CROSSWALK_GROUP |
| V018 | PASS | test_canonical_intake.py::test_v018_*, test_r_bundle_reader.py | manifest hashes match independent hashlib; exact roundtrip in Python and R | rerun collision; changed source bytes E_SOURCE_CHANGED; deleted mask / edited value rejected by Python verify and the R reader |
| V019 | PASS | test_canonical_intake.py::test_v019_* | grain, scale, biological/technical n, missing cells, provenance gaps; no asserted unknown facts | unknown original mask stays unknown, never all-observed |
| V020 | PASS | test_features_and_edges.py::test_v020_* | quoted TAB + Unicode identifier roundtrip | ragged row, duplicated ID, mismatched mask, Inf/NaN, non-UTF-8, unsupported assay: each with file/row/field |


## Notes

- Decisions taken on ambiguous points are D-01–D-38 in [decisions.md](../../../specs/001-downstream-proteomics/decisions.md).
- Independent audit: an external audit of fd8dfa9 (2026-10-02) returned ACCEPT WITH FIXES; its findings are fixed (ADR 0008). Private-study regression: **NOT_RUN** (Maintainer-only, out of scope). No private data or private-derived numbers were used.

## Acceptance

All listed gates exited 0 and every acceptance case above, including its negative case, passed on re-execution. R02 is recorded as accepted by the operator-authorized Claude route (self-verified); the independent audit of fd8dfa9 returned ACCEPT WITH FIXES, and the fixes are self-verified and not yet re-audited.
