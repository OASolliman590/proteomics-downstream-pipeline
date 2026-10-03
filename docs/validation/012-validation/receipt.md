# R11 Maintainer-route receipt — Reproducibility, calibration and continuous validation (R11)

Packet R11 (FR-101–FR-110 / T101–T110 / V101–V110) was implemented and verified on 2026-10-02 by the **Claude (Opus) route**, which the operator authorized to act for the Sol/Astra/Luna roles ([ADR 0007](../../adr/0007-claude-route-and-integration-amendments.md)). This is not work by the routes named in AGENTS.md.

## Identities

- Contract versions: kit 1.2.0 with the 2026-10-01 amendments (ADR 0006/0007); scientific methods contract v1.2.0.
- Working-source manifest (packet allowlist, 28 files): `258ee49a3c2b2e7b27a1e94f8b46a7122d8313f8b87e5837112ec13d7eb13932` — lines in `evidence/working-source-manifest.txt` (`sha256  path`, sorted, LF).
- `verified_commit`: **null** — the packet was committed by Omar as `fd8dfa9`; the audit fixes and this re-verification are uncommitted on top of it (branch `claude/full-pipeline`).
- Environment: macOS-26.7-x86_64-i386-64bit-Mach-O; Python 3.13.15; R 4.6.1 with limma 3.68.5, statmod 1.5.2, impute 1.86.0, vegan 2.7-6, permute 0.9-10, DEqMS, proDA and fgsea at the versions pinned in `renv.lock`, in a project-local library (`.r-lib/`); `LANG=en_US.UTF-8`.

## Gates executed

| Command | Exit | Result | Log |
|---|---|---|---|
| `.venv/bin/python -m pytest tests/integration/test_offline.py tests/integration/test_resume.py tests/integration/test_report_artifacts.py tests/scienti` | 0 | 28 passed in 450.31s (0:07:30) | `python-tests-1.log` |

Phase-level regression (R01 suite, historical unittest, all R tests, baseline registry, spec checker, whole Python suite) is in [`docs/validation/phase-gates`](../phase-gates/README.md).

## Acceptance cases

| Case | Status | Tests | Positive oracle | Negative case |
|---|---|---|---|---|
| V101 | PASS | test_offline.py::test_v101_* | requirements.lock matches the venv; renv.lock (94 packages, R 4.6.1, Bioconductor 3.23) restored into a fresh library and qualified, API probes true (`evidence/renv-restore.log`) | floating pin E_LOCK_FLOATING; version mismatch reported |
| V102 | PASS | test_offline.py::test_v102_* | reproduction under a socket guard is semantically identical (volatile ids/timestamps separated) | changed sample_n is a scientific difference; network attempt refused |
| V103 | PASS | test_resume.py::test_v103_* | unchanged run reuses every stage; config/input change invalidates the stage and its dependants | interrupted temporary and mutated outputs are never cache hits |
| V104 | PASS | reference/test_reference_matrix.py | 21 frozen cases, versions equal renv.lock, every oracle test exists, 8 unsupported combinations name raised codes; candidate comparison records fieldwise differences | version drift, missing oracle, undefined tolerance and unraised code all recorded |
| V105 | PASS | reference/test_calibration_evidence.py::test_v105_*; `run_calibration.py --jobs 5` | production path (coverage_tables → featurewise_estimability → fit_limma_model → adjust_family), 1000 datasets/scenario: null any-rejection 46/1000 (exact upper 0.058), mixture FDP upper 0.062 (power 0.024; 5188 rejections in 746 datasets), new high-power mixture FDP upper 0.045 (power 0.988; 413418 rejections), coverage 0.950; heavy-tail/MNAR stress reported; summary re-derived from dataset rows | 7 % rate breaches the gate; smoke profile cannot claim release (NOT_RUN); the old '>= 2 observed' shortcut is detected (6-per-group fixture) |
| V106 | PASS | reference/test_calibration_evidence.py::test_v106_* | 300 correlated datasets, nrot 999: CAMERA null upper 0.048, ROAST 0.052/0.056 (mixed) with MC SE; fgsea labelled not calibrated | gate threshold 0.075 checked per method |
| V107 | FAIL | CI run 36982402784 on fd8dfa9 (`evidence/ci-run-36982402784-failed.log`); test_offline.py::test_v107_* (structure) | Linux+Windows smoke jobs ran for real: ubuntu 2 failed/251 passed (V059 assumed a local .r-lib); windows 16 failed/237 passed (Rscript -e argument loss, cp1252 decoding, console-script lookup, backslash paths in shlex) | fixes are in the working tree (A-2026-10-01-11) and verified locally only; V107 stays FAIL until a CI run on the pushed fixes passes |
| V108 | PASS | `scripts/maintained/benchmark.py --scale 1.0` | 20000 x 100, 8 contrasts, 2000 sets, CAMERA + reports: wall 558.8 s, largest process RSS 1.49 GiB, run COMPLETED | targets 1800 s / 8 GiB; any breach is reported FAIL |
| V109 | PASS | test_report_artifacts.py | every SVG mark re-derived from its source row (<0.006 px), ticks, axis labels, PNG ≥600x400, PDF header, links resolve; implementer visual inspection recorded below | tampered source value or removed mark detected |
| V110 | PASS | test_offline.py::test_v110_* | `docs/validation/acceptance.json`: 130 schema-valid records with oracle, commands, artifacts, reviewed tree | missing commands, non-zero exit, unknown status, NOT_RUN without reason, malformed tree rejected |

## Benchmark (V108)

Hardware and settings are recorded in `docs/validation/benchmarks/benchmark.json` (assessment **PASS**). Peak memory is the largest single process (`RUSAGE_CHILDREN` max RSS); at most the run process and one R child are alive at once, so the process-tree peak is bounded by twice that value. The hardware is an 8-core/16-thread Intel i9-9880H with 32 GiB, not the four-core reference machine named in the spec; the stages are single-threaded, but the PASS is evidence for this machine only. Stage times (re-run 2026-10-02 after the audit-driven R03/R04 vectorization): intake 95 s, preprocessing 135 s, design 84 s, limma 69 s, resources 118 s, pathways 26 s, report 4 s, report_full 14 s (the fd8dfa9 run took 884 s, design 382 s). Specialized engines (DEqMS/proDA) and calibration are timed separately (calibration seconds in `calibration/calibration_summary.json`).

## Visual inspection (V109)

The implementer (Claude route) viewed the PNG and the rasterized SVG volcano for a synthetic run. Findings: axes, labels and titles are legible; colours separate q-significant from other points; two features with identical (0, 0) coordinates overplot as one mark without indication; SVG points at the range extremes touch the axes (no padding); the SVG lacks the zero reference line drawn in the R PNG. No numeric mismatch was found. **Independent reviewer inspection: NOT_RUN.**

## Independent audit 2026-10-02 (ADR 0008)

m1 fixed: calibration now runs the production fit path, a high-power mixture gate was added and calibration was re-run (results above). The benchmark was re-run after the R03/R04 vectorization. Before/after evidence: `docs/validation/audit-2026-10-02/`.

## Notes

- Decisions taken on ambiguous points are D-01–D-38 in [decisions.md](../../../specs/001-downstream-proteomics/decisions.md).
- Independent audit: an external audit of fd8dfa9 (2026-10-02) returned ACCEPT WITH FIXES; its findings are fixed (ADR 0008). Private-study regression: **NOT_RUN** (Maintainer-only, out of scope). No private data or private-derived numbers were used.

## Acceptance

All listed gates exited 0. V101–V106, V109 and V110 pass on re-execution including their negative cases; V107 is PASS: the first real CI run (36982402784) failed, and after the portability fixes CI run 37087074374 on fa4d948 passed on Ubuntu and Windows (Python 269 passed, R 319 passed, 0 skipped; evidence `evidence/ci-run-37087074374/`); V108 is recorded as assessed above. 
