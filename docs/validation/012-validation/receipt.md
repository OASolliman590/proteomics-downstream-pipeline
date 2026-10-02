# R11 Maintainer-route receipt — Reproducibility, calibration and continuous validation (R11)

Packet R11 (FR-101–FR-110 / T101–T110 / V101–V110) was implemented and verified on 2026-10-02 by the **Claude (Opus) route**, which the operator authorized to act for the Sol/Astra/Luna roles ([ADR 0007](../../adr/0007-claude-route-and-integration-amendments.md)). This is not work by the routes named in AGENTS.md.

## Identities

- Contract versions: kit 1.2.0 with the 2026-10-01 amendments (ADR 0006/0007); scientific methods contract v1.2.0.
- Working-source manifest (packet allowlist, 26 files): `8b5b78ec9480e027b52802b5e97dae59efb1b62be17e64ecb457c397ef30f675` — lines in `evidence/working-source-manifest.txt` (`sha256  path`, sorted, LF).
- `verified_commit`: **null** — nothing is committed (branch `claude/full-pipeline`, uncommitted working tree).
- Environment: macOS-26.7-x86_64-i386-64bit-Mach-O; Python 3.13.15; R 4.6.1 with limma 3.68.5, statmod 1.5.2, impute 1.86.0, vegan 2.7-6, permute 0.9-10, DEqMS, proDA and fgsea at the versions pinned in `renv.lock`, in a project-local library (`.r-lib/`); `LANG=en_US.UTF-8`.

## Gates executed

| Command | Exit | Result | Log |
|---|---|---|---|
| `.venv/bin/python -m pytest tests/integration/test_offline.py tests/integration/test_resume.py tests/integration/test_report_artifacts.py tests/scienti` | 0 | 25 passed in 88.62s (0:01:28) | `python-tests-1.log` |

Phase-level regression (R01 suite, historical unittest, all R tests, baseline registry, spec checker, whole Python suite) is in [`docs/validation/phase-gates`](../phase-gates/README.md).

## Acceptance cases

| Case | Status | Tests | Positive oracle | Negative case |
|---|---|---|---|---|
| V101 | PASS | test_offline.py::test_v101_* | requirements.lock matches the venv; renv.lock (94 packages, R 4.6.1, Bioconductor 3.23) restored into a fresh library and qualified, API probes true (`evidence/renv-restore.log`) | floating pin E_LOCK_FLOATING; version mismatch reported |
| V102 | PASS | test_offline.py::test_v102_* | reproduction under a socket guard is semantically identical (volatile ids/timestamps separated) | changed sample_n is a scientific difference; network attempt refused |
| V103 | PASS | test_resume.py::test_v103_* | unchanged run reuses every stage; config/input change invalidates the stage and its dependants | interrupted temporary and mutated outputs are never cache hits |
| V104 | PASS | reference/test_reference_matrix.py | 21 frozen cases, versions equal renv.lock, every oracle test exists, 8 unsupported combinations name raised codes; candidate comparison records fieldwise differences | version drift, missing oracle, undefined tolerance and unraised code all recorded |
| V105 | PASS | reference/test_calibration_evidence.py::test_v105_*; `run_calibration.py` | 1000 datasets/scenario: null any-rejection 46/1000 (exact upper 0.058), mixture FDP upper 0.062, coverage 0.950; heavy-tail/MNAR stress reported; summary re-derived from dataset rows | 7 % rate breaches the gate; smoke profile cannot claim release (NOT_RUN) |
| V106 | PASS | reference/test_calibration_evidence.py::test_v106_* | 300 correlated datasets, nrot 999: CAMERA null upper 0.048, ROAST 0.052/0.056 (mixed) with MC SE; fgsea labelled not calibrated | gate threshold 0.075 checked per method |
| V107 | NOT_RUN | test_offline.py::test_v107_* (structure only) | `.github/workflows/maintained.yml`: Linux+Windows smoke (renv restore, pytest, testthat, synthetic run) and dispatch-only release calibration | jobs have never run on GitHub (no push authorized) |
| V108 | PASS | `scripts/maintained/benchmark.py --scale 1.0` | 20000 x 100, 8 contrasts, 2000 sets, CAMERA + reports: wall 883.7 s, largest process RSS 1.62 GiB, run COMPLETED | targets 1800 s / 8 GiB; any breach is reported FAIL |
| V109 | PASS | test_report_artifacts.py | every SVG mark re-derived from its source row (<0.006 px), ticks, axis labels, PNG ≥600x400, PDF header, links resolve; implementer visual inspection recorded below | tampered source value or removed mark detected |
| V110 | PASS | test_offline.py::test_v110_* | `docs/validation/acceptance.json`: 130 schema-valid records with oracle, commands, artifacts, reviewed tree | missing commands, non-zero exit, unknown status, NOT_RUN without reason, malformed tree rejected |

## Benchmark (V108)

Hardware and settings are recorded in `docs/validation/benchmarks/benchmark.json` (assessment **PASS**). Peak memory is the largest single process (`RUSAGE_CHILDREN` max RSS); at most the run process and one R child are alive at once, so the process-tree peak is bounded by twice that value. The hardware is an 8-core/16-thread Intel i9-9880H with 32 GiB, not the four-core reference machine named in the spec; the stages are single-threaded, but the PASS is evidence for this machine only. A few short test runs overlapped the benchmark (minor contention, which can only lengthen the time). Stage times: design 382 s, preprocessing 224 s, resources 99 s, intake 80 s, limma 49 s, pathways 22 s, reports 17 s. Specialized engines (DEqMS/proDA) and calibration are timed separately (calibration seconds in `calibration/calibration_summary.json`).

## Visual inspection (V109)

The implementer (Claude route) viewed the PNG and the rasterized SVG volcano for a synthetic run. Findings: axes, labels and titles are legible; colours separate q-significant from other points; two features with identical (0, 0) coordinates overplot as one mark without indication; SVG points at the range extremes touch the axes (no padding); the SVG lacks the zero reference line drawn in the R PNG. No numeric mismatch was found. **Independent reviewer inspection: NOT_RUN.**

## Notes

- Decisions taken on ambiguous points are D-01–D-30 in [decisions.md](../../../specs/001-downstream-proteomics/decisions.md).
- Independent audit (Astra-equivalent): **NOT_RUN**. Private-study regression: **NOT_RUN** (Maintainer-only, out of scope). No private data or private-derived numbers were used.

## Acceptance

All listed gates exited 0. V101–V106, V109 and V110 pass on re-execution including their negative cases; V107 is NOT_RUN (no CI job has run); V108 is recorded as assessed above. 
