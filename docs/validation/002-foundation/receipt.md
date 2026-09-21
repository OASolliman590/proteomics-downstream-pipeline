# R01 Maintainer acceptance receipt

Packet R01 (foundation, FR-001–FR-010 / T001–T010 / V001–V010) was independently verified on 2026-09-21 in an isolated Linux environment. This receipt records foundation behavior only; it does not claim intake, QC, design, inference, reporting, private-data validation, licensing, or public-release readiness.

## Identities

- Frozen contract manifest: `0619003741ee23f93ed08c332c1bfce3e21b446f277073be7010fa68286a6ffa`.
- Prerequisite working-source manifest: `84149fd7526e07346ba914d349939f9945cd3f0f0e6653fc1591f8b807101120`.
- Accepted R01 source manifest: `ecfaf907a0e7194977729900642850e953dd2142d2d03eeea300c61e471da5f3` over 51 existing R01-allowlisted source, test, fixture, bootstrap, and stable validation files. The mutable `receipt.md` and `solved-versions.json` evidence records are excluded; lines are `sha256  repository-relative-path`, sorted lexicographically with LF termination.
- `verified_commit`: `bd64c42a3c652beb22986b5e2c8f9114542362a1`. This is the reviewed R01 source commit; the later evidence-only commit records this identity without changing the reviewed source manifest.

## Maintainer verification

- Editable installation completed in the project-local Python environment.
- Both CLI entry points reported `0.1.0.dev0`; help and schema-only validation succeeded.
- Environment doctor reported Python 3.13.15, R 4.5.3, `jsonlite` 2.0.0, `openssl` 2.4.2, `testthat` 3.3.2, and `proteomicsCore` 0.1.0 as available.
- `python -m unittest discover -s tests -v`: 10/10 historical tests passed unchanged.
- `python -m pytest tests/contract/foundation tests/unit/foundation -q -rs`: 69 passed, 0 skipped.
- `Rscript --vanilla scripts/maintained/test_r.R --suite foundation`: PASS with no warnings, including real UTF-8/NA I/O, literal metacharacter paths, authentic failure capture, unavailable dispatch, and connection lifecycle checks.
- Baseline registry verification returned zero protected-tree mismatches.
- Full `plan` and `run` remained unavailable with typed exit 3; missing-run verification exited 5; removing R from PATH made doctor report foundation NOT_AVAILABLE and exit 3.
- A clean copied checkout ran `scripts/maintained/bootstrap.sh` successfully using project-local Python and R libraries, then repeated the Python, R, baseline, doctor, and solved-version gates.
- Solved-version evidence contains redacted executable paths and no personal or private-study data.

## Independent review

- The bounded cloud-derived source/test delta received an Astra audit. Two P2 harness defects were corrected: the negative now exercises a real failing `testthat` assertion through the maintained runner, and child script/input paths are quoted and exercised under spaces/metacharacters.
- A fresh Astra re-audit found no gaps within the corrected harness scope; runtime verification remained the Maintainer's responsibility and was rerun afterward.

## Acceptance

`V001 PASS; V002 PASS; V003 PASS; V004 PASS; V005 PASS; V006 PASS; V007 PASS; V008 PASS; V009 PASS; V010 PASS.`

R01 is accepted as foundation only. R02 and every later packet remain unstarted and require separate authorization. Private regression, license/ownership, and public-release gates remain unresolved outside R01.
