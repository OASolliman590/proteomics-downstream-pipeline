# R10a report fixtures (synthetic)

The report tests build their runs with `tests/fixtures/design/builders.py`
(read-only reuse) so the report is exercised on real Python -> R stage
manifests: a zero-discovery run, a QC-only run, a run whose required limma
stage genuinely fails (every eligible feature overflows double precision),
and a technical-replicate run with a declared exclusion and unknown upstream
metadata. No study data are used.
