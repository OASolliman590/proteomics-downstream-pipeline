# Validation record

Executed on 2026-09-12 with Python 3.12.14, NumPy 2.3.5 and pandas 3.0.1.

| Check | Result |
|---|---|
| Original archive | 165 files; 59,355,573 uncompressed bytes; outer manifest 164/164 entries matched |
| Inner Method B manifest | 138 present entries match; two absent .DS_Store entries are stale metadata |
| Workbook | Four sheets, 3,714 rows each; 20 sample columns reconstructed; all observed/missing cells match exported matrix; eight missing drawing/VML relationships |
| Method A rerun | Executed successfully; 9 reproduced data/result/validation tables match archive within 1e-12 |
| Source adjusted P values | All four sheets' BH values agree within 1e-10 |
| Method A/Method B matrices | Equal; no numerical or missingness change |
| Prepared Method B inputs | All three restored inputs match source tables, including metadata; 1,345 missing values retained |
| Refit differential results | All 15 contrasts: mean-based logFC, BH adjustment and formal flags agree; maximum effect error below 2.2e-14 |
| Pathway adjustment | All 24 collection/contrast families match BH within 1e-10 using finite-P denominator; one nonfinite pathway P remains a reported finding |
| Reversal | Ratios and class counts reconstructed; exhaustive fixed-score arithmetic exposes +1 discrepancy; these P values remain selection-biased |
| Tabular shape | 75 archived CSV/TSV tables parsed with consistent row width |
| Python sources | 6 files parsed successfully; pipeline/legacy/runner/tests compileall also passed |
| Synthetic tests | 10 passed, covering missingness, IDs, groups, treatment labels, annotation alignment, infinity and no-overwrite/runtime failure |
| Source scan | No credential-like pattern or personal absolute path found in proposed source package; pattern scan is not a security guarantee |
| R preflight | Input validation passed; Rscript unavailable; full R analysis and R syntax/runtime tests NOT RUN |
| GitHub CI | Workflow prepared; NOT RUN |

The first unit-test attempt could not write temporary files in its execution environment. A compatible rerun passed all ten tests. The initial failure was environmental, not a passing test or a scientific result.

See `data_checks.json`, `validation_record.json`, `table_inventory.json` and `source_provenance.json` for structured evidence. The original study results remain immutable in the sibling private_evidence folder. No limma moderated statistics or fgsea enrichment was recomputed; observed output consistency must not be described as a full numerical reproduction.
