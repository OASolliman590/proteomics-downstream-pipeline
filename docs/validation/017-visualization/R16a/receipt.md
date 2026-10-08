# R16a Claude-route receipt: style system, figure registry and `report.figures` configuration (ADR 0012)

Packet R16a (FR-195–FR-198 / T195–T198 / V195–V198) was implemented on 2026-10-08 by the Claude route (Haiku builder, Opus review), at the operator's direction. This is not work by the routes named in AGENTS.md. The receipt means ready for review, not accepted. Independent audit: **NOT_RUN** (the Opus review follows this receipt).

## Identities

- Contract: `specs/017-visualization/contracts/figures.md` (configuration table, P-value rule, figure IDs, error codes). Decisions: D-68, D-72 (not in this stage), D-74 (120 mm default and the Okabe–Ito order).
- Working-source manifest: `evidence/working-source-manifest.txt`, 10 files, SHA-256 `9a7f826d366a8fe64b4e893e5207ff4dcf0c8fb71c2610dc0a43707fefa5a2ec`. It covers the packet files and the two schema copies and the reason-code catalogue that A-2026-10-01-29 names. The ledger, traceability, tasks and this receipt are excluded because they record the result.
- Shared-interface amendment made first: **A-2026-10-01-29** (`packet-ownership.json`). It covers the `report.figures` block in both `analysis.schema.json` copies, the traceability rows, the task checkboxes T195–T198, the regenerated ledger and the regenerated reason-code catalogue.
- `verified_commit`: **null**. The R16a commit is local on `claude/visualization` and is not pushed.
- Environment: macOS, Python 3.13.15 (`.venv`), Rscript 4.6.1, `LANG=en_US.UTF-8`. No R code was changed in R16a, so the R package was not reinstalled.

## Gates executed

| Command | Exit | Result | Log |
|---|---|---|---|
| `.venv/bin/python -m pytest tests/unit/test_figure_config.py -q -rs -p no:cacheprovider` | 0 | 16 passed | `evidence/python-tests-unit.log` |
| `.venv/bin/python -m pytest tests/scientific/test_figure_style.py -q -rs -p no:cacheprovider` | 0 | 26 passed | `evidence/python-tests-scientific.log` |
| `.venv/bin/python scripts/check_spec_kit.py` | 0 | PASS | `evidence/spec-kit-check.log` |
| `.venv/bin/python -m pytest tests/regression tests/unit -q -p no:cacheprovider` | 0 | 85 passed | `evidence/regression-unit.log` |

The regression and unit run includes the reason-code catalogue test (V115), the traceability reconciliation (V116) and the text-encoding lint.

## Acceptance cases

| Case | Status | Positive oracle | Negative case |
|---|---|---|---|
| V195 | PASS | the contract defaults written in the test (formats png, pzfx, html; dpi 300; width 120 mm; Okabe–Ito; Arial; exact). Absent and explicit-default configurations resolve to them, and the values are written to `registry.json` (`style`). | dpi 150 fails `E_FIGURE_DPI`; `formats: ["jpg"]` fails `E_FIGURE_FORMAT`; an unknown key fails `E_CONFIG_SCHEMA` (through `load_config` and the resolver). No figure or registry is written. |
| V196 | PASS | Okabe–Ito hex values and the declared group order written in the test. Three figure types in two orders get the same colour for each group. | an order-dependent assigner fails `E_FIGURE_COLOUR_UNSTABLE`. |
| V197 | PASS | the three registered catalogue IDs; the SHA-256 of each source table computed with `hashlib`. Each `sources/<id>.tsv` and every output (png, pzfx, html) equals its registered source. | a request without a source table fails `E_FIGURE_SOURCE_MISSING`; the renderer is never called and no file is written. A missing or changed source is refused. |
| V198 | PASS | the synthetic fixture `tests/fixtures/figures/style/p_values.tsv` (0.2, 0.04, 0.0004, 1e-7): stars ns, *, \*\*\*, \*\*\*\*; exact p = 0.20, p = 0.04, p = 0.0004, p < 0.0001. Thresholds at the contract boundaries. | stars shown in exact mode fail; a missing, NaN or out-of-range P gives `E_FIGURE_SOURCE_MISSING`, so no annotation appears without its P value. |

## Notes and decisions

- **Exact P rule.** Two significant digits, plain decimal, trailing zeros removed and then padded to two decimals. This reproduces the contract examples (0.20, 0.04, 0.0004). The floor is `p < 0.0001`.
- **Run settings.** The resolved style is written to `<run>/report/figure_catalogue/registry.json`, under `style`. `config.resolved.json` is written by `workflow.py`, which is outside the R16a allowlist, so the figure resolver is not yet called from `load_config` or the workflow. Wiring belongs to R16e.
- **Unknown-key and typed codes.** The contract names no code for width, font or annotation values, so those use `E_CONFIG_SCHEMA`. The schema accepts any integer `dpi` and any string format, so that `E_FIGURE_DPI` and `E_FIGURE_FORMAT` come from the resolver, as the negative cases require.
- **Registry.** Each source table is copied to `sources/<id>.tsv` and hashed. The original path is not recorded, so no private path enters the registry.
- **No R file.** `r/proteomicsCore/R/figure_style.R` (on the allowlist) was not created. R16a needs no R code: the palette and annotation rule are in Python, and the ggprism theme belongs to R16c.
- **Traceability text.** The FR-195 criterion in `traceability.json` still says `journal_width_mm` default 85 (the pre-D-74 wording). The spec and the contract say 120 mm. The traceability text was not edited; the Maintainer should reconcile it.
- **Reason codes.** `docs/methods/reason-codes.md` was regenerated with `scripts/maintained/build_release.py codes`. The generator lists `config.py` twice for `E_CONFIG_SCHEMA` because two files share that basename. This is cosmetic.
- No private data and no study output was used. All fixtures are synthetic.

## Acceptance

All four gates exited 0, and every case above, including its negative case, passed. R16a is recorded as ready for review (PASS by the Claude route; independent audit NOT_RUN). The ledger was regenerated with `write_acceptance_ledger()`.
