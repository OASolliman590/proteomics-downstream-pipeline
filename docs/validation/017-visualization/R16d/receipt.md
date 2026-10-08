# R16d Claude-route receipt: standalone interactive HTML with vendored plotly.js (ADR 0012)

Packet R16d (FR-204, FR-205 / T204, T205 / V204, V205) was implemented on 2026-10-08 by the Claude route (Haiku builder, Opus review), at the operator's direction. This is not work by the routes named in AGENTS.md. The receipt means ready for review, not accepted. Independent audit: **NOT_RUN** (the Opus review follows this receipt).

## Identities

- Contract: `specs/017-visualization/contracts/figures.md` ("Vendored plotly.js", the error codes). Decisions: D-68, D-72 (not in this stage). The plan (`plan.md`) names the vendored files and `PINNED.json`.
- Working-source manifest: `evidence/working-source-manifest.txt`, 7 files, SHA-256 `35662ccd930b63c8b1befb425e7a59ee53930d9d677c1fba8e2f3fb61e01eac0`. It covers the writer, its test, the synthetic fixture, the vendored library, its LICENSE and PINNED.json, and the reason-code catalogue (A-2026-10-01-33).
- Shared-interface amendment: **A-2026-10-01-33** (`packet-ownership.json`): the reason-code catalogue (`E_PLOTLY_HASH`, `E_PLOTLY_HOVER`, `E_PLOTLY_NETWORK`), traceability, task ticks T204–T205, and the ledger.
- `verified_commit`: **null**. The R16d commit is local on `claude/visualization` and is not pushed.
- Environment: macOS, Python 3.13.15 (`.venv`), R 4.6.1 (the plotly package was installed in the project-local `.r-lib` only to obtain the vendored file). No R code changed in R16d.

## Vendored plotly.js (FR-205)

- `vendor/plotly/plotly.min.js`: plotly.js **2.25.2**, 3,575,010 bytes, SHA-256 `419deebe105f4c993feb7e6d1fe94fb1c2f7f98716db5d9b8fe3e37773838898`. It was copied from `htmlwidgets/lib/plotlyjs/plotly-latest.min.js` of the installed CRAN package plotly 4.12.1 in `.r-lib`. No network download was used.
- `vendor/plotly/LICENSE`: the MIT notice shipped with that package, unmodified. The library banner reads "Copyright 2012-2025, Plotly, Inc."; the LICENSE file reads "Copyright (c) 2021 Plotly, Inc". Both are kept as shipped.
- `vendor/plotly/PINNED.json`: version, source, SHA-256, licence, `vendored_by_packet: R16d`, and **approval: pending**.

## Implementation

- `src/proteomics_pipeline/figures/plotly_html.py`
  - `verified_library()` hashes the vendored file and compares it with the pin (`E_PLOTLY_HASH`).
  - `volcano_traces()` builds one scatter trace per group, in declared order. x is the effect, y is −log10 P, and `customdata` holds (feature, group, P, q). The hover template shows the feature identifier, the group, P and q.
  - `check_hover()` refuses a hover template that lacks the feature identifier or the group or P/q (`E_PLOTLY_HOVER`).
  - `render_html()` writes a standalone document: the MIT notice in the head, the vendored library inline, the figure JSON in a data element (`<`, `>` and `&` escaped), and the `Plotly.newPlot` call. No network reference.
  - `scan_html()` runs over the complete document before it is written.
- The figure JSON is written by Python, without the Python plotly package.

## Static scan: how the library block is handled

The minified library contains URL strings of its own (map tiles, documentation links, and the XML namespaces). The contract's allowlist applies to the document outside the library. The scan therefore does two things:

1. The library block (between `<script id="vendored-plotly-js">` and its closing tag) must equal the pinned file byte for byte (SHA-256 equal to the pin). Any other text in the block fails `E_PLOTLY_HASH`.
2. Outside the block, the document may contain no `<script src>`, no `<link href>`, no `fetch(`, no `XMLHttpRequest`, no dynamic `import(`, and every `http(s)://` string must start with one of `http://www.w3.org/2000/svg`, `http://www.w3.org/1999/xhtml` or `http://www.w3.org/1999/xlink`. Otherwise `E_PLOTLY_NETWORK`.

The library itself contains no `</script` and no `<!--`, so inlining it is safe.

## Gates executed

| Command | Exit | Result | Log |
|---|---|---|---|
| `.venv/bin/python -m pytest tests/scientific/test_plotly_html.py -q -rs -p no:cacheprovider` | 0 | 9 passed | `evidence/python-tests-plotly-html.log` |
| `.venv/bin/python scripts/check_spec_kit.py` | 0 | PASS | `evidence/spec-kit-check.log` |
| `.venv/bin/python -m pytest tests/regression tests/unit -q -p no:cacheprovider` | 0 | 85 passed | `evidence/regression-unit.log` |
| Full Python suite, `.venv/bin/python -m pytest tests -q -p no:cacheprovider` (background, once at the end) | 0 | 488 passed, 1 skipped in 30 min 52 s | `evidence/full-python-suite.log` |

## Acceptance cases

| Case | Status | Positive oracle | Negative case |
|---|---|---|---|
| V204 | PASS | the trace arrays parsed from the HTML figure JSON with `json.loads` equal the source rows of the synthetic volcano table: x = log2FC, y = −log10 P, customdata = (feature, group, P, q), one point per row. The hover template names the feature identifier, the group, P and q. | a hover template without the feature identifier fails `E_PLOTLY_HOVER` (direct and through the writer). |
| V205 | PASS | the SHA-256 of `vendor/plotly/plotly.min.js` equals the pin, and the MIT notice is in the vendored LICENSE and in the HTML. The static scan outside the library finds no external reference, and every URL is one of the three namespace strings. The library block equals the pinned file by hash. | a CDN script tag fails `E_PLOTLY_NETWORK`; an external `fetch` outside the library fails `E_PLOTLY_NETWORK`; a modified vendored file fails `E_PLOTLY_HASH` and no HTML is written. |

## Open points for the Maintainer

- **Pin approval is pending.** `contracts/figures.md` says `PINNED.json` is finalised when the Maintainer approves the pin, and that R16d is NOT_RUN until then. The receipt records the automated checks as PASS by the Claude route, but the pin's `approval.status` remains `pending`.
- **Version.** plotly.js 2.25.2 is the version bundled in the installed plotly R package 4.12.1. The plan asks for the latest stable release. No network access was used to check it, so the Maintainer must confirm the version or supply a newer file.
- **Browser check.** The HTML was checked by its byte structure and by the scan. It was not opened in a browser in this packet, so the rendering and the hover behaviour are not yet verified visually.
- No private data was used. The fixture is synthetic.

## Notes on the full suite

- The full run ended before the R16c review fixes were applied (they are in a later commit), so it tests the R16d tree.
- One test is skipped in the full run. The log does not name it (no `-rs` was used for that run); the R16d-specific tests are not skipped.

## Acceptance

All gates exited 0, and every case, including its negative case, passed. R16d is recorded as ready for review (PASS by the Claude route; independent audit NOT_RUN; the pin awaits Maintainer approval).
