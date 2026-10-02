# Example configurations

| File | What it shows |
|---|---|
| `example-independent.json` | Phase 1: three independent groups, limma zero-null, primary and secondary families. |
| `example-paired.json` | Phase 1: paired before/after design with fixed subject terms. |
| `example-effect-threshold.json` | Phase 1: TREAT (effect-threshold) primary hypothesis. |
| `../../tests/fixtures/reports/demo-phase2.json` | Phase 2 demonstration: limma + descriptive response + PERMANOVA on the same synthetic data. |

All inputs under `fixtures/` are invented (see `fixtures/README.md`). Paths in a configuration are resolved relative to the configuration file.
