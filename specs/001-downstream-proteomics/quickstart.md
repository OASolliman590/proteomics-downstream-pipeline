# Quickstart status

The maintained CLI is not implemented at the recovered baseline. The commands below are **future Phase 1 acceptance targets**, not instructions claiming they work today:

```bash
proteomics doctor --json --config specs/001-downstream-proteomics/contracts/example-independent.json
proteomics validate --config specs/001-downstream-proteomics/contracts/example-independent.json --json
proteomics plan --config specs/001-downstream-proteomics/contracts/example-independent.json --output plan.json
proteomics run --config specs/001-downstream-proteomics/contracts/example-independent.json --output runs/example
proteomics verify --run runs/example --json
```

Paths inside each config resolve from its file directory. [Tiny fixtures](contracts/fixtures/README.md) are actual public synthetic specification examples, not private study data or fitted references. R01 creates only the real foundation; it must not pretend to fit these examples before R02–R05/R10a exist. Current recovery usage is in the root [README](../../README.md) and unchanged [pipeline documentation](../../docs/PIPELINE.md).

Read [START_HERE.md](START_HERE.md) and [IMPLEMENTATION_BRIEF_R01.md](IMPLEMENTATION_BRIEF_R01.md). Independent review and Maintainer freeze precede implementation. No slash-command installation or specify init --force is part of this workflow.
