# Release/status checklist

Current state: recovered baseline plus frozen v1.2.0 Spec Kit. No maintained phase has shipped and no implementation dispatch is authorized without an explicit packet GO. This checklist governs status/scope only; it does not rewrite the scientific audit.

- [ ] Independently audit the actual candidate tree and resolve findings; record explicit Maintainer freeze before R01.
- [ ] Verify protected pipeline/, legacy/, docs/audit/ bytes and recovered provenance hashes unchanged on the actual checkout.
- [ ] Accept only the phase whose tests actually executed: Phase 1=R01–R05/R10a; Phase 2=R06–R09; Phase 3=R10b/R11/R12. R10a never closes full R10 or SYS-03.
- [ ] Match README and GitHub description to the actual baseline/candidate/released phase. Do not describe planned paths, adapters, locks or calibration as shipped.
- [ ] Preserve NOT_RUN/unresolved for unavailable R, private archive, unchosen license and unexecuted gates; never invent evidence.
- [ ] Obtain owner decisions on license/ownership/resource terms and explicit public-release/disclosure authorization; exclude private data/credentials/logs and preserve originals.

Proposed GitHub description: **Recovered proteomics study baseline and frozen v1.2.0 Spec Kit for a maintained successor; implementation has not started.** Updating repository metadata is a separate Maintainer action and was not performed by this uncommitted working tree. Until it is corrected, public-description consistency remains outstanding.

The future proteomics CLI is not currently implemented; present recovery usage is in [PIPELINE.md](PIPELINE.md). [PHASES.md](../specs/001-downstream-proteomics/PHASES.md) and [packet-index.md](../specs/001-downstream-proteomics/packet-index.md) are authoritative for future release/dispatch scope. The hardening report is [outside the protected audit tree](SPEC_HARDENING_REPORT.md).
