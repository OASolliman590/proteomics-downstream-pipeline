# Release/status checklist

Current state (2026-10-01): recovered baseline, frozen v1.2.0 Spec Kit with operator-authorized amendments, R01 accepted, and R02–R05, R10a and R13 passing their gates in an uncommitted working tree (Claude route, independent audit NOT_RUN). No maintained phase has shipped or been tagged. This checklist governs status/scope only; it does not rewrite the scientific audit.

- [x] Independently audit the actual candidate tree and resolve findings; record explicit Maintainer freeze before R01 (done 2026-09-20).
- [ ] Independently audit the R02–R05/R10a/R13 diff (not yet done; the 2026-10-01 work was self-verified by the Claude route).
- [ ] Verify protected pipeline/, legacy/, docs/audit/ bytes and recovered provenance hashes unchanged on the actual checkout.
- [ ] Accept only the phase whose tests actually executed: Phase 1=R01–R05/R10a; Phase 2=R06–R09; Phase 3=R10b/R11/R12. R10a never closes full R10 or SYS-03.
- [ ] Match README and GitHub description to the actual baseline/candidate/released phase. Do not describe planned paths, adapters, locks or calibration as shipped.
- [ ] Preserve NOT_RUN/unresolved for unavailable R, private archive, unchosen license and unexecuted gates; never invent evidence.
- [ ] Obtain owner decisions on license/ownership/resource terms and explicit public-release/disclosure authorization; exclude private data/credentials/logs and preserve originals.

Proposed GitHub description (once committed): **Recovered proteomics study baseline and a maintained successor in progress: Phase 1 limma core and PERMANOVA implemented on synthetic fixtures; not yet released or independently audited.** Updating repository metadata is a separate Maintainer action and was not performed by this uncommitted working tree. Until it is corrected, public-description consistency remains outstanding.

The maintained `proteomics` CLI implements validate/plan/run/verify/doctor for the Phase 1 scope and PERMANOVA; recovery usage is in [PIPELINE.md](PIPELINE.md). [PHASES.md](../specs/001-downstream-proteomics/PHASES.md) and [packet-index.md](../specs/001-downstream-proteomics/packet-index.md) are authoritative for future release/dispatch scope. The hardening report is [outside the protected audit tree](SPEC_HARDENING_REPORT.md).
