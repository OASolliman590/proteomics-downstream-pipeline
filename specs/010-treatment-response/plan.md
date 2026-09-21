# Implementation plan — R09

**Phase:** 2. **Status:** frozen; dispatch requires Phase 1 acceptance and separate authorization.

Read [spec.md](spec.md), [tasks.md](tasks.md), the [shared execution contract](../001-downstream-proteomics/contracts/cli-and-artifacts.md), and the exact [ownership record](../001-downstream-proteomics/packet-ownership.json). No shared file may be reopened by this packet. Dependency dispatch and release scope are controlled by [packet-index.md](../001-downstream-proteomics/packet-index.md), not by a worker's interpretation of this plan.

Use each V case's independent fixture/oracle before implementing its behavior. Keep fixture authorship/reference calculations separate from the production adapter; do not import the adapter into its own oracle. Consume preceding packet artifacts through the typed interfaces and expose capabilities()/execute(request) only when the owned behavior really exists. All planned source/test/evidence paths are creation targets, not current files.

The per-task order below is serial unless the packet index explicitly permits independent packets. Keep scientific eligibility, tests and result types unchanged. If an interface cannot satisfy a requirement, return a concrete contract blocker to the Maintainer; do not change a frozen spec or another packet's files. Record actual command, expected/observed result, version, input/output hashes and exit status. The Maintainer reruns and reviews evidence before marking any task complete.

Acceptance is all assigned V cases, including negative cases, plus the phase's cross-cutting scenarios. Missing R or private evidence is NOT_RUN. Returning a receipt, producing a directory or showing green Python-only tests is not completion.
