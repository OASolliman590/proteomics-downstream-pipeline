# Roadmap — v1.2.0-frozen

The twelve original slices and all 120 FR/T/V identities remain. [PHASES.md](PHASES.md) is the release-scope authority; [packet-index.md](packet-index.md) is the dispatch authority and links the exact disjoint allowlists. The machine-readable [roadmap](roadmap.json) records full-slice completion dependencies, including R10 subpackets.

```text
Phase 0: baseline-recovered (current study code)
Phase 1: R01 → R02 → R03 → R04 → R05 → R10a
Phase 2:                         R06 ∥ R07 ∥ R09 → barrier → R08
Phase 3: R10a + R08 + R09 → R10b → R11 → R12
```

The parallel group starts only after R05 and serial R10a are accepted. R08 depends scientifically on R06/R07 and has R09 as an explicit scheduling dependency, so it starts only after the whole group is accepted. Full R10 is complete only after R10b; R11 depends on full R10. Phase names do not change eligibility requirements or turn absent adapters into optional successes.

No phase of the maintained successor is shipped. Independent candidate audit and Maintainer freeze precede the [R01 brief](IMPLEMENTATION_BRIEF_R01.md). [Progress](progress.md) and [traceability](traceability.json) hold actual state, not anticipated implementation evidence.
