# Canonical pipeline rework checklist

Updated 2026-11-09. This is the current 13-category checklist; older checkpoint
narratives are acceptance evidence, not current instructions to repeat investigations.
The supported core is defined by [interface 1.0](../docs/core-interface.md).

| Category | Status | Evidence / next step |
| --- | --- | --- |
| 1. Assembly/QC acceptance | Accepted | Ten finalized assemblies; final QC, Hi-C retention/cis-trans and contact maps reconciled in the accepted integration |
| 2. General-purpose inputs | Accepted | HiFi, phased/collapsed shotgun, identity/ploidy/species isolation and reasoned skips; [identity acceptance](core-identity-reference-01-assessment.md) |
| 3. Failure and organelle handling | Accepted | Controlled recovery/failure and native HiFi/shotgun routes; [shotgun acceptance](core-shortread-native-01-completion-assessment.md) |
| 4. Shared references/provenance | Accepted | Reference selection independent of graph admission; final identity fingerprints; observed tool/environment provenance |
| 5. Core pangenome | Accepted | [Production graph acceptance](core-graph-production-02-assessment.md) and [tool-change rationale](../docs/pangenome-tool-changes.md); five known representation flags retained |
| 6. Unified Markdown reporting and freeze | Accepted | [Core-freeze-01](core-freeze-01-assessment.md): native PASS, 875 cached/64 completed tasks, 543 report links, zero version-record issues and unchanged selected evidence |
| 7. Operational cleanup | Accepted | [Combined integration](core-integration-02-assessment.md); ordinary Nextflow cache handling and exact YaHS replay |
| 8. PSMC | Accepted | Five individuals, 55 verified fits, core integration and read-support audits passed. Auto routes paired haplotypes to Dipcall and a sole assembly from a diploid individual to HiFi. Method-history disagreement is a deferred open question, not an implementation blocker; [assessment and decision](psmc-methods-and-cohort-acceptance-2026-11-09.md) |
| 9. BlobToolKit evidence / cleaning | Planned | Rebuild evidence first; separately validate explicit cleaning and its downstream invalidation |
| 10. Graph QC | Planned | Bounded read-to-graph self-mapping plus separately scheduled missing intrinsic graph metrics; [agreed QC/MultiQC scope](pangenome-qc-growth-followup.md) |
| 11. Private/accessory evidence | Planned | Distinguish individual/haplotype privacy, missing representation and detection power |
| 12. Graph SV interpretation / supplementary SyRI | Planned | Known-answer method validation and bounded pilot; active SyRI route remains a pilot |
| 13. Additional graph products / presentation | Planned | Restore pangenome MultiQC; haplotype/individual growth fits, open/closed summary and confidence-band exploration; [agreed scope](pangenome-qc-growth-followup.md) |

## Core freeze and cleanup

The main/detailed reports, both chimera evidence figures, per-library junction maps,
sharing/representation figures and input-aware missingness are integrated. Four cuts
and two retained joins remain verified and unchanged. No further join investigation
or biological purging optimization is required by the core reporting gate.

The interface/documentation freeze removes obsolete implementations and their tests,
report mockups and investigation plans, without deleting scientific outputs. General
launchers no longer embed specimen-specific settings. Integration tests explicitly
recover supported settings from a supplied accepted run's recorded provenance.

Cooler/tidk version producers now emit labelled single-line records; the real
Nextflow producer-statement/collector test passes with zero issues. Existing cached
files are immutable. The [native core-freeze check](core-freeze-01-assessment.md) confirmed the new
producer records with zero issues, without rebuilding assemblies or graphs. All
84 selected baseline evidence files were byte-identical.

## Extension requirements retained

The [optional analysis scope](post-core-optional-analyses-plan.md) and
[PSMC design](completed-run-review-and-psmc-plan.md) remain authoritative for future
work, subject to the completed core contracts above. They are designs, not claims
that a flag already enables a working implementation.

PSMC is first: support paired same-individual haplotypes and one assembly plus
complete HiFi reads. Validate callability/masks, fitting, independent scheduled
bootstraps and mutation-rate/generation-time scaling. Preserve species and individual
identities. It must work independently of graph construction.

Each extension is independently enabled, has an explicit unsupported/skip/failure
policy, contributes evidence and report outputs, and requires relevant known-answer
and native acceptance. Automatic chimera cutting remains deferred until explicitly
resumed. Retired code remains in Git history. At the user's request, the variant/
rearrangement, private/accessory, progressive-graph and 2D-layout prototypes are
retained in their original source locations for review and reuse at categories
11-13. They remain uncalled and require validation before production integration;
this does not create a second supported interface.

The next pangenome work includes the [agreed graph-QC, MultiQC and growth-fit
additions](pangenome-qc-growth-followup.md). Growth fitting/open-closed summaries
are requested deliverables, not undecided functionality. Confidence-band exploration is also explicitly included in that visit, covering
observed-cohort variability and fitted-model uncertainty. Any published band must
use a defined, validated interpretation; the old independent-base band is not retained.
