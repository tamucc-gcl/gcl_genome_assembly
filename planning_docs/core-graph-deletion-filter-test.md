# Cached deletion-filter correction test — 2026-10-23

The filter log confirms that `filter-paf-deletions` removed 59,116 and 62,194
alignment records across the missing chr7 portions of the retained CTlk composites.
The v1.1.8 implementation evaluates a dominant reference path but can delete other
chromosome mappings inside the chosen record range without including their support
in its removal cost. Assembly cuts and retained joins remain unchanged.

## Removable correction and scope

`--pangenome_deletion_filter_replay true` enables a separate diagnostic workflow.
It consumes the existing SV rGFA, exported PAF, pre-deletion PAF and filter log from
normal Nextflow construction outputs. It does not change construction code or
production graph policy. There are no dataset names, chromosome pairs or biological
break coordinates in the diagnostic implementation.

The `chromosome-isolated-v1` guard wraps the existing pinned executable; it is not
a rebuilt upstream binary or a claim of an identical minimal upstream patch.
It separates each source query's mappings into chromosome-specific virtual lanes,
applies native deletion filtering to each lane, and restores original names.
Original source lengths, query/target coordinates and alignment tags remain intact.
It also filters within non-dominant chromosome lanes, a policy difference that must
be assessed in the comparison rather than assumed harmless.

Chromosome membership comes from the SV graph's rank-zero reference paths. Only
non-reference nodes are connected when finding reference-boundary components;
reference nodes are boundaries, so a link between reference chromosomes does not
merge their identities. Used nodes bordering conflicting reference chromosomes or
having no reference anchor cause an explicit error rather than a guessed label.

The guard version/source SHA256, Cactus commit, executable hashes, effective XML,
lane ledger and filtering commands are retained. This branch can be removed without
changing the ordinary graph workflow if a verified upstream fix replaces it.

## Test matrix and cache contract

One preparation task converts the cached graph for the pinned native filter and
builds node ownership. Three filtering cases then run independently:

| Case | Deletion filtering |
|---|---|
| native | Pinned executable and unchanged thresholds |
| chromosome_isolated | Same executable/thresholds, separated chromosome lanes |
| disabled | Bypass deletion filtering only; retain subsequent PAF filters |

All cases receive the complete original pre-deletion PAF. Keeping its original
ordering avoids changing equal-query-start sort ties in the native control.
Native replay must reproduce the complete cached exported PAF byte for byte, or
the batch fails before accepting the comparison. After deletion filtering, whole
selected queries are retained for the common pinned PAF filter and assignments.
Selection includes missing sources, composites, low-retention chromosome paths
and ordinary represented chromosome controls using the existing generic screen.

Each filtered result then runs whole-query and 1 Mb regional assignment in separate
tasks: six assignment tasks across the three filtering cases. Large mappings,
graph conversion products and node ownership tables stay in the work cache.
Compact reports, intervals, commands, lane ledgers and logs are published and
included in the review archive. The source assemblies are never split by this test.

## Run on Crest after syncing

From the project directory, with no other controller using this resume session:

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_chimera_pipeline.sbatch \
  data/assembly_samplesheet.csv data/hic_readsets.csv \
  comparisons/core-graph-deletion-01 \
  --chimera_break "$PWD/data/chimera_reviews/Sde-CTlk_104.v8.reviewed.tsv" \
  --chimera_hifi_bam_manifest "$PWD/comparisons/chimera-observable-contact-v7/results/assembly/chimeras/mapping_manifest.json" \
  --run_pangenome true --pangenome_validate_only false \
  --pangenome_assignment_replay false \
  --pangenome_deletion_filter_replay true
```

Normal Nextflow resume should reuse assembly, Cactus construction and existing
graph summaries. Only this diagnostic branch and affected report consumers need
new work. No compiler, replacement container or new conda environment is required.

Start with `results/pangenome/373251/deletion_filter/deletion_filter_report.md` and
return `core-graph-deletion-01.review.tar.gz`.

## Acceptance before production integration

1. Native replay exactly matches the existing exported PAF.
2. The guard preserves the missing chromosome portions through the subsequent
   quality/overlap filters, rather than merely increasing an assigned span.
3. Regional assignment places both recovered portions on credible chromosome
   groups; inspect assigned regions and surviving mapping intervals.
4. Review every ordinary-control decrease and cohort-wide source decrease. Compare
   the disabled case to establish which problematic mappings the guard still removes.
5. Only after accepting this batch, select the production integration and perform
   one graph rebuild. Verify final FULL/CLIP representation and downstream summaries;
   intermediate recovery alone is not final graph acceptance.

Local validation uses focused Python regression tests plus the actual diagnostic
Nextflow workflow with tiny external tool doubles. A two-species construction
fixture checks the graph/pre-deletion/exported-PAF/log channels. These validate
calculation, routing and publication; real pinned-tool behavior on the biological
inputs remains the purpose of the Crest batch.
