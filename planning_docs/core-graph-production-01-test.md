# Production graph sequence-policy acceptance

The supported chunk comparison selected chromosome-isolated deletion filtering,
the pinned regional grouping correction and a 25 kb regional chunk. These are now
the normal graph-construction path, independent of all optional replay/trace flags.
The implementation contains no sample, scaffold, chromosome-name or breakpoint rules.

## What the production run changes

- Build the checksum-verified corrected `rgfa-split` once in a separate cached task.
  A native CLI regression must pass before graph construction can start.
- Use task-local wrappers around Cactus 3.2.1 preprocessing. Deletion filtering
  compares mappings within each source/reference-chromosome lane, then restores
  the original source identifiers. The installed container stays unchanged.
- Apply `pangenome_regional_chunk_bp` (default 25000) with regional softclipping.
  Keep all original assembly sequences unchanged.
- Record every source interval assigned to a chromosome, assigned ambiguously or
  omitted from partitions. Do not invent an assignment for softclipped ends or
  claim they enter an "other" bin. All excluded sequence remains in the assemblies.
- Compare assigned intervals with actual FULL and CLIP source-coordinate spans;
  fail on unknown biological identities, out-of-bounds or overlapping source spans,
  and duplicate graph segment IDs. Report representation discrepancies explicitly.
  Count each graph segment once, separately from haplotype source spans.

These changes intentionally invalidate the Cactus construction cache. Unchanged
assembly, chimera, QC and upstream work uses ordinary Nextflow `-resume` behavior.
There is no cache recovery system. The first run rebuilds the graph; subsequent
unchanged runs should reuse the corrected tool, graph and audits.

## Crest run

Sync the updated repository, then launch from the project directory containing
`gcl_genome_assembly`, `data`, `work` and `comparisons`:

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_chimera_pipeline.sbatch \
  data/assembly_samplesheet.csv data/hic_readsets.csv \
  comparisons/core-graph-production-01 \
  --chimera_break "$PWD/data/chimera_reviews/Sde-CTlk_104.v8.reviewed.tsv" \
  --chimera_hifi_bam_manifest "$PWD/comparisons/chimera-observable-contact-v7/results/assembly/chimeras/mapping_manifest.json" \
  --run_pangenome true --pangenome_validate_only false \
  --pangenome_regional_chunk_bp 25000
```

The BAM manifest is an existing test-run reuse optimization, not a required
production pipeline input. All assignment/filter/gate/chunk replay flags default
to false and should remain off for this acceptance run.

## Review artifacts and acceptance

Return `comparisons/core-graph-production-01.review.tar.gz`.
Under `results/pangenome/373251/`, examine:

- `sequence_policy/sequence_accounting.md`: assigned, ambiguous and omitted totals;
  `split_*.excluded.tsv.gz` records exact original-source exclusion intervals.
- `representation/sequence_graph_audit.md` and `source_stage_accounting.tsv`:
  every source through preprocessing, FULL and CLIP. Exact differences are in
  `source_stage_discrepancies.tsv.gz`.
- `representation/distinct_graph_nodes.tsv`: distinct graph segment/base totals,
  separately from biological source-coordinate spans.
- Existing source, clipping, graph-role and sharing reports: chromosome-arm
  recovery must agree with the supported-chunk comparison. A source path name
  being present is not sufficient evidence of complete representation.

Acceptance requires successful real construction, recovery of the previously
missing arms, no duplicate source assignment, and an explanation of substantial
new representation discrepancies. Source-coordinate span coverage does not prove
that every base inside a walk survives alignment or clipping; retain the existing
clipping/path-length comparisons. Category 5 stays open until this review passes.

Local validation uses the real corrected compiler/CLI and production adapters,
with external Cactus tool doubles for workflow, publication, failure and cache
checks. It does not substitute for the real Crest graph build.

Completed local checks: 138 pangenome unit/regression tests; all 11 input,
reference, graph-role and failure integration cases; and a two-launch production
workflow check confirming cached reuse of the corrected tool, graph construction
and source audits. Both fixture species publish the new accounting outputs.
