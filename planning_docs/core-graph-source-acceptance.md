# Core graph source representation acceptance

The core-graph-02 archive establishes that the CTlk hap1 source name
`chr4_1+chr14_1` is absent already in FULL, while it is present in the finalized
name map. The earlier FULL-to-CLIP name-matching errors are corrected. The archive
does not contain the upstream GAF/PAF or input FASTA dictionaries needed to establish
where this source path was excluded. Do not call the core graph accepted yet.

Local comparison using all uploaded final name maps (sequence counts and total lengths
reconciled against the input identity ledger) covers 6,180 source sequences. Exactly
one chromosome-scale source name is absent from FULL/CLIP: the 74,760,027 bp CTlk hap1
chr4/chr14 composite. There are also 1,516 other source names absent from both stages,
retained in the detailed accounting rather than interpreted as biological deletions.
This local preview derives its dictionaries from name maps; the Crest process uses the
actual finalized FAI files and existing mapping records.

## Implemented next check

`PANGENOME_SOURCE_AUDIT` runs inside the normal PANGENOME workflow, after the compact
FULL/CLIP summary. It:

- Checks every source FAI dictionary against the already audited FASTA record count
  and total length, with explicit species/assembly identities.
- Accounts for every source sequence, including those absent from both graph stages.
- Uses final eligibility chromosome scope, with no sample IDs or genomic coordinates
  built into detection rules.
- Streams existing Cactus GAF/PAF once each; measures unioned query coverage for missing
  source paths and preserves bounded mapping examples. It does not remap reads or genomes.
- Distinguishes absent files from measured zero; excludes synthetic mapping queries;
  treats MAPQ 255 as unavailable. Reference query absence is not a mapping failure.
- Publishes a Markdown source report, all-sequence table, missing-path table, mapping
  examples and machine-readable audit in the existing representation directory.

A chromosome-scale source name missing from CLIP yields **REVIEW**, not biological
absence or automatic success. Non-chromosome source omissions remain in the tables.
Upstream mapped queries missing from FULL localize the discrepancy after mapping,
but do not by themselves identify the exact split/alignment/filter decision.

The Cactus 3.2.1 configuration uses chromosome-assignment coverage and uniqueness
filters. Ambiguous contigs are excluded from the chromosome alignment list. Thus a
composite mapping to two chromosomes is a plausible partitioning issue; that is a
hypothesis, not an established cause for this particular source sequence.

Sources: [pinned split implementation](https://github.com/ComparativeGenomicsToolkit/cactus/blob/v3.2.1/src/cactus/refmap/cactus_graphmap_split.py),
[pinned configuration](https://github.com/ComparativeGenomicsToolkit/cactus/blob/v3.2.1/src/cactus/cactus_progressive_config.xml).

## Cached run

From the Crest project directory, after syncing the repository:

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_chimera_pipeline.sbatch \
  data/assembly_samplesheet.csv data/hic_readsets.csv \
  comparisons/core-graph-03 \
  --chimera_hifi_bam_manifest "$PWD/comparisons/chimera-observable-contact-v7/results/assembly/chimeras/mapping_manifest.json" \
  --chimera_break "$PWD/data/chimera_reviews/Sde-CTlk_104.v7.reviewed.tsv" \
  --run_pangenome true --pangenome_validate_only false
```

This uses ordinary Nextflow resume, the existing work directory and reviewed cuts.
The optional HiFi manifest reuses the prior investigation's mapping; it is not a
production requirement. Do not run another controller against the same cache concurrently.

Review `results/pangenome/373251/representation/source_report.md` and
`source_audit.json` in the selected output directory; the small outputs are included
in `comparisons/core-graph-03.review.tar.gz` by the existing launcher.

## Completion gate

1. Confirm input dictionaries and identities reconcile; classify every missing
   chromosome-scale path using the actual mapping evidence.
2. For mapped source paths absent from FULL, establish the split/alignment cause before
   changing construction settings. If needed, use a bounded replay of chromosome
   assignment; do not blindly relax filters or launch a complete graph rebuild.
3. Correct an established input/naming/partitioning defect; change construction only
   when its evidence warrants doing so. Rebuild only the affected graph if required.
4. Reconcile source representation, FULL/CLIP accounting, sharing and plots across the
   whole cohort before checking off core graph acceptance.

Existing construction cleanup removed `chrom-subproblems/minigraph.split.log` and
`contig_sizes.tsv`. When construction next changes, preserve those compact diagnostics
before deleting scratch, along with the effective Cactus configuration. That cannot
recover the old log. The expensive CACTUS_PANGENOME command is deliberately unchanged
in this audit-only patch so the current graph's cache remains valid.

## Local verification

- 49 focused pangenome tests passed, including 12 source-audit tests covering overlap,
  missing versus zero evidence, GAF-to-PAF differences, synthetic queries, reference
  mapping behavior, invalid coordinates and dictionary reconciliation.
- All 10 input/reference/graph/organelle Nextflow contract checks passed.
- Two-species integration exercises actual GAF/PAF staging with tiny tool doubles;
  ordinary Nextflow resume caches both Cactus construction and the source audit.
  These orchestration checks do not validate the biological graph.
