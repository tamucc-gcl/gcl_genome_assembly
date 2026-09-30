# Targeted chimera exploration, outside the production workflow

Purpose: distinguish a technical adjacency from a genuine rearrangement and locate
a defensible breakpoint. Neither the global Hi-C minimum nor absence of an N-gap
settles those questions. No changes to production chimera policy in this checkpoint.

## Run on the cluster

Sync the new scripts and test. Keep the current results and work directories; do
not move them or launch the assembly pipeline for this check. From the project root:

```bash
export PYTHONDONTWRITEBYTECODE=1
python3 -m unittest discover -s gcl_genome_assembly/tests -p 'test_*.py'
```

The suite now contains 65 tests, including an end-to-end scratch-cleanup fixture.
Jobs 1501347/1501348 exposed the invalid assumption that staged inputs and temporary
pairs survive scratch execution. No tests or diagnostic analysis were run locally.
After tests pass, submit this read-only investigation of the 1501091 work files:

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_chimera_region.sbatch \
  --joins-work work/18/6be99359df0e14b22942c5ad2c0c2a \
  --evidence-work work/e0/5fb28ccc2dc140819584098058ef01 \
  --assessment-fasta genome_assembly/assembly/scaffold/yahs_round2/Sde-CTlk_104_hap2_round2_scaffolds.fa \
  --last-agp genome_assembly/assembly/scaffold/yahs_round2/Sde-CTlk_104_hap2_round2_scaffolds_final.agp \
  --source-pairs genome_assembly/bam/hic/scaffold/filtered/Sde-CTlk_104_hap2.pairs.gz \
  --results genome_assembly \
  --work-root work \
  --log logs/nextflow_assembly_1501091.log \
  --assembly Sde-CTlk_104_hap2 --scaffold scaffold_3 \
  --start 30000000 --end 36000000 \
  --comparison-assembly Sde-CTlk_104_hap1 \
  --comparison-assembly Sde-CMat_203_hap2
```

These task paths come from the completed run's log. They are invocation arguments,
not sample-specific rules in the script. Requires Python 3.9+ and retained work
files; no new aligner or plotting environment is required. The Slurm request is
one CPU, 8 GB and two hours; this is a bounded initial allocation, not a benchmark.

The command above uses retained, declared outputs. It checks the FASTA SHA256,
validates AGP lengths and N-gaps against that FASTA, and rebuilds only the selected
scaffold's projected pairs using chimera_hic_pairs.py. Source pairs must be mapped
to the last-round AGP input: scaffold/filtered for this round-2 checkpoint.
No read mapping, assembly or scaffolding is rerun. Large regenerated .pairs files
stay in comparisons/ and are excluded from the transfer archive; projection audits
and regional summaries are included. Inspector task wrappers are copied as text
to recover input-path provenance when staged inputs have been cleaned up.

Without an explicit FASTA, the resolver checks files in the task roots and one
directory level below, accepting only content matching the call-table SHA256.
It supports flat, renamed and nested staging. An optional --assessment-fasta path
also requires that exact checksum; never substitute the finished assembly merely
because its scaffold name matches. Job 1501347 stopped before regional analysis.
Its partial comparison directory can remain: each submission creates a new path.

Upload comparisons/chimera-region-<timestamp>-<jobid>.tar.gz and its
logs/chimera-region-<jobid>.out. If it fails, send that log; do not rerun assembly.

## Outputs and interpretation

- provenance.json verifies the exact pre-finishing FASTA checksum against both
  detection and evidence calls. Coordinates are not taken from the finished FASTA.
- regional_alignments.tsv retains PAF records, mapping quality and tags across
  30–36 Mb, without hiding ambiguous alignments. PAF coordinates are 0-based,
  half-open; AGP coordinates retain their original 1-based convention.
- target_context.paf includes all target-haplotype mappings to reference sequences
  encountered in that interval. Companion PAFs show the sister haplotype and an
  independent assembly against those same reference sequences. These allow us to
  compare overlapping versus complementary chr7 coverage and chr12 continuity.
- region.agp, sequence_gaps.tsv and regional_older_joins.tsv show last-round
  components, actual N-gaps and exact recovered older joins. Source history is
  contextual; it must not be arithmetically projected across sequence corrections.
- local_hic.tsv scans the selected region in 100 kb bins using full 2 Mb flanks.
  Cross-boundary and within-flank contacts are compared within three bin-separation
  bands (100–500 kb, 500 kb–1 Mb, 1–2 Mb). Counts and bin-pair opportunities are both
  exported. The statistic uses unbalanced counts and binned separation, not exact
  read distance. It is exploratory, not a significance test or automatic threshold.
- hic_bin_coverage.tsv and hic_matrix.tsv help distinguish loss of cross-boundary
  contacts from low-mappability/low-coverage stripes. Distance matching alone does
  not remove coverage, restriction-site, compartment or mapping biases.
- inspector_context/ contains existing stage-labelled reports and task commands
  as text only. inspector_work_inputs.tsv and inspector_inventory.tsv locate the
  BAMs, original mapping references and reads; BAMs and raw reads are not packaged.

## HiFi follow-up

This first job does not pretend the pre-correction Inspector BAM is mapped to the
current assessment sequence. Review its exact reference and the region's provenance
first. If the sequence correspondence is demonstrably unchanged, inspect existing
read alignments at the matched locus. Otherwise remap HiFi reads against the exact
assessment assembly with competing homologous sequence retained, then inspect
spanning reads, clipping/supplementary alignments, coverage and mapping ambiguity.
Do not map only to an isolated short suspect interval and interpret forced mappings
as support. A lack of spanning reads is informative only with sufficient read length,
coverage and unique flanking sequence.

No-read-support, technical-error and biological-variant conclusions are not emitted
by this collector. The next decision depends on the returned evidence. A positive
technical diagnosis may require a within-contig correction rather than an AGP-gap
cut; that is a separate policy from the current conservative gap-only cutter.
