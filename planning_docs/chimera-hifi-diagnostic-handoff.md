# HiFi diagnostic after regional review 1501349

This is a separate comparison job, not a Nextflow change. It maps all of this individual's HiFi reads against the complete, checksum-verified hap2 assessment assembly, then exports regional evidence. Assembly files are never modified. No local execution or tests were performed; the first cluster run remains the runtime check.

## Prepare tools once on the cluster

Run from the project root. This dedicated environment avoids depending on a transient Nextflow environment. If an existing environment has these tools, pass its absolute prefix instead.

```bash
module load miniconda3
mkdir -p comparisons/envs logs
conda create -y -p "$PWD/comparisons/envs/chimera-hifi" \
  --override-channels -c conda-forge -c bioconda \
  minimap2=2.28 samtools=1.21
```

## Submit

The read FASTQ below is the retained input recorded in Inspector's task wrapper. The diagnostic checks that it exists before mapping. If it has been cleaned up, supply the equivalent complete HiFi FASTQ for the individual; do not supply the HiFi BAM directly to minimap2.

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_chimera_hifi.sbatch \
  --tool-env "$PWD/comparisons/envs/chimera-hifi" \
  --fasta genome_assembly/assembly/scaffold/yahs_round2/Sde-CTlk_104_hap2_round2_scaffolds.fa \
  --sha256 2cbed25444ba9e50e49aea7751d6068a22da49c89602413f4c5976ecf21acf23 \
  --reads work/99/b12984f10ed0e298e6cebe24b617a8/Sde-CTlk_104.fastq.gz \
  --region scaffold_3:30000001-36000000
```

The region is 1-based inclusive, equivalent to [30,000,000, 36,000,000) in the earlier tables. Defaults request 16 CPUs, 64 GB and 24 hours; these are starting resource estimates, not measured runtime requirements. The script uses the shared comparison directory rather than ephemeral scratch, so outputs survive job completion. Allow space for a complete sorted BAM and sorting intermediates. Each submission creates a new directory; do not resubmit a successful mapping merely to change the reporting interval.

## Return the results

The last lines of the job log give the exact archive path. List outputs with:

```bash
ls -ltrh comparisons/chimera-hifi-*-reports.tar.gz logs/chimera-hifi-*.out
```

Download the new `comparisons/chimera-hifi-<timestamp>-<jobid>-reports.tar.gz` and matching `logs/chimera-hifi-<jobid>.out`. On failure, send the log and retain the comparison directory rather than rerunning assembly.

The full `hifi.bam` and index stay in that comparison directory. The archive contains mapping provenance, tool versions, mapping log, alignment statistics, all reported placements of regional reads as compressed TSV, and 1-kb primary-alignment depth summaries at MAPQ 0 and 20. Alignment positions are 1-based; depth windows are 0-based, half-open. Secondary placements are limited to minimap2's configured maximum of 20 and its reporting heuristics; they are not an exhaustive enumeration.

## Interpretation and limits

Inspect CIGAR and supplementary alignment tags for clustered clipping, indels and continuity near 33.7–33.8 Mb, and compare with surrounding sequence. Multiple distinct well-anchored reads crossing a localized candidate support continuity. No crossing reads is not decisive without adequate depth, read length and unique flanks. A single read need not span the whole broad chromosome-transition interval. Depth excludes supplementary, secondary, duplicate, QC-failed and unmapped records; supplementary records remain in the alignment table.

The reference contains the complete target haplotype, not both haplotypes. The individual's reads contain both, so alternate-haplotype mismatches or split mappings must not automatically be called errors. If ambiguity persists, comparison against both haplotypes is a follow-up experiment; mapping-quality changes under that competing reference require separate interpretation.

Do not infer an automatic cut from a depth trough, alignment count, or Inspector call alone. This diagnostic emits no biological verdict and makes no pipeline policy changes. Reuse the retained BAM for follow-up rather than repeating mapping unnecessarily.
