# Resolve the chromosome transition and repeat anchoring

## Purpose

Use the existing HiFi mapping, not another mapping run, to test whether reads crossing the candidate interval are anchored outside telomeric repeats. Compare the entire interval with three complete assemblies, including the original harmonization reference. This diagnostic makes no cuts and no automatic biological classification. Scripts contain no sample-specific decisions.

## Run on the cluster after syncing

From the project root, add tidk to the environment already used for the mapping. Tool versions are captured in each output directory.

```bash
module load miniconda3
conda install -y -p "$PWD/comparisons/envs/chimera-hifi" \
  --override-channels -c conda-forge -c bioconda tidk

sbatch gcl_genome_assembly/scripts/comparisons/run_chimera_sequence_context.sbatch \
  --tool-env "$PWD/comparisons/envs/chimera-hifi" \
  --mapping-run comparisons/chimera-hifi-20260930-082608-1501351 \
  --region scaffold_3:33500001-34000000 \
  --motif CCCTAA \
  --comparison sister genome_assembly/assembly/scaffold/yahs_round2/Sde-CTlk_104_hap1_round2_scaffolds.fa \
  --comparison independent genome_assembly/assembly/scaffold/yahs_round2/Sde-CMat_203_hap2_round2_scaffolds.fa \
  --comparison reference genome_assembly/assembly/scaffold/yahs_round2/Sde-CMat_203_hap1_round2_scaffolds.fa
```

All comparison files must be the pre-finishing round-two assemblies. Missing inputs fail before work starts. The job requests 8 CPUs, 32 GB and four hours as initial estimates. Existing BAM/FASTA stay unchanged; output is written under comparisons. No Nextflow or assembly rerun is required. This script was statically reviewed, not executed locally; cluster execution remains the runtime check.

Return the resulting archive and log:

```bash
ls -ltrh comparisons/chimera-sequence-*.tar.gz logs/chimera-sequence-*.out
```

Download the new `comparisons/chimera-sequence-<timestamp>-<jobid>.tar.gz` and matching `logs/chimera-sequence-<jobid>.out`. On failure, send the log and retain the directory. The previous complete HiFi BAM remains reusable.

## Evidence collected

- tidk profiles at 1-kb resolution across the complete target scaffold and comparison assemblies, including terminal controls.
- A 100-bp tidk profile across the 500-kb interval and across the full sequences of reads having any alignment there. Counts are window-based, not a contiguous-repeat-tract call. Both motif orientations must be interpreted with the read alignment strand.
- Full primary read sequences, including reads whose primary placement lies elsewhere, plus all reported BAM placements for selected read names. This supports checking sequence outside motifs rather than accepting MAPQ alone.
- Base-level assembly-to-assembly PAF records with CIGAR and cs tags against complete comparison references, retaining secondary matches. Index lengths are capped to avoid multi-part mapping qualities. The asm5 choice follows the current within-species comparison and is not a universal divergent-species default.
- FASTA index lengths, checksums, paths and versions for provenance and determining whether homologous segments actually approach scaffold ends.

Coordinates: SAM positions are 1-based; PAF coordinates are 0-based half-open. PAF query positions and regional tidk positions refer to the extracted region, not the original scaffold. Add 33,500,000 to zero-based regional positions for this invocation. The complete-scaffold tidk profile needs no region offset. tidk window-position conventions must be checked against its output/version before combining tables. Read motif coordinates refer to original read orientation, not BAM reference orientation.

## Adjudication

1. Locate motif-rich sequence and determine whether the proposed read anchors are outside it. Absence of the selected motif does not demonstrate sequence uniqueness: subtelomeric and other repeats also matter.
2. Review overlapping reads throughout the full chromosome-transition interval, not just its midpoint. Require distinct molecules and inspect substantial anchors on both sides of any localized candidate; record alternative placements, clipping, indels and coverage.
3. Compare corresponding sequences and chromosome-end positions in sister and independent assemblies. Do not assume scaffold labels or membership imply equivalent coordinates.
4. Separate outcomes: supported continuity; localized unsupported connection; repeat/coverage-limited uncertainty; or chromosome-assignment ambiguity. Continuous support alone cannot distinguish a genuine rearrangement from a haplotype switch elsewhere. No cut follows solely from a telomere peak, missing chromosome label, or local Hi-C minimum.

## Pipeline integration after reviewing this result

Keep mapping and evidence collection as reusable modules with assembly checksum, input stage, and software versions. Candidate-level summaries should report the interval, coordinate validity, telomere context, primary coverage and read anchoring evidence, with concise reasons for REVIEW. Reuse an existing BAM only when its reference is verified compatible. Retain raw tidk windows as declared outputs so scratch cleanup cannot remove needed evidence.

Before enabling any new cutting rule, validate against both supported-continuity and known unsupported-join examples. Preserve the existing gap-join route; add a separately validated contig-internal route only if this experiment supports it. Unresolved evidence must not block assembly output or force a split. No sample-specific exception is allowed.

Once this locus is either resolved or explicitly classified as evidence-limited, implement the useful reporting components with targeted tests and then return to the remaining pipeline revamp. Do not prolong the assembly rewrite until every biological rearrangement can be proven.
