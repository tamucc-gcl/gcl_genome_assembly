# Bounded rearrangement method comparison

## Current investigation: anchor uniqueness and resource cost

Corrected jobs 1503144/1503145 completed all five eligible cases in each suite. Event labels and reported intervals match the previous maxmatch results reviewed: clean controls, INV, CPG/TDM, INVDP and TRANS, including the same extra 3 bp deletion in the translocation case. Nucmer took 0.10–0.18 seconds and approximately 19–25 MiB peak RSS on these small inputs; this is not evidence of chromosome-scale resource improvement. No observed loss on these fixtures justifies advancing to the real pilot, not general repeat sensitivity claims.

After syncing, rerun the same chromosome pair with the new sixth argument (no environment changes):

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_sv_real_region.sbatch \
  "$PWD/comparisons/environments/syri-1.7.1" \
  genome_assembly/assembly/final/Sde-CMat_203_hap1.fasta \
  genome_assembly/assembly/final/Sde-CBau_104_hap1.fasta \
  chr15_1 chr15_1 mumreference
```

Keep the eight-hour/32 GB pilot request initially. While running, inspect `sstat -j JOBID.batch --format=JobID,AveCPU,MaxRSS -P`; after completion collect `sacct -j JOBID --format=JobID,State,Elapsed,AllocCPUS,TotalCPU,MaxRSS,ExitCode -P`. These are scheduler snapshots, not proof of specific algorithmic behavior. Return the archive, Slurm log and accounting. Retain the prior failed comparison; no production rerun is involved. The real runner preserves maxmatch as the default for old commands and records the selected mode explicitly.

Correction after jobs 1503142/1503143: all eligible cases failed at command parsing because MUMmer 4.0.1 does not accept --mumreference. No alignment/caller sensitivity or resource conclusions can be drawn from those runs. In the pinned version, reference-unique anchors are the default when neither --maxmatch nor --mum is specified (verified in the v4.0.1 src/umd/nucmer_cmdline.yaggo source). The runner now maps its user-facing mumreference label to an empty anchor-option array, while maxmatch still supplies --maxmatch. The submission commands below remain unchanged. Repeat both small jobs after syncing; no environment changes needed.

Real pilot 1502945 timed out during nucmer after eight hours, with an empty delta and no SyRI results. Slurm reported 121096712K MaxRSS (approximately 115.5 GiB), despite the 32 GB request. Its near-zero TotalCPU cannot establish idle behavior: signal-interrupted accounting may omit child CPU. The original zero exit_status was erroneous; the real runner now requires explicit completion and handles termination signals. No evidence of an OOM kill was supplied.

The user reports approximately 50% repeat content from previous annotations. This makes excessive nonunique anchors a plausible explanation, not a demonstrated cause. The next controlled comparison changes only nucmer's anchor mode to reference-unique (the MUMmer 4 default). This requires uniqueness in the reference, not in both assemblies, so multiple query copies can still be anchored. Events confined to repeated reference sequence may lose support. Do not promote this setting based on speed alone.

After syncing, use the existing environment:

```bash
for suite in basic repeats; do
  sbatch gcl_genome_assembly/scripts/comparisons/run_sv_benchmark.sbatch \
    syri-mummer "$PWD/comparisons/environments/syri-1.7.1" "$suite" mumreference
done
```

Return both printed archives and their Slurm logs. Existing maxmatch outputs provide the baseline; no repeat of those or of the real chromosome run is needed yet. Fixture hashes should match the corresponding earlier suite. The fourth argument is optional and defaults to maxmatch, preserving prior commands. Output names/settings record the anchor mode, and separate GNU time reports capture alignment and SyRI resources. These tiny fixtures establish event retention, not performance at 50% repeat content. Review lost/retained raw blocks, filtering, event labels, source/destination recovery and false calls before choosing a monitored real-data rerun.

## Current next step: repeat/boundary gate

After syncing, use the existing successful MUMmer/SyRI environment; no reinstall is needed:

```bash
export PYTHONDONTWRITEBYTECODE=1
python3 -m unittest discover -s gcl_genome_assembly/tests -p 'test_*.py'
# Submit after tests pass.
sbatch gcl_genome_assembly/scripts/comparisons/run_sv_benchmark.sbatch \
  syri-mummer "$PWD/comparisons/environments/syri-1.7.1" repeats
```

The original basic fixtures remain unchanged. The repeats suite plants identical 2 kb segments on both chromosomes at five positions near event boundaries, then applies the same four engineered changes. Includes an unchanged repeat-bearing control. Some repeat origins and breakpoint placements are inherently ambiguous: inspect event recovery, source/destination, orientation, false positives and alternative equivalent boundaries rather than demanding a single exact breakpoint. This is still a small synthetic stress test, not a realistic repeat landscape or validation of nested graph alleles. Return the printed archive and Slurm log, including failures.

## Prepared next gate: real homologous chromosome pair

Repeat job 1502943 completed all five eligible cases. The repeat-bearing unchanged control had no calls; inversion and inverted duplication had exact engineered intervals. Tandem duplication changed from CPG to TDM, with an expanded repeat-associated interval (reference 22,000 bp, query 42,002 bp). The translocation linked the correct chromosomes but extended to 22,001 bp around the 20,000 bp engineered move; the same extra 3 bp source-junction deletion persisted. Thus event recovery is encouraging, but repeat boundaries and local false calls remain unresolved. Report native CPG/TDM classes and interval spans without treating them as exact net gain, or blindly combining overlapping blocks. Proceed to a bounded real chromosome pair for practical assessment, not production acceptance.

After reviewing the repeat gate, select one complete homologous chromosome from two chromosome-scale assemblies of the same species. Use the graph reference assembly for the reference side and a different individual's eligible assembly for the query side. Confirm homolog identity from harmonization/alignment evidence; do not infer it merely from matching labels. Selection is explicit, with no hard-coded samples. For a first pilot, prefer a smaller complete chromosome with unambiguous correspondence; keep its full sequence rather than cropping unequal coordinate windows.

The runner accepts five positional arguments: existing environment, reference FASTA, query FASTA, reference sequence ID, query sequence ID. Template (replace all uppercase placeholders):

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_sv_real_region.sbatch \
  "$PWD/comparisons/environments/syri-1.7.1" \
  REFERENCE_FASTA QUERY_FASTA REFERENCE_CHROMOSOME_ID QUERY_CHROMOSOME_ID
```

Existing adjacent .fai indexes are required; production FASTAs/indexes are not changed. Local copies of the selected chromosomes are renamed chr1 and the original IDs recorded. A 100 Mb cap per sequence bounds this pilot's workload; it is not a biological chromosome eligibility threshold. The 32 GB/8 h request is provisional. Alignment, filtering and calling record separate elapsed time and maximum resident memory through GNU time. Retain the full local comparison folder; return its printed archive and logs/sv-real-JOBID.out. FASTAs and delta files are excluded from the return archive; coordinate tables, calls, resource reports, selected-sequence hashes and environment build records are included. A scheduler kill may prevent archive creation; return the Slurm log in that case.

A single chromosome pair cannot assess interchromosomal translocations and has no independent event truth. Treat unmatched sequence as unresolved in this restricted search space, not absence from the other whole genome. Inspect alignment coverage, repeat-associated calls, boundary artifacts and resource requirements. This gate evaluates practical behavior, not precision/recall or biological confirmation. A later multi-chromosome comparison is required before accepting translocation reporting on real data.

Scripts and tests were authored and statically reviewed locally; execution remains on the cluster. Neither gate changes production workflows.

This comparison does not change the production graph, VCF or assembly processes. Keep results and work; no baseline move or Nextflow rerun is required.

The generator creates two 300 kb synthetic chromosomes with a fixed seed and records FASTA checksums. Each independent query has one known inversion, tandem duplication, dispersed inverted duplication, or interchromosomal cut-and-paste event. Unchanged and fragmented unchanged assemblies are negative controls. Truth coordinates are zero-based half-open; compare caller coordinates using their documented conventions, not literal equality across formats.

This is a basic capability test, not a realistic sensitivity estimate. Random sequence deliberately avoids most repeat ambiguity. Later gates must include repeated sequence, nested/nonreference graph alleles, divergence and selected real loci. Neither assembly caller tests GREF allele decomposition. Do not select a production method from these simple fixtures alone.

## Methods

SyRI is the first candidate for chromosome-scale rearrangements. Follow its [official alignment example](https://schneebergerlab.github.io/syri/pipeline.html): minimap2 asm5 with eqx and SAM input. It requires chromosome correspondence; the fragmented control is explicitly not eligible.

SVIM-asm provides a contrasting assembly-to-reference caller using its [documented haploid workflow](https://github.com/eldariont/svim-asm): minimap2 asm5/cs/r2k and sorted indexed BAM. Haploid refers to each synthetic query assembly, not the biology of the user's samples. The authors say it is no longer actively maintained; this weighs against adopting it permanently, even if it performs well here.

The environments below pin candidate/tool versions; the runner also captures installed package build metadata. Environment solves and runtime compatibility remain untested locally. Separate environments avoid coupling their dependencies. SyRI 1.7.1 is a deliberate baseline with Python 3.11 support documented in its release notes, not a claim that it is the newest release.

## Cluster commands

From the project root, after syncing, run tests first:

```bash
export PYTHONDONTWRITEBYTECODE=1
python3 -m unittest discover -s gcl_genome_assembly/tests -p 'test_*.py'
```

Create these comparison-only environments once, on a cluster host where environment installation is permitted:

```bash
module load miniconda3
mkdir -p comparisons/environments
conda create -y --strict-channel-priority -c conda-forge -c bioconda \
  -p "$PWD/comparisons/environments/syri-1.7.1" \
  python=3.11 syri=1.7.1 pandas=2.2.3 minimap2=2.28 samtools=1.21
conda create -y --strict-channel-priority -c conda-forge -c bioconda \
  -p "$PWD/comparisons/environments/svim-asm-1.0.3" \
  python=3.10 svim-asm=1.0.3 minimap2=2.28 samtools=1.21
```

If a solve fails, retain its message; do not silently substitute versions. Submit each job after its environment exists:

```bash
sbatch gcl_genome_assembly/scripts/comparisons/run_sv_benchmark.sbatch \
  syri "$PWD/comparisons/environments/syri-1.7.1"
sbatch gcl_genome_assembly/scripts/comparisons/run_sv_benchmark.sbatch \
  svim-asm "$PWD/comparisons/environments/svim-asm-1.0.3"
sbatch gcl_genome_assembly/scripts/comparisons/run_gref_inspection.sbatch \
  --pangenome-dir genome_assembly/pangenome/373251 --taxid 373251 \
  --conda-cache /work/birdlab/.conda_builds
```

Each job prints its archive path. Return both comparisons/sv-benchmark-*.tar.gz archives, the new comparisons/gref-inspection-*.tar.gz and matching Slurm logs. Alignments and synthetic FASTAs stay in their comparison folders; caller outputs, truth, checksums, commands, package builds and per-case elapsed seconds are archived. A failed caller does not prevent collection of other cases; the job still exits unsuccessfully if any eligible case failed. Setup failures appear in the Slurm log.

## First results: jobs 1502925–1502927

GREF: all 4,808 standard and 5,544 raw window records resolved their expected slots without called placeholders or unresolved assignments. The reference individual's partner column has one missing placeholder in every record. Its biological slot was called in 4,642 standard records and missing in 166; the separate reference slot was called in 4,688 and missing in 120. Raw corresponding counts were 5,098/446 and 5,270/274. No AC/AN, GT-range or AT-cardinality flags were raised. These selected windows do not establish genomewide missingness or phase accuracy.

SVIM-asm completed all six cases in 2–3 seconds each including alignment. Neither unchanged control produced a call. The inversion had the correct 40 kb interval but FILTER=incomplete_inversion. Tandem and inverted duplications were represented as 20 kb insertions, without duplication/source labels. The translocation was a 20 kb deletion plus 20 kb insertion, without a linking translocation/BND call. These are partial representations, not four fully classified events. Small boundary shifts and the tandem copy's placement require allele equivalence checks before judging coordinate accuracy.

SyRI failed all five eligible cases, including the unchanged control, with `ValueError: buffer source array is read-only`. The captured environment installed pandas 3.0.6 and NumPy 2.4.6. The accompanying chained-assignment warning and [pandas Copy-on-Write documentation](https://pandas.pydata.org/docs/user_guide/copy_on_write.html) support a pandas 3 compatibility diagnosis; this is not evidence of rearrangement sensitivity failure. The initial environment specification left pandas unconstrained. It now pins pandas 2.2.3; this repair still needs cluster confirmation.

Repair the existing comparison environment and rerun only SyRI (no GREF/SVIM-asm/full-pipeline rerun needed):

```bash
module load miniconda3
conda install -y --strict-channel-priority -c conda-forge -c bioconda \
  -p "$PWD/comparisons/environments/syri-1.7.1" \
  python=3.11 syri=1.7.1 pandas=2.2.3 minimap2=2.28 samtools=1.21
sbatch gcl_genome_assembly/scripts/comparisons/run_sv_benchmark.sbatch \
  syri "$PWD/comparisons/environments/syri-1.7.1"
```

Submit only after the environment update succeeds. Return the new comparisons/sv-benchmark-syri-*.tar.gz and logs/sv-benchmark-JOBID.out. Keep the failed-run archive as environment provenance. No production method is selected yet.

## Interpretation criteria

## MUMmer comparison after SyRI job 1502932

The pandas pin restored successful execution for all five eligible cases. SyRI correctly reported the inversion (with 206 bp adjacent NOTAL sequence) and no control variants. Both duplications remained insertions. The translocation produced only a source deletion, with neither a destination insertion nor a linked rearrangement. Successful completion is therefore not sufficient for acceptance.

Use the same generator and compare manifest FASTA hashes against job 1502932. The new `syri-mummer` mode follows the official SyRI MUMmer example: nucmer --maxmatch -c 100 -b 500 -l 50, delta-filter -m -i 90 -l 100, show-coords -THrd, then SyRI with coordinates and delta. It also repeats only the tiny minimap2 alignment for diagnostic comparison, not its SyRI call. All SAM alignment flags and CIGAR strings are retained in minimap2_alignments.tsv; sequence/quality strings are omitted. Raw and filtered MUMmer coordinates and delta files are included in the return archive. The original modes remain available.

After syncing, add MUMmer to the working comparison environment, retaining the core tested versions:

```bash
module load miniconda3
conda install -y --strict-channel-priority -c conda-forge -c bioconda \
  -p "$PWD/comparisons/environments/syri-1.7.1" \
  python=3.11 syri=1.7.1 pandas=2.2.3 numpy=2.4.6 \
  minimap2=2.28 samtools=1.21 mummer4=4.0.1
# Submit only after the installation succeeds.
sbatch gcl_genome_assembly/scripts/comparisons/run_sv_benchmark.sbatch \
  syri-mummer "$PWD/comparisons/environments/syri-1.7.1"
```

Return the printed comparisons/sv-benchmark-syri-mummer-*.tar.gz archive and logs/sv-benchmark-JOBID.out, including on caller failure. No production rerun, baseline move, or GREF/SVIM rerun is required. The comparison runner has been statically reviewed, not executed locally. This mode's elapsed time includes the diagnostic minimap2 run as well as MUMmer and SyRI and is not directly comparable to the earlier timing.

Review whether each true source-to-destination block exists in raw MUMmer coordinates, survives filtering, and appears in SyRI's rearrangement calls. Compare minimap2 CIGAR insertions/deletions and secondary/supplementary alignments to locate where information was lost. Distinguish missing alignment, filtered alignment and unclassified retained alignment. Do not infer method superiority from a successful exit or a change in raw call count.

## Final review gate

### MUMmer result: job 1502940

All five eligible cases completed; the fragmented control remained explicitly ineligible. The unchanged control had no variants. The 40 kb inversion was recovered exactly, without the previous 206 bp NOTAL flank. Tandem duplication was represented as CPG (copy gain), inverted duplication as INVDP/copygain with reversed alignment orientation, and the interchromosomal move as TRANS with the correct source and destination chromosomes.

Raw and filtered MUMmer coordinate rows were identical in the three duplication/translocation cases inspected. They preserve the overlapping, reversed and cross-chromosomal blocks needed for classification. The diagnostic minimap2 records instead encode these changes inside long alignments as 20 kb I/D CIGAR operations, with no separate source-to-destination block. This supports alignment representation as the principal cause of the earlier classification differences for these fixtures; it does not prove that all minimap2 settings would behave this way.

Boundary accuracy is not perfect. INVDP and TRANS span 20,001 bp rather than the engineered 20,000 bp because alignments extend into adjacent matching bases. The translocation output also contains an extra 3 bp deletion near the source junction, absent from the engineered truth, associated with an overlapping gapped flank alignment. Treat this as an apparent alignment/calling artifact to investigate, not a new true event.

The tandem CPG summary reports 40,002 bp of query interval (120000–160001 inclusive), covering the duplicated region with flanking bases; this is not 40,002 newly gained bases. The engineered gain is 20,000 bp. Preserve event classes, source/destination intervals, orientation, hierarchy and explicit length definitions; do not sum parent and alignment-child rows or convert summary lengths directly into net gains.

Decision: advance MUMmer+SyRI as the leading chromosome-scale rearrangement candidate, not yet a production default. Next gate: repeat/boundary ambiguity fixtures, then a bounded real chromosome comparison with measured alignment and calling resources. Retain explicit ineligibility for fragmented assemblies. Keep GREF variant catalog and assembly rearrangement results distinct. SVIM-asm need not become a parallel default branch on these results. No production changes or local analysis execution were made for this review.

Check each expected event against call coordinates/orientation and inspect false positives in controls. Record partial, missed, ambiguous and ineligible outcomes explicitly. Do not equate VCF/BND row counts with event counts: a translocation may be represented by several breakends or separate loss/gain calls. Review duplications for source and destination recovery, not just inserted length. Runtime includes alignment plus calling; these tiny jobs do not predict whole-genome resource requirements.

No automatic accuracy score is provided before reviewing the actual caller formats. A successful exit only means the tools ran. After inspecting outputs, decide which candidate warrants the repeat/real-locus gate and whether any useful legacy analysis remains uncovered.
