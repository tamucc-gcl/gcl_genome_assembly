# Bounded rearrangement method comparison

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

Check each expected event against call coordinates/orientation and inspect false positives in controls. Record partial, missed, ambiguous and ineligible outcomes explicitly. Do not equate VCF/BND row counts with event counts: a translocation may be represented by several breakends or separate loss/gain calls. Review duplications for source and destination recovery, not just inserted length. Runtime includes alignment plus calling; these tiny jobs do not predict whole-genome resource requirements.

No automatic accuracy score is provided before reviewing the actual caller formats. A successful exit only means the tools ran. After inspecting outputs, decide which candidate warrants the repeat/real-locus gate and whether any useful legacy analysis remains uncovered.
