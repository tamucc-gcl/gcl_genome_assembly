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
  python=3.11 syri=1.7.1 minimap2=2.28 samtools=1.21
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

## Review gate

Check each expected event against call coordinates/orientation and inspect false positives in controls. Record partial, missed, ambiguous and ineligible outcomes explicitly. Do not equate VCF/BND row counts with event counts: a translocation may be represented by several breakends or separate loss/gain calls. Review duplications for source and destination recovery, not just inserted length. Runtime includes alignment plus calling; these tiny jobs do not predict whole-genome resource requirements.

No automatic accuracy score is provided before reviewing the actual caller formats. A successful exit only means the tools ran. After inspecting outputs, decide which candidate warrants the repeat/real-locus gate and whether any useful legacy analysis remains uncovered.
