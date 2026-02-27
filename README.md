# Genome Assembly and Scaffolding Pipeline

**Author:** Jason Selwyn  
**Framework:** Nextflow DSL2  
**Executor:** SLURM (HPC)

A comprehensive, modular Nextflow pipeline for diploid genome assembly and scaffolding from PacBio HiFi and Hi-C sequencing data. The pipeline takes raw reads through phased contig assembly, iterative scaffolding, gap filling, decontamination, and multi-stage quality assessment, producing chromosome-level haplotype-resolved assemblies with extensive QC reporting.

---

## Table of Contents

1. [Pipeline Overview](#pipeline-overview)
2. [Requirements](#requirements)
3. [Quick Start](#quick-start)
4. [Input](#input)
5. [Pipeline Steps](#pipeline-steps)
6. [Parameters](#parameters)
7. [Output](#output)
8. [Execution Profiles](#execution-profiles)
9. [Resource Configuration](#resource-configuration)
10. [Caching and Resume](#caching-and-resume)

---

## Pipeline Overview

```
HiFi BAM + Hi-C FASTQ
        │
        ├──► BAM → FASTQ conversion
        ├──► Hi-C adapter trimming (fastp)
        │
        ▼
   HIFIASM (phased diploid assembly)
        │
        ├──► [Optional] Purge duplicates (purge_dups)
        ├──► [Optional] Misassembly correction — contigs (Inspector)
        ├──► [Optional] Decontamination — contigs (NCBI FCS-GX)
        │
        ▼
   Hi-C mapping (BWA-MEM2) + filtering (pairtools)
        │
        ▼
   Scaffolding Round 1 (YaHS)
        │
        ├──► [Optional] Misassembly correction — scaffolds (Inspector)
        ├──► [Optional] Decontamination — scaffolds (NCBI FCS-GX)
        ├──► [Optional] Scaffolding Round 2 (YaHS)
        │
        ▼
   Gap Filling (TGSGapCloser)
        │
        ▼
   Finalization & Reporting
        ├── Snail plots (BlobToolKit)
        ├── Contact maps (cooler / HiCExplorer)
        ├── Coverage books
        ├── Pairwise dotplots (minimap2 / pafr)
        ├── Telomere detection
        ├── Compartment & TAD calling
        └── Compiled QC report
```

Assembly QC (QUAST, MERQURY, BUSCO, HiFi mapping statistics) is run at every major checkpoint: initial contigs, purged contigs, corrected contigs, decontaminated contigs, scaffolds round 1, corrected scaffolds, decontaminated scaffolds, scaffolds round 2, and gap-filled assemblies.

---

## Requirements

**Software (managed automatically via Conda / Singularity):**

| Category | Tools |
|---|---|
| Assembly | hifiasm, gfatools |
| Scaffolding | YaHS, BWA-MEM2, samtools |
| Gap filling | TGSGapCloser |
| Duplicate purging | purge_dups, minimap2 |
| Misassembly correction | Inspector |
| Hi-C processing | pairtools, cooler, HiCExplorer |
| Decontamination | NCBI FCS-GX (Singularity), FCS-adaptor |
| Assembly QC | QUAST, MERQURY (meryl), BUSCO |
| Read QC | FastQC, fastp |
| Visualization | BlobToolKit (snail plots), pafr (dotplots), R/ggplot2 |
| Alignment | minimap2 |

**Infrastructure:**

- Nextflow ≥ 22.10
- SLURM scheduler (or local execution)
- Conda and/or Singularity
- Sufficient storage for large databases (FCS-GX ~500 GB)

---

## Quick Start

```bash
# Minimal run
nextflow run main.nf \
  --sample_sheet samples.csv \
  --outdir ./results \
  -profile slurm

# With custom options
nextflow run main.nf \
  --sample_sheet samples.csv \
  --outdir ./results \
  --busco_lineage insecta_odb10 \
  --decon_source_taxid 7041 \
  --hifiasm_primary false \
  --inspector_run_on_contigs true \
  --inspector_run_on_scaffolds true \
  -profile slurm \
  -resume
```

---

## Input

### Sample Sheet

A CSV file with a header row provided via `--sample_sheet`. Rows with any empty column are automatically skipped.

| Column | Description |
|---|---|
| `sample_id` | Unique sample identifier (used throughout all output naming) |
| `hifi_bam` | Path to PacBio HiFi reads in BAM format |
| `hic_r1` | Path to Hi-C forward reads (FASTQ/FASTQ.GZ) |
| `hic_r2` | Path to Hi-C reverse reads (FASTQ/FASTQ.GZ) |

**Example `samples.csv`:**

```csv
sample_id,hifi_bam,hic_r1,hic_r2
Sde-CBau_104,data/raw_bam/m84066_251218_031205_s4.hifi_reads.bc2064.bam,data/raw_fastq/hic/Sde_CBau_104_R1.fq.gz,data/raw_fastq/hic/Sde_CBau_104_R2.fq.gz
Sde_CPla_115,data/raw_bam/m84066_260114_212201_s2.hifi_reads.bc2086.bam,data/raw_fastq/hic/Sde_CPla_115_R1.fq.gz,data/raw_fastq/hic/Sde_CPla_115_R2.fq.gz
Sde_CTlk_101,,data/raw_fastq/hic/Sde_CTlk_101_R1.fq.gz,data/raw_fastq/hic/Sde_CTlk_101_R2.fq.gz
```

In the example above, `Sde_CTlk_101` would be skipped because its `hifi_bam` field is empty.

### Databases

The pipeline can automatically download and cache the following databases. Paths default to `${params.db_base}` (default: `/work/birdlab/databases`).

| Database | Parameter | Description |
|---|---|---|
| BUSCO lineage | `--busco_lineage` | Lineage dataset name (e.g., `actinopterygii_odb10`, `insecta_odb10`) |
| FCS-GX | `--gxdb_dir` | NCBI Foreign Contamination Screen database (~500 GB) |

---

## Pipeline Steps

### Step 0 — Database Setup

Downloads and prepares the FCS-GX decontamination database in parallel with read processing. Only runs if decontamination is enabled.

### Step 1 — BAM to FASTQ Conversion

Converts PacBio HiFi BAM files to FASTQ using `samtools fastq`. The FASTQ is reused by multiple downstream steps (assembly, Inspector, purge_dups, gap filling, coverage analysis).

### Step 2 — Hi-C Read Trimming

Trims Hi-C reads with `fastp` to remove adapters and low-quality bases. Both raw and trimmed reads are quality-checked with FastQC.

### Step 3 — Read Merging

Combines trimmed Hi-C reads with HiFi FASTQ into a unified channel for assembly.

### Step 4 — Phased Assembly (HIFIASM)

Produces a phased diploid assembly from HiFi reads with optional Hi-C integration for phasing. Outputs two haplotype assemblies (hap1, hap2) as FASTA and GFA.

### Step 4a — Purge Duplicates *(optional)*

Removes haplotig duplications from HIFIASM output using `purge_dups`. Runs per-haplotype. Enable with `--run_purge_dups true`. Recommended when HIFIASM's built-in purging (`-l`) is insufficient.

### Step 5 — Misassembly Correction on Contigs *(optional)*

Uses Inspector to map HiFi reads back to the assembly, identify structural and small-scale errors, and break contigs at misassembly sites. Controlled by `--inspector_run_on_contigs`.

### Step 5.5 — Decontamination of Contigs *(optional)*

Screens contigs against the NCBI FCS-GX database and removes contaminant sequences. Optionally includes adapter/vector screening with FCS-adaptor. Controlled by `--decon_run_on_contigs`.

### Step 6–7 — Hi-C Mapping and Filtering

Maps trimmed Hi-C reads to the assembly (per-haplotype) with BWA-MEM2, then filters the BAM with pairtools to remove unmapped, multi-mapped, and duplicate read pairs. Produces filtered BAM, sorted pairs.gz, and parse/dedup statistics.

### Step 8 — Hi-C Scaffolding Round 1 (YaHS)

Scaffolds each haplotype using filtered Hi-C reads with YaHS. Outputs scaffolded FASTA and AGP file.

### Step 8.5 — Misassembly Correction on Scaffolds *(optional)*

Applies Inspector to scaffolded assemblies. Uses higher default thresholds than contig correction to avoid breaking legitimate scaffold joins. Controlled by `--inspector_run_on_scaffolds`.

### Step 9 — Decontamination of Scaffolds *(optional)*

Re-screens scaffolds for contamination. Particularly useful if scaffolding merged clean and contaminated contigs. Controlled by `--decon_run_on_scaffolds`.

### Steps 10–12 — Scaffolding Round 2 *(conditional)*

If scaffold correction or decontamination was performed, the pipeline automatically re-maps Hi-C reads to the corrected/decontaminated scaffolds and runs a second round of YaHS scaffolding. This can be explicitly enabled/disabled with `--run_scaffold_round2`.

### Step 13 — Gap Filling (TGSGapCloser)

Fills gaps in the final scaffolded assemblies using HiFi reads. Operates on whichever scaffold output is most downstream (round 2 > decontaminated > corrected > round 1).

### Step 14 — Finalization and Reporting

Produces final visualizations and reports:

- **Snail plots** — BlobToolKit assembly visualization (BUSCO + assembly stats)
- **Contact maps** — Hi-C contact matrices at multiple resolutions (cooler/HiCExplorer)
- **Coverage books** — HiFi coverage visualization across the assembly
- **Pairwise dotplots** — Hap1 vs hap2 (or all-vs-all) genome alignments via minimap2 + pafr
- **Telomere scanning** — Detection of telomeric repeats at sequence ends
- **Compartment calling** — A/B compartment identification from Hi-C data
- **TAD calling** — Topologically associated domain boundaries
- **Compiled QC report** — Aggregated TSV + HTML with metrics across all assembly stages
- **Summary report** — Markdown + HTML with embedded snail plots, contact maps, and dotplots

### Assembly QC (runs at every checkpoint)

The `ASSEMBLY_QC` subworkflow runs at up to 9 checkpoints throughout the pipeline:

| Checkpoint | Label |
|---|---|
| Initial HIFIASM contigs | `initial` |
| After purge_dups | `purged` |
| After contig correction | `contig_corrected` |
| After contig decontamination | `contig_decontam` |
| After scaffolding round 1 | `scaffold` |
| After scaffold correction | `scaffold_corrected` |
| After scaffold decontamination | `scaffold_decontam` |
| After scaffolding round 2 | `scaffold_round2` |
| After gap filling | `gap_filled` |

Each QC run produces QUAST statistics, MERQURY k-mer QV and completeness, BUSCO gene completeness, and HiFi read mapping statistics.

---

## Parameters

All parameters use a flattened naming convention for easy command-line override. Default values are shown.

### Core Parameters

| Parameter | Default | Description |
|---|---|---|
| `--sample_sheet` | *required* | Path to input CSV sample sheet |
| `--outdir` | `./results` | Output directory |
| `--publish_dir_mode` | `link` | Nextflow publishDir mode (`link`, `copy`, `symlink`) |

### HIFIASM Assembly

| Parameter | Default | Description |
|---|---|---|
| `--hifiasm_useHiC` | `true`* | Use Hi-C reads for phasing |
| `--hifiasm_primary` | `false` | Output primary + alternate instead of hap1 + hap2 |
| `--hifiasm_k` | `51` | K-mer length (must be <64) |
| `--hifiasm_w` | `51` | Minimizer window size |
| `--hifiasm_f` | `37` | Bloom filter bits; 0 to disable |
| `--hifiasm_D` | `5.0` | Drop k-mers occurring >FLOAT×coverage times |
| `--hifiasm_N` | `100` | Max overlaps per oriented read |
| `--hifiasm_r` | `3` | Rounds of error correction |
| `--hifiasm_z` | `0` | Adapter length to remove |
| `--hifiasm_maxKOCC` | `2000` | K-mer occurrence cap for repetitive overlap rescue |
| `--hifiasm_hgSize` | `auto` | Estimated haploid genome size (e.g., `1g`, `500m`); `auto` to infer |
| `--hifiasm_a` | `4` | Rounds of assembly cleaning |
| `--hifiasm_m` | `10000000` | Max bubble size in contig graphs |
| `--hifiasm_p` | `0` | Max bubble size in unitig graphs |
| `--hifiasm_n` | `3` | Remove tip unitigs with ≤N reads |
| `--hifiasm_x` | `0.8` | Max overlap drop ratio |
| `--hifiasm_y` | `0.2` | Min overlap drop ratio |
| `--hifiasm_u` | `1` | Post-join step (1=enable, 0=disable) |
| `--hifiasm_homCov` | `auto` | Homozygous read coverage; `auto` to infer |
| `--hifiasm_lowQ` | `70` | Output regions with ≥N% inconsistency in BED; 0 to disable |
| `--hifiasm_bCov` | `0` | Break contigs at positions with <N-fold coverage; 0 to disable |
| `--hifiasm_hCov` | `-1` | Break contigs at positions with >N-fold coverage; -1 to disable |
| `--hifiasm_mRate` | `0.75` | Works with `--hifiasm_bCov` / `--hifiasm_hCov` |
| `--hifiasm_ctgN` | `3` | Remove tip contigs with ≤N reads |
| `--hifiasm_l` | `3` | Purge level (0–3) |
| `--hifiasm_s` | `0.55` | Similarity threshold for purge |
| `--hifiasm_O` | `1` | Min number of overlaps for discard |
| `--hifiasm_nHaplotypes` | `2` | Number of haplotypes |
| `--hifiasm_dualScaf` | `false` | Enable dual scaffolding |
| `--hifiasm_scafGap` | `300000` | Scaffold gap size |

#### HIFIASM Telomere Parameters

| Parameter | Default | Description |
|---|---|---|
| `--telomere_motif` | `CCCTAA` | Telomeric repeat motif (used by both HIFIASM and telomere scanning) |
| `--hifiasm_teloP` | `1` | Non-telomeric penalty |
| `--hifiasm_teloD` | `2000` | Max drop |
| `--hifiasm_teloS` | `500` | Min score for telomere reads |

### Purge Duplicates

| Parameter | Default | Description |
|---|---|---|
| `--run_purge_dups` | `false` | Enable purge_dups after HIFIASM |

### Misassembly Correction (Inspector)

**Contig-level parameters:**

| Parameter | Default | Description |
|---|---|---|
| `--inspector_run_on_contigs` | `true` | Enable Inspector on contigs |
| `--inspector_contig_skip_baseerror` | `true` | Skip base-level error correction |
| `--inspector_contig_min_depth` | `null` (auto) | Min read depth; default ~20% of average |
| `--inspector_contig_min_contig_length` | `10000` | Min contig length to evaluate |
| `--inspector_contig_min_contig_length_assemblyerror` | `1000000` | Min contig length for structural error detection |
| `--inspector_contig_min_assembly_error_size` | `50` | Min error size (bp) |
| `--inspector_contig_max_assembly_error_size` | `4000000` | Max error size (bp) |

**Scaffold-level parameters:**

| Parameter | Default | Description |
|---|---|---|
| `--inspector_run_on_scaffolds` | `true` | Enable Inspector on scaffolds |
| `--inspector_scaffold_skip_baseerror` | `true` | Skip base-level error correction |
| `--inspector_scaffold_min_depth` | `null` (auto) | Min read depth |
| `--inspector_scaffold_min_contig_length` | `10000` | Min contig length to evaluate |
| `--inspector_scaffold_min_contig_length_assemblyerror` | `10000000` | Min contig length for structural errors (higher for scaffolds) |
| `--inspector_scaffold_min_assembly_error_size` | `50` | Min error size (bp) |
| `--inspector_scaffold_max_assembly_error_size` | `10000000` | Max error size (bp; higher for scaffolds) |

### Hi-C Mapping and QC

| Parameter | Default | Description |
|---|---|---|
| `--bwa_mem2_hic_args` | `""` | Extra arguments for BWA-MEM2 |
| `--hic_coverage_window` | `100000` | Window size for coverage calculation |
| `--hic_min_mapq` | `30` | Min MAPQ for filtered BAM |
| `--hic_resolutions` | `2500000,...,10000` | Comma-separated cooler resolutions |
| `--hic_base_bin` | `10000` | Base bin size for cooler |
| `--hic_plot_resolutions` | `1000000,500000,250000,100000` | Resolutions to plot contact maps |
| `--hic_balance` | `true` | Balance contact matrices (ICE) |
| `--scaffold_min_size` | `0` | Min scaffold size for contact maps (0 = all) |

### Hi-C Scaffolding (YaHS)

Both round 1 and round 2 have independent parameter sets.

| Parameter (Round 1 / Round 2) | Default R1 / R2 | Description |
|---|---|---|
| `--yahs_round{1,2}_min_contig_length` | `20000` / `100000` | Min contig length for scaffolding |
| `--yahs_round{1,2}_min_mapq` | `10` / `20` | Min MAPQ for Hi-C reads |
| `--yahs_round{1,2}_resolutions` | *(see defaults)* | Resolution ladder for YaHS |
| `--yahs_round{1,2}_rounds_per_resolution` | `null` | Rounds per resolution (null = YaHS default) |
| `--yahs_round{1,2}_enzyme` | `null` | Restriction enzyme (null = enzyme-free) |
| `--yahs_round{1,2}_no_contig_ec` | `false` / `false` | Disable contig error correction in YaHS |
| `--yahs_round{1,2}_no_scaffold_ec` | `false` / `true` | Disable scaffold error correction in YaHS |

| Parameter | Default | Description |
|---|---|---|
| `--run_scaffold_round2` | *auto* | Explicitly enable/disable round 2; auto-enabled when scaffold correction or decontamination is active |

### Decontamination

| Parameter | Default | Description |
|---|---|---|
| `--decon_run_on_contigs` | `true` | Run FCS-GX on contigs |
| `--decon_run_on_scaffolds` | `false` | Run FCS-GX on scaffolds |
| `--decon_source_taxid` | `7898` | Source organism NCBI taxonomy ID (default: Actinopterygii) |
| `--decon_run_fcs_adaptor` | `false` | Run FCS-adaptor for adapter/vector screening |
| `--decon_fcsadaptor_mode` | `euk` | FCS-adaptor mode (`euk` or `prok`) |
| `--decon_container_engine` | `singularity` | Container engine for FCS tools |

### Database Paths

| Parameter | Default | Description |
|---|---|---|
| `--db_base` | `/work/birdlab/databases` | Base directory for all databases |
| `--gxdb_dir` | `${db_base}/fcs-gx` | FCS-GX database directory |
| `--gxdb_profile` | `all` | FCS-GX profile (`all` or `test-only`) |
| `--gxdb_force` | `false` | Force re-download of FCS-GX database |
| `--busco_lineage` | `actinopterygii_odb10` | BUSCO lineage dataset |
| `--busco_downloads` | `/work/birdlab/databases/busco` | Local BUSCO dataset cache |
| `--merqury_k` | `21` | K-mer size for meryl database |

### Finalization and Visualization

| Parameter | Default | Description |
|---|---|---|
| `--make_final_contact_maps` | `true` | Generate final Hi-C contact maps |
| `--run_pairwise_alignments` | `true` | Generate pairwise dotplots |
| `--pairwise_alignment_preset` | `asm5` | minimap2 preset for assembly alignments |
| `--pairwise_alignment_min_mapq` | `5` | Min MAPQ for alignments |
| `--pairwise_alignment_min_aln_bp` | `10000` | Min alignment length (bp) |
| `--pairwise_alignment_mode` | `within_sample` | `within_sample` (hap1 vs hap2) or `all` (all pairs) |
| `--pairwise_dotplot_width` | `10` | Dotplot width (inches) |
| `--pairwise_dotplot_height` | `10` | Dotplot height (inches) |

### Telomere Detection

| Parameter | Default | Description |
|---|---|---|
| `--telomere_window` | `10000` | Window size (bp) at each sequence end |
| `--telomere_min_repeats` | `10` | Min consecutive motif repeats required |

### Compartments and TADs

| Parameter | Default | Description |
|---|---|---|
| `--compartment_resolution` | `250000` | Resolution for A/B compartment calling |
| `--compartment_min_contig_bp` | `5000000` | Min contig size for compartment analysis |
| `--compartment_max_contigs` | `30` | Max contigs to include |
| `--tad_resolution` | `50000` | Resolution for TAD calling |
| `--tad_window_bp` | `500000` | Window size for insulation score |
| `--tad_min_contig_bp` | `5000000` | Min contig size for TAD analysis |
| `--tad_max_contigs` | `0` | Max contigs (0 = no limit) |

### Coverage Book

| Parameter | Default | Description |
|---|---|---|
| `--bin_size` | `1000` | Bin size for coverage calculation |
| `--min_len` | `1000000` | Min scaffold length to include |
| `--min_mapq` | `5` | Min MAPQ for coverage reads |

---

## Output

All outputs are written to `--outdir` (default `./results`) with the following directory structure:

```
results/
├── assembly/
│   ├── contig/
│   │   ├── hifiasm/                    # Raw HIFIASM output (FASTA, GFA, logs)
│   │   ├── purge_dups/                 # Purged assemblies (if enabled)
│   │   ├── misassembly_correction/     # Inspector-corrected contigs
│   │   └── decontamination/            # FCS-GX cleaned contigs
│   └── scaffold/
│       ├── round1/                     # YaHS round 1 scaffolds + AGP
│       ├── misassembly_correction/     # Inspector-corrected scaffolds
│       ├── decontamination/            # FCS-GX cleaned scaffolds
│       └── round2/                     # YaHS round 2 scaffolds + AGP
│
├── assemblies/
│   └── gap_filled/                     # Final gap-filled assemblies (FASTA)
│
├── qc/
│   ├── reads/
│   │   ├── hifi/                       # HiFi read QC (FastQC)
│   │   ├── hic_raw/                    # Raw Hi-C QC (FastQC)
│   │   └── hic_trimmed/               # Trimmed Hi-C QC (FastQC, fastp)
│   ├── assembly/
│   │   ├── intermediate/              # Per-checkpoint QC summaries (TSV)
│   │   └── compiled/                  # Aggregated final QC report
│   └── hic/                           # Hi-C mapping metrics (BAM + pairs stats)
│
├── snail_plots/                        # BlobToolKit snail plots (SVG)
├── contact_maps/                       # Hi-C contact maps (PNG, .mcool)
├── coverage_books/                     # HiFi coverage visualizations
├── pairwise_alignments/                # Dotplot PNGs and PAF files
├── telomeres/                          # Telomere scan results
│
├── reports/
│   ├── pipeline_summary_report.md      # Markdown summary with images
│   └── pipeline_summary_report.html    # HTML summary with embedded images
│
└── pipeline/
    ├── pipeline_timeline.html          # Nextflow execution timeline
    ├── pipeline_report.html            # Nextflow execution report
    ├── pipeline_trace.txt              # Process-level resource trace
    └── pipeline_dag.png                # Workflow DAG visualization
```

### Key Output Files

| File | Description |
|---|---|
| `assemblies/gap_filled/{sample}_hap{1,2}_filled.fasta` | Final gap-filled haplotype assemblies |
| `qc/assembly/compiled/final_qc_report.tsv` | All QC metrics across all checkpoints |
| `snail_plots/{sample}_hap{1,2}_gap_filled_snail.svg` | Assembly quality visualization |
| `contact_maps/{sample}_hap{1,2}_final_*bp_contact_map.png` | Hi-C contact maps at various resolutions |
| `pairwise_alignments/{sample}_hap1_vs_{sample}_hap2_dotplot.png` | Haplotype comparison dotplots |
| `reports/pipeline_summary_report.html` | Visual summary of all samples |

---

## Execution Profiles

The pipeline ships with three profiles defined in `nextflow.config`:

```bash
# SLURM cluster (default for HPC)
nextflow run main.nf -profile slurm ...

# Local execution
nextflow run main.nf -profile local ...

# Docker
nextflow run main.nf -profile docker ...
```

The SLURM profile targets a cluster with `normal` nodes (64 CPUs, 343 GB RAM), `bigmem` nodes (64 CPUs, 700 GB RAM), and `ultramem` nodes (64 CPUs, 1416 GB RAM). Processes that may exceed 343 GB (e.g., FCS-GX screening, DIAMOND BLASTX) are automatically routed to higher-memory queues. Dynamic retry logic escalates memory allocation on failure.

---

## Resource Configuration

Resource allocation is tiered by process label in `nextflow.config`. Key allocations:

| Process | CPUs | Memory (attempt 1 → retry) | Queue |
|---|---|---|---|
| HIFIASM | 48 | 100 GB → 300 GB | normal → bigmem |
| Inspector | 48 | 110 GB → 300 GB | normal → bigmem |
| BWA-MEM2 (Hi-C) | 48 | 120 GB | normal |
| BUSCO | 48 | 55 GB → 100 GB | normal → bigmem |
| BUILD_MERYL_DB | 16 | 100 GB → 200 GB | normal → bigmem |
| Gap Filling | 12 | 150 GB → 343 GB | normal → bigmem |
| FCS-GX | 32 | 500 GB | ultramem/bigmem |
| YaHS | 1 | 8 GB | normal |
| QUAST | 4 | 16 GB | normal |
| MERQURY | 8 | 20 GB | normal |
| fastp | 16 | 4 GB | normal |

Global limits: `--max_cpus 64`, `--max_memory '700 GB'`, `--max_time '240.h'`.

---

## Caching and Resume

The pipeline is designed for robust `-resume` behavior:

- Each sample processes completely independently — adding new samples does not invalidate caches for existing samples.
- The meryl database is built once per sample and reused across all assembly QC stages.
- Decontamination databases are downloaded once and shared across all samples.
- All parameters are flattened to avoid Nextflow caching issues with nested maps.

```bash
# Resume after failure or with additional samples
nextflow run main.nf \
  --sample_sheet updated_samples.csv \
  --outdir ./results \
  -profile slurm \
  -resume
```

---

## License

*Add license information here.*