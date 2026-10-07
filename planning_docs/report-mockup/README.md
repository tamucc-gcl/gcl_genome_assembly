# Genome Assembly Pipeline Report

> **Layout example:** existing run-report sections are preserved below, with example first-pass chimera evidence integrated. General QC values and plots illustrate the report layout; they are not results from the current evidence job.

> Generated: **2026-10-07**

## Run Summary

- **Pipeline:** gcl_genome_assembly unknown (unknown)
- **Nextflow:** 23.10.1  ·  **Profile:** slurm  ·  **Run:** astonishing_visvesvaraya  ·  **Started:** 2026-08-19T09:39:27.648509263-05:00
- **Samples:** 5  ·  **Inputs (samples with):** HiFi 5 · Hi-C 5 · short-read 0
- **Assemblers:** hifiasm (5)
- **QC stages present:** ctg.base → ctg.purged → ctg.cor → ctg.deco → gap_fill → teloclip → final → ctg.org → scaf.base → scaf.cor → scaf2

---

> **Review needed — possible chromosome joins:** CTlk hap1 scaffold_1 (chr9/chr7/chr12) and scaffold_5 (chr4/chr14); CTlk hap2 scaffold_5 (chr7/chr12) and scaffold_7 (chr14/chr4). Eight sampled boundaries require review; some sample the same chromosome combination. [See candidate locations, evidence and cut instructions](#7-chimera-detection-and-review).

## Table of Contents

1. [Final Genome Assemblies](#1-final-genome-assemblies)
2. [Assembly QC Summary](#2-assembly-qc-summary)
3. [Visual Summary](#3-visual-summary)
4. [Assembly QC Across Pipeline Stages](#4-assembly-qc-across-pipeline-stages)
5. [Mitochondrial Genome](#5-mitochondrial-genome)
6. [Telomere Detection](#6-telomere-detection)
7. [Chimera Detection and Review](#7-chimera-detection-and-review)
8. [Pairwise Alignment Summary](#8-pairwise-alignment-summary)
9. [Pangenome](#9-pangenome)
10. [Methods and Citations](#10-methods-and-citations)

---

## 1. Final Genome Assemblies

Final gap-filled, chromosome-level assemblies produced by the [gcl_genome_assembly](https://github.com/tamucc-gcl/gcl_genome_assembly) pipeline.

### Key Resources

- **Compiled QC data (CSV):** [assembly_qc_metrics.csv](qc/assembly/assembly_qc_metrics.csv)
- **Interactive QC report (HTML):** [assembly_qc_report.html](reports/assembly_qc_report.html) *(download and open in browser)*

### Assembly Files

| Sample | Haplotype 1 | Haplotype 2 |
| --- | --- | --- |
| Sde-CBau_104 | [Sde-CBau_104_hap1.fasta](https://crest-files.tamucc.edu/genomics/sde_genome_assembly/Sde-CBau_104_hap1.fasta) | [Sde-CBau_104_hap2.fasta](https://crest-files.tamucc.edu/genomics/sde_genome_assembly/Sde-CBau_104_hap2.fasta) |
| Sde-CLim_110 | [Sde-CLim_110_hap1.fasta](https://crest-files.tamucc.edu/genomics/sde_genome_assembly/Sde-CLim_110_hap1.fasta) | [Sde-CLim_110_hap2.fasta](https://crest-files.tamucc.edu/genomics/sde_genome_assembly/Sde-CLim_110_hap2.fasta) |
| Sde-CMat_203 | [Sde-CMat_203_hap1.fasta](https://crest-files.tamucc.edu/genomics/sde_genome_assembly/Sde-CMat_203_hap1.fasta) | [Sde-CMat_203_hap2.fasta](https://crest-files.tamucc.edu/genomics/sde_genome_assembly/Sde-CMat_203_hap2.fasta) |
| Sde-CPla_115 | [Sde-CPla_115_hap1.fasta](https://crest-files.tamucc.edu/genomics/sde_genome_assembly/Sde-CPla_115_hap1.fasta) | [Sde-CPla_115_hap2.fasta](https://crest-files.tamucc.edu/genomics/sde_genome_assembly/Sde-CPla_115_hap2.fasta) |
| Sde-CTlk_104 | [Sde-CTlk_104_hap1.fasta](https://crest-files.tamucc.edu/genomics/sde_genome_assembly/Sde-CTlk_104_hap1.fasta) | [Sde-CTlk_104_hap2.fasta](https://crest-files.tamucc.edu/genomics/sde_genome_assembly/Sde-CTlk_104_hap2.fasta) |

### Sample Taxonomy & Genome Profile

Per-sample organism identity (resolved from the NCBI taxid) and k-mer-based
estimates from [GenomeScope2](https://github.com/tbenavi1/genomescope2.0): haploid
genome size, heterozygosity, and repeat content.

| Sample | Species | Taxid | Kingdom | BUSCO Lineage | Est. Genome Size | Heterozygosity | Repeat Content |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Sde-CBau_104 | Spratelloides delicatulus | 373251 | animal | actinopterygii_odb10 | 1.39 Gb | 1.48% | 70.2% |
| Sde-CLim_110 | Spratelloides delicatulus | 373251 | animal | actinopterygii_odb10 | 0.75 Gb | 2.95% | 37.7% |
| Sde-CMat_203 | Spratelloides delicatulus | 373251 | animal | actinopterygii_odb10 | 1.44 Gb | 1.12% | 71.0% |
| Sde-CPla_115 | Spratelloides delicatulus | 373251 | animal | actinopterygii_odb10 | 0.77 Gb | 2.91% | 36.7% |
| Sde-CTlk_104 | Spratelloides delicatulus | 373251 | animal | actinopterygii_odb10 | 0.80 Gb | 2.83% | 36.6% |

## 2. Assembly QC Summary

### Overview (Final Assemblies)

| Sample | Status | Total Length | Est. Genome Size | % of Estimate | Largest Scaffold | auN | auN Pieces | L90 | GC (%) | Coverage | BUSCO Complete | QV | K-mer Completeness |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Sde-CBau_104 | ✅ | 1.08 Gb / 1.08 Gb | 1.39 Gb | 78% / 78% | 95.7 Mb / 95.0 Mb | 69.2 Mb / 64.6 Mb | 15.6 / 16.7 | 15 / 16 | 43.2% / 43.2% | 40.3 / 40.2 | 92.7% / 93.1% | 62.8 | 98.6% |
| Sde-CLim_110 | ✅ | 1.00 Gb / 1.04 Gb | 0.75 Gb | 134% / 139% | 83.7 Mb / 89.5 Mb | 60.5 Mb / 63.8 Mb | 16.6 / 16.4 | 15 / 15 | 43.2% / 43.2% | 25.8 / 25.3 | 90.6% / 91.8% | 60.4 | 97.5% |
| Sde-CMat_203 | ✅ | 1.06 Gb / 1.06 Gb | 1.44 Gb | 74% / 73% | 91.6 Mb / 94.7 Mb | 74.4 Mb / 74.7 Mb | 14.3 / 14.1 | 13 / 13 | 43.2% / 43.2% | 39.5 / 39.7 | 92.6% / 92.7% | 63.0 | 98.6% |
| Sde-CPla_115 | ⚠️ | 1.06 Gb / 0.98 Gb | 0.77 Gb | 137% / 126% | 14.3 Mb / 15.8 Mb | 2.9 Mb / 2.9 Mb | 360.2 / 334.1 | 739 / 654 | 43.2% / 43.2% | 23.3 / 24.5 | 90.9% / 89.5% | 60.2 | 97.1% |
| Sde-CTlk_104 | ✅ | 1.02 Gb / 1.03 Gb | 0.80 Gb | 128% / 129% | 111.6 Mb / 88.9 Mb | 61.2 Mb / 55.6 Mb | 16.7 / 18.6 | 17 / 18 | 43.3% / 43.2% | 26.4 / 26.5 | 91.0% / 91.9% | 60.8 | 98.0% |

> **Status:** ✅ all checks pass · ⚠️ below threshold (BUSCO complete ≥ 90%, QV ≥ 40, k-mer completeness ≥ 90%, auN pieces ≤ 30 (2 × 15 chromosomes) for hifi+hic, largest scaffold ≥ 0.5 × mean chromosome). Thresholds are per-run tunable.

Flagged: **Sde-CPla_115** (BUSCO<90%, auN-pieces 360 > 30, largest 0.20× mean chrom, cohort outlier 20.0 MAD)

> **BUSCO lineage:** completeness is scored per sample against the lineage resolved from its NCBI taxid (section 1); where a specific clade can't be resolved it falls back to `params.busco_lineage` (here `eukaryota_odb10`). All samples scored against `actinopterygii_odb10`.

### Contiguity Detail

| Sample | Tier | Expected chrom | auN Pieces (worst hap) | Gate | Largest / mean chrom | Cohort MAD z |
| --- | --- | --- | --- | --- | --- | --- |
| Sde-CBau_104 | hifi+hic | 15 (name map) | 16.7 | ≤ 30 | 1.32× | 0.0 |
| Sde-CLim_110 | hifi+hic | 15 (name map) | 16.6 | ≤ 30 | 1.25× | -0.1 |
| Sde-CMat_203 | hifi+hic | 15 (name map) | 14.3 | ≤ 30 | 1.29× | -1.0 |
| Sde-CPla_115 | hifi+hic | 15 (name map) | 360.2 | ≤ 30 | 0.20× | 20.0 |
| Sde-CTlk_104 | hifi+hic | 15 (name map) | 18.6 | ≤ 30 | 1.29× | 0.7 |

> **auN pieces** = total length / auN, the length-weighted effective piece count (exactly n for n equal pieces); worst haplotype per sample. Chromosome-scale gates apply only to the HiFi+Hi-C tier; cohort bar within (taxid, tier) = median + max(3 × MAD, log10(1.5)) of log10(auN pieces) = 26.5 pieces.

Advisory — largest scaffold above 1.5× mean chromosome, i.e. a possible over-join that a worst-haplotype gate cannot see: **Sde-CTlk_104** (1.64×).

## 3. Visual Summary

Each row shows one sample with snail plots, contact maps, the within-sample
haplotype-vs-haplotype dotplot, and riparian (ribbon) synteny plot.

<table>
<tr>
  <th>Sample</th>
  <th>Hap1 Snail</th>
  <th>Hap2 Snail</th>
  <th>Hap1 Contact Map</th>
  <th>Hap2 Contact Map</th>
  <th>Hap1 vs Hap2 Dotplot</th>
  <th>Hap1 vs Hap2 Riparian</th>
</tr>
<tr>
  <td><b>Sde-CBau_104</b></td>
  <td><img src="snail_plots/Sde-CBau_104_hap1_final_snail.svg" alt="Sde-CBau_104_hap1" width="500"></td>
  <td><img src="snail_plots/Sde-CBau_104_hap2_final_snail.svg" alt="Sde-CBau_104_hap2" width="500"></td>
  <td><img src="contact_maps/Sde-CBau_104_hap1_final_1000000bp_contact_map.png" alt="Sde-CBau_104_hap1" width="500"></td>
  <td><img src="contact_maps/Sde-CBau_104_hap2_final_1000000bp_contact_map.png" alt="Sde-CBau_104_hap2" width="500"></td>
  <td><img src="pairwise_alignments/Sde-CBau_104_hap1_vs_Sde-CBau_104_hap2_dotplot.png" alt="Sde-CBau_104_hap1 vs Sde-CBau_104_hap2" width="500"></td>
  <td><img src="pairwise_alignments/Sde-CBau_104_hap1_vs_Sde-CBau_104_hap2_riparian.png" alt="Sde-CBau_104_hap1 vs Sde-CBau_104_hap2 riparian" width="500"></td>
</tr>
<tr>
  <td><b>Sde-CLim_110</b></td>
  <td><img src="snail_plots/Sde-CLim_110_hap1_final_snail.svg" alt="Sde-CLim_110_hap1" width="500"></td>
  <td><img src="snail_plots/Sde-CLim_110_hap2_final_snail.svg" alt="Sde-CLim_110_hap2" width="500"></td>
  <td><img src="contact_maps/Sde-CLim_110_hap1_final_1000000bp_contact_map.png" alt="Sde-CLim_110_hap1" width="500"></td>
  <td><img src="contact_maps/Sde-CLim_110_hap2_final_1000000bp_contact_map.png" alt="Sde-CLim_110_hap2" width="500"></td>
  <td><img src="pairwise_alignments/Sde-CLim_110_hap1_vs_Sde-CLim_110_hap2_dotplot.png" alt="Sde-CLim_110_hap1 vs Sde-CLim_110_hap2" width="500"></td>
  <td><img src="pairwise_alignments/Sde-CLim_110_hap1_vs_Sde-CLim_110_hap2_riparian.png" alt="Sde-CLim_110_hap1 vs Sde-CLim_110_hap2 riparian" width="500"></td>
</tr>
<tr>
  <td><b>Sde-CMat_203</b></td>
  <td><img src="snail_plots/Sde-CMat_203_hap1_final_snail.svg" alt="Sde-CMat_203_hap1" width="500"></td>
  <td><img src="snail_plots/Sde-CMat_203_hap2_final_snail.svg" alt="Sde-CMat_203_hap2" width="500"></td>
  <td><img src="contact_maps/Sde-CMat_203_hap1_final_1000000bp_contact_map.png" alt="Sde-CMat_203_hap1" width="500"></td>
  <td><img src="contact_maps/Sde-CMat_203_hap2_final_1000000bp_contact_map.png" alt="Sde-CMat_203_hap2" width="500"></td>
  <td><img src="pairwise_alignments/Sde-CMat_203_hap1_vs_Sde-CMat_203_hap2_dotplot.png" alt="Sde-CMat_203_hap1 vs Sde-CMat_203_hap2" width="500"></td>
  <td><img src="pairwise_alignments/Sde-CMat_203_hap1_vs_Sde-CMat_203_hap2_riparian.png" alt="Sde-CMat_203_hap1 vs Sde-CMat_203_hap2 riparian" width="500"></td>
</tr>
<tr>
  <td><b>Sde-CPla_115</b></td>
  <td><img src="snail_plots/Sde-CPla_115_hap1_final_snail.svg" alt="Sde-CPla_115_hap1" width="500"></td>
  <td><img src="snail_plots/Sde-CPla_115_hap2_final_snail.svg" alt="Sde-CPla_115_hap2" width="500"></td>
  <td><img src="contact_maps/Sde-CPla_115_hap1_final_1000000bp_contact_map.png" alt="Sde-CPla_115_hap1" width="500"></td>
  <td><img src="contact_maps/Sde-CPla_115_hap2_final_1000000bp_contact_map.png" alt="Sde-CPla_115_hap2" width="500"></td>
  <td><img src="pairwise_alignments/Sde-CPla_115_hap1_vs_Sde-CPla_115_hap2_dotplot.png" alt="Sde-CPla_115_hap1 vs Sde-CPla_115_hap2" width="500"></td>
  <td><img src="pairwise_alignments/Sde-CPla_115_hap1_vs_Sde-CPla_115_hap2_riparian.png" alt="Sde-CPla_115_hap1 vs Sde-CPla_115_hap2 riparian" width="500"></td>
</tr>
<tr>
  <td><b>Sde-CTlk_104</b></td>
  <td><img src="snail_plots/Sde-CTlk_104_hap1_final_snail.svg" alt="Sde-CTlk_104_hap1" width="500"></td>
  <td><img src="snail_plots/Sde-CTlk_104_hap2_final_snail.svg" alt="Sde-CTlk_104_hap2" width="500"></td>
  <td><img src="contact_maps/Sde-CTlk_104_hap1_final_1000000bp_contact_map.png" alt="Sde-CTlk_104_hap1" width="500"></td>
  <td><img src="contact_maps/Sde-CTlk_104_hap2_final_1000000bp_contact_map.png" alt="Sde-CTlk_104_hap2" width="500"></td>
  <td><img src="pairwise_alignments/Sde-CTlk_104_hap1_vs_Sde-CTlk_104_hap2_dotplot.png" alt="Sde-CTlk_104_hap1 vs Sde-CTlk_104_hap2" width="500"></td>
  <td><img src="pairwise_alignments/Sde-CTlk_104_hap1_vs_Sde-CTlk_104_hap2_riparian.png" alt="Sde-CTlk_104_hap1 vs Sde-CTlk_104_hap2 riparian" width="500"></td>
</tr>
</table>

All contact map resolutions: [contact_maps/](contact_maps/)

Additional cross-sample dotplots (40): [pairwise_alignments/](pairwise_alignments/)

### GenomeScope Profiles

K-mer spectra with the fitted [GenomeScope2](https://github.com/tbenavi1/genomescope2.0)
model — visual support for the haploid genome-size estimate (peak position), plus
heterozygosity and repeat content. One profile per sample.

<table>
<tr><th>Sample</th><th>GenomeScope Profile</th></tr>
<tr><td><b>Sde-CBau_104</b></td><td><img src="est_genome_size/Sde-CBau_104.linear_plot.png" alt="Sde-CBau_104" width="500"></td></tr>
<tr><td><b>Sde-CLim_110</b></td><td><img src="est_genome_size/Sde-CLim_110.linear_plot.png" alt="Sde-CLim_110" width="500"></td></tr>
<tr><td><b>Sde-CMat_203</b></td><td><img src="est_genome_size/Sde-CMat_203.linear_plot.png" alt="Sde-CMat_203" width="500"></td></tr>
<tr><td><b>Sde-CPla_115</b></td><td><img src="est_genome_size/Sde-CPla_115.linear_plot.png" alt="Sde-CPla_115" width="500"></td></tr>
<tr><td><b>Sde-CTlk_104</b></td><td><img src="est_genome_size/Sde-CTlk_104.linear_plot.png" alt="Sde-CTlk_104" width="500"></td></tr>
</table>

## 4. Assembly QC Across Pipeline Stages

<details>
<summary><b>Click to expand: QC trend plots and cross-stage comparison</b></summary>

### QC Trend Plots

#### BUSCO Completeness

<img src="qc/assembly/busco.png" alt="BUSCO Completeness" width="800">

#### K-mer QV & Completeness

<img src="qc/assembly/kmer.png" alt="K-mer QV & Completeness" width="800">

#### QUAST Assembly Metrics

<img src="qc/assembly/quast_misc.png" alt="QUAST Assembly Metrics" width="800">

#### Contig/Scaffold Counts

<img src="qc/assembly/contig_count.png" alt="Contig/Scaffold Counts" width="800">

#### Assembly Size

<img src="qc/assembly/contig_length.png" alt="Assembly Size" width="800">

#### Hi-C Trans:Cis Ratio

<img src="qc/assembly/trans_cis.png" alt="Hi-C Trans:Cis Ratio" width="800">


</details>

## 5. Mitochondrial Genome

<table>
<tr>
  <th>Sample</th>
  <th>Gene Map</th>
  <th>Length (bp)</th>
  <th>Circular</th>
  <th>Genes</th>
  <th>tRNAs</th>
  <th>rRNAs</th>
  <th>GenBank</th>
</tr>
<tr><td><b>Sde-CLim_110</b></td><td><img src="organelle/Sde-CLim_110_mito_circular.png" alt="Mito gene map: Sde-CLim_110" width="500"></td><td>16,618</td><td>Yes</td><td>26</td><td>44</td><td>2</td><td><a href="organelle/Sde-CLim_110_mitogenome.gb">Sde-CLim_110_mitogenome.gb</a></td></tr>
<tr><td><b>Sde-CBau_104</b></td><td><img src="organelle/Sde-CBau_104_mito_circular.png" alt="Mito gene map: Sde-CBau_104" width="500"></td><td>16,618</td><td>Yes</td><td>26</td><td>44</td><td>2</td><td><a href="organelle/Sde-CBau_104_mitogenome.gb">Sde-CBau_104_mitogenome.gb</a></td></tr>
<tr><td><b>Sde-CMat_203</b></td><td><img src="organelle/Sde-CMat_203_mito_circular.png" alt="Mito gene map: Sde-CMat_203" width="500"></td><td>16,616</td><td>Yes</td><td>26</td><td>44</td><td>2</td><td><a href="organelle/Sde-CMat_203_mitogenome.gb">Sde-CMat_203_mitogenome.gb</a></td></tr>
<tr><td><b>Sde-CPla_115</b></td><td><img src="organelle/Sde-CPla_115_mito_circular.png" alt="Mito gene map: Sde-CPla_115" width="500"></td><td>16,617</td><td>Yes</td><td>26</td><td>44</td><td>2</td><td><a href="organelle/Sde-CPla_115_mitogenome.gb">Sde-CPla_115_mitogenome.gb</a></td></tr>
<tr><td><b>Sde-CTlk_104</b></td><td><img src="organelle/Sde-CTlk_104_mito_circular.png" alt="Mito gene map: Sde-CTlk_104" width="500"></td><td>16,619</td><td>Yes</td><td>26</td><td>44</td><td>2</td><td><a href="organelle/Sde-CTlk_104_mitogenome.gb">Sde-CTlk_104_mitogenome.gb</a></td></tr>
</table>

Mitochondrial contigs were identified and removed from each haplotype
assembly prior to purge\_dups, Inspector, decontamination, and scaffolding.

## 6. Telomere Detection

#### Telomere Presence Summary (tidk)

Telomeric repeat detection across scaffold ends using
[tidk](https://github.com/tolkit/telomeric-identifier).

| Haplotype | Scaffolds | Total Length | 5' Telo | 3' Telo | Both | 5' Only | 3' Only | None | % with Telomere |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Sde-CBau_104_hap1 |  172 | 1,078,241,830 | 14 | 14 | 6 |  8 |  8 |  150 | 12.79% |
| Sde-CBau_104_hap2 |  174 | 1,082,713,338 | 16 | 15 | 4 | 12 | 11 |  147 | 15.52% |
| Sde-CLim_110_hap1 | 1050 | 1,012,798,437 | 30 |  8 | 4 | 26 |  4 | 1016 | 3.24% |
| Sde-CLim_110_hap2 |  999 | 1,054,354,184 | 28 | 19 | 2 | 26 | 17 |  954 | 4.5% |
| Sde-CMat_203_hap1 |  115 | 1,062,920,163 | 12 | 14 | 8 |  4 |  6 |   97 | 15.65% |
| Sde-CMat_203_hap2 |  119 | 1,057,380,794 | 12 |  6 | 1 | 11 |  5 |  102 | 14.29% |
| Sde-CPla_115_hap1 | 2177 | 1,082,639,583 | 42 | 18 | 0 | 42 | 18 | 2117 | 2.76% |
| Sde-CPla_115_hap2 | 1908 | 995,810,437 | 28 | 12 | 4 | 24 |  8 | 1872 | 1.89% |
| Sde-CTlk_104_hap1 |  532 | 1,027,194,400 | 21 | 13 | 1 | 20 | 12 |  499 | 6.2% |
| Sde-CTlk_104_hap2 |  441 | 1,035,959,456 | 27 | 14 | 2 | 25 | 12 |  402 | 8.84% |

## 7. Chimera Detection and Review

### Assessment Summary

All ten haplotype assemblies were assessed for chromosome-scale composite scaffolds and local join support. The other assemblies in this run supply comparison evidence for each assessed assembly. Eight candidate boundaries were flagged in four CTlk scaffolds; six are in hap1 and two in hap2. No candidates were flagged in CBau, CLim, CMat or CPla in this example.

A candidate is a boundary to examine, not a confirmed misassembly. Multiple boundaries on one scaffold can investigate the same suspected chromosome combination. An absence of detected candidates does not establish that an assembly is free of misjoins.

**First-pass status:** evidence review is pending. All cut-file rows are selected=NO; no cuts have been applied. These assemblies retain their assessed joins.

### Candidate Boundaries

| Candidate | Assembly | Scaffold | Chromosomes left → right | Region to review, bp | Exact cut, bp | Review priority | For cutting | Against cutting |
| --- | --- | --- | --- | --- | ---: | --- | --- | --- |
| [C01](Sde-CTlk_104_hap1.review/report.md#candidate-c01) | Sde-CTlk_104_hap1 | scaffold_1 | chr9 → chr7 | 62618168–62941917 | Not assigned | Investigate chromosome transition; exact cut not localized | Separate chromosomes: Sde-CMat_203: chr9 → chr7 | No opposing chromosome evidence observed |
| [C02](Sde-CTlk_104_hap1.review/report.md#candidate-c02) | Sde-CTlk_104_hap1 | scaffold_1 | chr9 → chr7 | 63051225–63051325 | 63051325 | Prioritize gap-cut review | Separate chromosomes: Sde-CBau_104: chr9 → chr7; Sde-CMat_203: chr9 → chr7 | No opposing chromosome evidence observed |
| [C03](Sde-CTlk_104_hap1.review/report.md#candidate-c03) | Sde-CTlk_104_hap1 | scaffold_1 | Unresolved from qualified local alignments | 95877814–96047526 | Not assigned | Insufficient evidence to propose a break | No informative independent chromosome evidence supporting a break | No opposing chromosome evidence observed |
| [C04](Sde-CTlk_104_hap1.review/report.md#candidate-c04) | Sde-CTlk_104_hap1 | scaffold_1 | chr12 → chr12 | 116057447–116057547 | 116057547 | Evidence favors retaining this sampled boundary | No informative independent chromosome evidence supporting a break | Same chromosome: Sde-CBau_104: chr12 → chr12; Sde-CMat_203: chr12 → chr12 |
| [C05](Sde-CTlk_104_hap1.review/report.md#candidate-c05) | Sde-CTlk_104_hap1 | scaffold_1 | chr12 → chr12 | 116171488–116171588 | 116171588 | Evidence favors retaining this sampled boundary | No informative independent chromosome evidence supporting a break | Same chromosome: Sde-CBau_104: chr12 → chr12; Sde-CLim_110: chr12 → chr12; Sde-CMat_203: chr12 → chr12 |
| [C06](Sde-CTlk_104_hap1.review/report.md#candidate-c06) | Sde-CTlk_104_hap1 | scaffold_5 | chr4 → chr14 | 35443804–35581476 | Not assigned | Investigate chromosome transition; exact cut not localized | Separate chromosomes: Sde-CLim_110: chr4 → chr14; Sde-CMat_203: chr4 → chr14 | No opposing chromosome evidence observed |
| [C01](Sde-CTlk_104_hap2.review/report.md#candidate-c01) | Sde-CTlk_104_hap2 | scaffold_5 | Unresolved from qualified local alignments | 43055156–43233354 | Not assigned | Insufficient evidence to propose a break | No informative independent chromosome evidence supporting a break | No opposing chromosome evidence observed |
| [C02](Sde-CTlk_104_hap2.review/report.md#candidate-c02) | Sde-CTlk_104_hap2 | scaffold_7 | chr14 → chr4 | 39549503–39688365 | Not assigned | Investigate chromosome transition; exact cut not localized | Separate chromosomes: Sde-CBau_104: chr14 → chr4; Sde-CLim_110: chr14 → chr4; Sde-CMat_203: chr14 → chr4 | No opposing chromosome evidence observed |

IDs are scoped to an assembly. Review priorities organize evidence, not cut selections. Different chromosome labels mean the sampled sides map to separate chromosomes in other individuals. Same chromosome labels oppose a fusion at this boundary.

### Evidence Available for Review

The [detailed Markdown evidence report](chimera-evidence.md) presents each candidate's peer chromosome assignments, aligned coverage and ambiguity, immediate HiFi flank/spanning measurements, graph context, and Hi-C contact counts and control plots for each library.

Independent comparisons are grouped by individual. The second haplotype of an individual is not a second independent individual. The other CTlk haplotype provides within-individual context; matching chromosome combinations are biologically relevant but do not prove the local join. A poorly observed boundary or absent alignment is reported as uninformative rather than as evidence of a misjoin.

For the sampled chr9/chr7 boundary, two informative independent individuals place its sides on different chromosomes. Two farther sampled boundaries in scaffold_1 map to chr12 on both sides in informative peers. The chr4/chr14 transitions have informative chromosome-disagreement evidence but lack an exact localized internal cut. These findings require different review decisions despite appearing in the same candidate inventory.

### Record Decisions and Request Cuts

- **Keep or leave unresolved:** leave selected=NO in the [review TSV](chimera_review.tsv) and record reviewer and reason when the boundary has been reviewed. A completed keep decision differs from a pending review even though neither requests a cut.
- **Cut a detected boundary:** set selected=YES, verify the action and exact coordinates, and supply reviewer and reason. Internal cuts also require localization_status=localized.
- **Add an undetected cut:** append a new row for any assessed assembly using its original scaffold ID and checksum from the [assembly registry](assembly-registry.tsv). Supply an exact coordinate and supporting evidence. [Instructions and row template](cut-interface.md#add-a-cut-the-detector-did-not-find).

Copy the TSV outside the generated results before editing. Repeat the original launch with `--chimera_break /absolute/path/reviewed-cuts.tsv -resume`. The pipeline validates the source FASTA and cut coordinates, preserves every base, writes a coordinate lift and reconstruction audit, and repeats chromosome assignment for affected cohorts.

[Full evidence report](chimera-evidence.md) · [Editable review file](chimera_review.tsv) · [Cut-file instructions](cut-interface.md)

## 8. Pairwise Alignment Summary

<details>
<summary><b>Click to expand: Pairwise alignment metrics</b></summary>

| ref_id | qry_id | preset | min_mapq | min_aln_bp | n_align | sum_aln_bp | sum_match_bp | mean_identity | mean_dv | mean_mapq |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Sde-CBau_104_hap1 | Sde-CBau_104_hap2 | asm5 | 5 | 150000 | 227 | 1187204940 | 669473560 | 0.988421 | 0.0115789 | 58.1366 |
| Sde-CBau_104_hap1 | Sde-CLim_110_hap1 | asm5 | 5 | 150000 | 437 | 1048418641 | 584873765 | 0.988603 | 0.0113971 | 59.3272 |
| Sde-CBau_104_hap1 | Sde-CLim_110_hap2 | asm5 | 5 | 150000 | 446 | 1096073744 | 613815745 | 0.98818 | 0.0118196 | 58.9664 |
| Sde-CBau_104_hap1 | Sde-CMat_203_hap1 | asm5 | 5 | 150000 | 221 | 1173642403 | 656319169 | 0.988091 | 0.0119094 | 58.1629 |
| Sde-CBau_104_hap1 | Sde-CMat_203_hap2 | asm5 | 5 | 150000 | 215 | 1169652324 | 655665529 | 0.988352 | 0.011648 | 57.1721 |
| Sde-CBau_104_hap1 | Sde-CPla_115_hap1 | asm5 | 5 | 150000 | 1232 | 1052974637 | 605382160 | 0.988686 | 0.0113137 | 59.6437 |
| Sde-CBau_104_hap1 | Sde-CPla_115_hap2 | asm5 | 5 | 150000 | 1120 | 988326926 | 579828915 | 0.98992 | 0.0100798 | 59.5402 |
| Sde-CBau_104_hap1 | Sde-CTlk_104_hap1 | asm5 | 5 | 150000 | 343 | 1106197009 | 605451791 | 0.987756 | 0.0122441 | 59.2595 |
| Sde-CBau_104_hap1 | Sde-CTlk_104_hap2 | asm5 | 5 | 150000 | 315 | 1129076066 | 621603061 | 0.988033 | 0.0119671 | 58.9143 |
| Sde-CBau_104_hap2 | Sde-CLim_110_hap1 | asm5 | 5 | 150000 | 449 | 1050839120 | 585326797 | 0.988231 | 0.011769 | 59.118 |
| Sde-CBau_104_hap2 | Sde-CLim_110_hap2 | asm5 | 5 | 150000 | 449 | 1090285943 | 609906937 | 0.988343 | 0.0116572 | 58.902 |
| Sde-CBau_104_hap2 | Sde-CMat_203_hap1 | asm5 | 5 | 150000 | 211 | 1167797173 | 655424791 | 0.988314 | 0.0116862 | 58.8531 |
| Sde-CBau_104_hap2 | Sde-CMat_203_hap2 | asm5 | 5 | 150000 | 197 | 1164115176 | 651772524 | 0.988361 | 0.0116391 | 59.3858 |
| Sde-CBau_104_hap2 | Sde-CPla_115_hap1 | asm5 | 5 | 150000 | 1222 | 1052642721 | 606267016 | 0.988645 | 0.0113552 | 59.8642 |
| Sde-CBau_104_hap2 | Sde-CPla_115_hap2 | asm5 | 5 | 150000 | 1105 | 987180393 | 576660909 | 0.989748 | 0.0102519 | 59.8914 |
| Sde-CBau_104_hap2 | Sde-CTlk_104_hap1 | asm5 | 5 | 150000 | 353 | 1105746718 | 604323321 | 0.987682 | 0.0123185 | 58.9717 |
| Sde-CBau_104_hap2 | Sde-CTlk_104_hap2 | asm5 | 5 | 150000 | 317 | 1118463501 | 626136546 | 0.988356 | 0.0116436 | 59.5268 |
| Sde-CLim_110_hap1 | Sde-CLim_110_hap2 | asm5 | 5 | 150000 | 500 | 996949229 | 564823077 | 0.988291 | 0.0117092 | 59.028 |
| Sde-CLim_110_hap1 | Sde-CMat_203_hap1 | asm5 | 5 | 150000 | 439 | 1045735559 | 591024681 | 0.987772 | 0.0122279 | 59.2916 |
| Sde-CLim_110_hap1 | Sde-CMat_203_hap2 | asm5 | 5 | 150000 | 438 | 1045124100 | 592419895 | 0.988169 | 0.0118312 | 59.1142 |
| Sde-CLim_110_hap1 | Sde-CPla_115_hap1 | asm5 | 5 | 150000 | 1120 | 951234026 | 554734675 | 0.988917 | 0.0110834 | 59.8786 |
| Sde-CLim_110_hap1 | Sde-CPla_115_hap2 | asm5 | 5 | 150000 | 1032 | 904811401 | 535043497 | 0.989916 | 0.0100842 | 59.874 |
| Sde-CLim_110_hap1 | Sde-CTlk_104_hap1 | asm5 | 5 | 150000 | 465 | 1010013806 | 550206210 | 0.987102 | 0.0128977 | 59.557 |
| Sde-CLim_110_hap1 | Sde-CTlk_104_hap2 | asm5 | 5 | 150000 | 467 | 1020875851 | 558989511 | 0.987278 | 0.0127221 | 59.394 |
| Sde-CLim_110_hap2 | Sde-CMat_203_hap1 | asm5 | 5 | 150000 | 469 | 1075052204 | 614157166 | 0.988141 | 0.0118594 | 59.3198 |
| Sde-CLim_110_hap2 | Sde-CMat_203_hap2 | asm5 | 5 | 150000 | 443 | 1072962672 | 611883427 | 0.988135 | 0.0118646 | 59.2867 |
| Sde-CLim_110_hap2 | Sde-CPla_115_hap1 | asm5 | 5 | 150000 | 1129 | 984110266 | 578284734 | 0.988928 | 0.011072 | 59.8149 |
| Sde-CLim_110_hap2 | Sde-CPla_115_hap2 | asm5 | 5 | 150000 | 1045 | 927587416 | 551017266 | 0.990139 | 0.00986098 | 59.8268 |
| Sde-CLim_110_hap2 | Sde-CTlk_104_hap1 | asm5 | 5 | 150000 | 485 | 1032031406 | 565157476 | 0.987238 | 0.0127617 | 59.101 |
| Sde-CLim_110_hap2 | Sde-CTlk_104_hap2 | asm5 | 5 | 150000 | 470 | 1051833904 | 575764587 | 0.987275 | 0.0127252 | 59.7277 |
| Sde-CMat_203_hap1 | Sde-CMat_203_hap2 | asm5 | 5 | 150000 | 187 | 1154337930 | 663310292 | 0.989066 | 0.0109338 | 59.6471 |
| Sde-CMat_203_hap1 | Sde-CPla_115_hap1 | asm5 | 5 | 150000 | 1224 | 1039243217 | 612232193 | 0.989284 | 0.010716 | 59.8676 |
| Sde-CMat_203_hap1 | Sde-CPla_115_hap2 | asm5 | 5 | 150000 | 1118 | 972432100 | 579473149 | 0.990374 | 0.00962603 | 59.78 |
| Sde-CMat_203_hap1 | Sde-CTlk_104_hap1 | asm5 | 5 | 150000 | 340 | 1094939406 | 600483107 | 0.987883 | 0.0121167 | 59.3147 |
| Sde-CMat_203_hap1 | Sde-CTlk_104_hap2 | asm5 | 5 | 150000 | 352 | 1115134072 | 612429244 | 0.987937 | 0.0120628 | 59.1193 |
| Sde-CMat_203_hap2 | Sde-CPla_115_hap1 | asm5 | 5 | 150000 | 1195 | 1037140339 | 608596444 | 0.989229 | 0.0107712 | 59.7967 |
| Sde-CMat_203_hap2 | Sde-CPla_115_hap2 | asm5 | 5 | 150000 | 1091 | 975177991 | 582487372 | 0.99042 | 0.0095801 | 59.9368 |
| Sde-CMat_203_hap2 | Sde-CTlk_104_hap1 | asm5 | 5 | 150000 | 331 | 1108179612 | 601316653 | 0.987576 | 0.012424 | 59.4018 |
| Sde-CMat_203_hap2 | Sde-CTlk_104_hap2 | asm5 | 5 | 150000 | 314 | 1120692261 | 611281189 | 0.987877 | 0.0121231 | 59.6242 |
| Sde-CPla_115_hap1 | Sde-CPla_115_hap2 | asm5 | 5 | 150000 | 1351 | 894668332 | 545832655 | 0.990605 | 0.00939451 | 59.7654 |
| Sde-CPla_115_hap1 | Sde-CTlk_104_hap1 | asm5 | 5 | 150000 | 1183 | 994763018 | 558417378 | 0.987719 | 0.0122809 | 59.8233 |
| Sde-CPla_115_hap1 | Sde-CTlk_104_hap2 | asm5 | 5 | 150000 | 1191 | 1007747117 | 568665998 | 0.988074 | 0.0119259 | 59.6809 |
| Sde-CPla_115_hap2 | Sde-CTlk_104_hap1 | asm5 | 5 | 150000 | 1093 | 933509232 | 535097406 | 0.989094 | 0.0109061 | 59.8866 |
| Sde-CPla_115_hap2 | Sde-CTlk_104_hap2 | asm5 | 5 | 150000 | 1103 | 957994909 | 551845438 | 0.989273 | 0.0107271 | 59.6954 |
| Sde-CTlk_104_hap1 | Sde-CTlk_104_hap2 | asm5 | 5 | 150000 | 376 | 1072689769 | 596217101 | 0.988398 | 0.011602 | 59.0931 |

</details>


## 9. Pangenome

Minigraph-Cactus pangenome graph for *Spratelloides delicatulus*, over 10 haplotypes.

### Downstream files

Graph products for downstream use, split by graph: **clip** is the reference-anchored graph (the default for mapping and variant work); **full** retains sequence that had no reference alignment. Paths are relative to this report; [`pangenome_manifest.tsv`](pangenome/Spratelloides_delicatulus/pangenome_manifest.tsv) carries the same set machine-readably.

| product | clip graph | full graph | description |
|---|---|---|---|
| Graph (GBZ) | [`Spratelloides_delicatulus.gbz`](https://crest-files.tamucc.edu/genomics/sde_genome_assembly/Spratelloides_delicatulus.gbz) | [`Spratelloides_delicatulus.full.gbz`](https://crest-files.tamucc.edu/genomics/sde_genome_assembly/Spratelloides_delicatulus.full.gbz) | Compressed graph + haplotype paths — the mapping substrate (`vg giraffe`, `vg call`). |
| Graph (GFA) | [`Spratelloides_delicatulus.gfa.gz`](https://crest-files.tamucc.edu/genomics/sde_genome_assembly/Spratelloides_delicatulus.gfa.gz) | [`Spratelloides_delicatulus.full.gfa.gz`](https://crest-files.tamucc.edu/genomics/sde_genome_assembly/Spratelloides_delicatulus.full.gfa.gz) | Text graph interchange (gzipped) — panacus, odgi, Bandage, `vg convert`. |
| Snarls | [`Spratelloides_delicatulus.snarls`](https://crest-files.tamucc.edu/genomics/sde_genome_assembly/Spratelloides_delicatulus.snarls) | [`Spratelloides_delicatulus.full.snarls`](https://crest-files.tamucc.edu/genomics/sde_genome_assembly/Spratelloides_delicatulus.full.snarls) | Snarl (bubble) decomposition — the genotyping units for `vg call`. |
| Haplotype index | [`Spratelloides_delicatulus.hapl`](https://crest-files.tamucc.edu/genomics/sde_genome_assembly/Spratelloides_delicatulus.hapl) | — | Haplotype-sampling index — builds personalized reference graphs for `vg giraffe`. |
| Variant catalog | [`Spratelloides_delicatulus.variants.vcf.gz`](https://crest-files.tamucc.edu/genomics/sde_genome_assembly/Spratelloides_delicatulus.variants.vcf.gz) | — | Variant catalog vs the reference path: top-level bubbles, one allele per row, SNP/indel/SV. |
| Variant catalog index | [`Spratelloides_delicatulus.variants.vcf.gz.tbi`](https://crest-files.tamucc.edu/genomics/sde_genome_assembly/Spratelloides_delicatulus.variants.vcf.gz.tbi) | — | Tabix index for the variant catalog. |
| Reference FASTA | [`Spratelloides_delicatulus.reference.fa`](https://crest-files.tamucc.edu/genomics/sde_genome_assembly/Spratelloides_delicatulus.reference.fa) | — | Reference-path FASTA — surjection target and linear-coordinate seam. |
| Reference FASTA index | [`Spratelloides_delicatulus.reference.fa.fai`](https://crest-files.tamucc.edu/genomics/sde_genome_assembly/Spratelloides_delicatulus.reference.fa.fai) | — | faidx index for the reference-path FASTA. |

*The whole-graph odgi graphs (`.og`, `.full.og`) are published alongside these but not linked above — they are far larger than every other product. Resolve them from the manifest by the `odgi` role, or rebuild either from its GFA with `odgi build -g`. Per-chromosome graphs sit under `<label>.chroms/`. Read-mapping indexes (`.dist` / `.min`) are not built here: haplotype sampling makes them sample-specific, so the consuming pipeline builds them from the GBZ + `.hapl` at map time.*

### Graph

| property | value |
|---|---:|
| Haplotypes | 10 |
| Total length | 1.927 Gb |
| Nodes | 151,700,713 |
| Edges | 210,093,543 |
| Paths | 9,019 |
| Acyclic | cyclic |

### Variant catalog

| class | count |
|---|---:|
| SNP | 23,122,105 |
| Indel | 16,741,186 |
| SV | 1,567,272 (INS 953,489 / DEL 613,783) |

Total: 41,430,563 variants relative to the reference path.

### Openness / growth

Pangenome 1,926.9 Mb · core 307.7 Mb · accessory 806.6 Mb · private 812.6 Mb. Heaps' γ = 0.32 (open).

![Growth and core curves](pangenome/Spratelloides_delicatulus/Spratelloides_delicatulus.growth_curves.png)

![Coverage histogram](pangenome/Spratelloides_delicatulus/Spratelloides_delicatulus.coverage_histogram.png)

### Structural variants

![SV size spectrum](pangenome/Spratelloides_delicatulus/Spratelloides_delicatulus.sv_size_histogram.png)

### Population structure

Ordination (PCoA) and neighbour-joining trees from graph-similarity distances (shared node content — SNPs, indels and SVs), at two levels: per haploid assembly (10 haplotypes, incl. the reference) and per diploid individual (haplotypes aggregated).

**Per haplotype**

![Haplotype PCoA](pangenome/Spratelloides_delicatulus/Spratelloides_delicatulus.pca_haplotype.png)

![Haplotype NJ tree](pangenome/Spratelloides_delicatulus/Spratelloides_delicatulus.njtree_haplotype.png)

**Per individual**

![Individual PCoA](pangenome/Spratelloides_delicatulus/Spratelloides_delicatulus.pca_individual.png)

![Individual NJ tree](pangenome/Spratelloides_delicatulus/Spratelloides_delicatulus.njtree_individual.png)

### Graph quality

| metric | value |
|---|---:|
| Re-alignment identity | 0.931367 (edit rate 0.068633) |
| Realigned | 10.503 Gb over 42,040 alignments |
| Mean node degree | 2.76984 (max 17) |
| Mean links length | 398.014 bp |

*Re-alignment identity is minigraph-GAF level (coarse); per-haplotype assembly BUSCO/QV remain the authoritative per-assembly completeness/accuracy metrics.*

## 10. Methods and Citations

This run processed 5 samples (HiFi (5), Hi-C (5)). Genome assembly used hifiasm. Haplotypic duplication was removed with purge_dups. Haploid genome size and heterozygosity were estimated from k-mer spectra (Jellyfish + GenomeScope 2). Assemblies were screened for contaminants with NCBI FCS (FCS-Adaptor and FCS-GX). Contigs were scaffolded against Hi-C data with YaHS. Organelle genomes were assembled with MitoHiFi. Assembly quality was assessed with BUSCO (per-sample lineage; see section 2), Merqury (consensus QV and k-mer completeness) and QUAST (contiguity), with read coverage from minimap2/SAMtools alignments; telomeric repeats were surveyed with tidk. Scaffold ends were extended into telomeric repeats with teloclip. Synteny was visualised from minimap2 alignments (gggenomes). Per-step parameters and exact software versions are recorded in the pipeline's Nextflow execution reports.

### Tool References

Primary references for the tools used in this run. Exact versions are recorded in the pipeline's Nextflow execution reports (not reproduced here).

- **gcl_genome_assembly** — TAMU-CC Genomics Core Lab. https://github.com/tamucc-gcl/gcl_genome_assembly
- **hifiasm** — Cheng et al. (2021) *Nat. Methods* 18:170-175.
- **purge_dups** — Guan et al. (2020) *Bioinformatics* 36:2896-2898.
- **NCBI FCS** (FCS-GX, FCS-Adaptor) — Astashyn et al. (2024) *Genome Biol.* 25:60.
- **YaHS** — Zhou et al. (2023) *Bioinformatics* 39:btac808.
- **MitoHiFi** — Uliano-Silva et al. (2023) *BMC Bioinformatics* 24:288.
- **Jellyfish** — Marçais & Kingsford (2011) *Bioinformatics* 27:764-770.
- **GenomeScope 2** — Ranallo-Benavidez et al. (2020) *Nat. Commun.* 11:1432.
- **BUSCO** — Simão et al. (2021) *Bioinformatics* 31:3210-3212.
- **Merqury** — Rhie et al. (2020) *Genome Biol.* 21:245.
- **QUAST** — Gurevich et al. (2013) *Bioinformatics* 29:1072-1075.
- **minimap2** — Li (2018) *Bioinformatics* 34:3094-3100.
- **SAMtools** — Danecek et al. (2021) *GigaScience* 10:giab008.
- **tidk** — Brown et al. (2021) *Bioinformatics* 41:btaf049.
- **MultiQC** — Ewels et al. (2016) *Bioinformatics* 32:3047-3048.
- **teloclip** — Taranto. https://github.com/Adamtaranto/teloclip
- **gggenomes** — Hackl et al. (2024) *arXiv* arXiv:2411.13556
- **Minigraph-Cactus** — Hickey et al. (2024) *Nat. Biotechnol.* 42:663-673.
- **minigraph** — Li et al. (2020) *Genome Biol.* 21:265.
- **vg** (variation-graph toolkit) — Garrison et al. (2018) *Nat. Biotechnol.* 36:875-879.
- **odgi** — Guarracino et al. (2022) *Bioinformatics* 38:3319-3326.
- **Panacus** — Parmigiani et al. (2024) *Bioinformatics* 40:btae720.
- **vcflib/vcfbub** — Garrison et al. (2022) *PLoS Comput. Biol.* 18:e1009123.
- **BCFtools** — Danecek et al. (2021) *GigaScience* 10:giab008.
- **ape** (neighbour-joining) — Paradis & Schliep (2019) *Bioinformatics* 35:526-528.

### Software Versions

| Tool | Version |
| --- | --- |
| 0.2.65 | NA |
| branch: | NA |
| build_env:rustc 1.97.0 (2d8144b78 2026-07-07), | NA |
| build_time:2026-07-20 15:27:13 +00:00 | NA |
| BUSCO | 5.3.2 |
| commit_hash: | NA |
| FastQC | v0.12.1 |
| FCS-GX | v0.5.5 |
| GenomeScope | 2.0 |
| HARMONIZE_SCAFFOLDS:HARMONIZE_CANDIDATES | python	3.10.20 |
| HARMONIZE_SCAFFOLDS:HARMONIZE_SCORE | minimap2	2.31-r1302 |
| HARMONIZE_SCAFFOLDS:HARMONIZE_SELECT | python	3.10.20 |
| HARMONIZE_SCAFFOLDS:HARMONIZE_SPECIES | minimap2	2.31-r1302 |
| hifiasm | 0.25.0-r726 |
| Inspector | 1.3.1 |
| Jellyfish | 2.3.1 |
| meryl | branch HEAD +0 changes (r977 5f3f1eb8c2392656d591a47a872c2efd3706e3ac) |
| minimap2 | 2.30-r1287 |
| MitoHiFi | Matplotlib created a temporary config/cache directory at /tmp/matplotlib-522u0f_k because the default path (/home/jselwyn/.config/matplotlib) is not a writable directory; it is highly recommended to set the MPLCONFIGDIR environment variable to a writable directory, in particular to speed up the import of Matplotlib and to better support multiprocessing. |
| PANGENOME:CACTUS_PANGENOME | cactus	9.1.2 |
| PANGENOME:MULTIQC_PANGENOME | multiqc	1.31 |
| PANGENOME:PANGENOME_2D_VIZ | odgi	v0.9.2-0-gbe6a0202 |
| PANGENOME:PANGENOME_GROWTH | panacus	0.5.2 |
| PANGENOME:PANGENOME_ODGI_STATS_MQC | odgi	v0.9.2-0-gbe6a0202 |
| PANGENOME:PANGENOME_PCA_NJ | odgi	v0.9.2-0-gbe6a0202 |
| PANGENOME:PANGENOME_QC | odgi	v0.9.2-0-gbe6a0202 |
| PANGENOME:PANGENOME_REPORT | Rscript	(2024-02-29) |
| PANGENOME:PANGENOME_VARIANTS | bcftools	1.24 |
| process | tool	version |
| purge_dups | 1.2.6 |
| QUAST | 5.3.0 |
| samtools | 1.23 |
| teloclip | 0.3.4 |
| TGSGapCloser | 1.2.1 |
| tidk | Warning! No clades found in the database. Run 'build' to fetch the latest data. |
| YaHS | 1.2.2 |

---

*Report generated on 2026-08-19 09:40:04 by the Genome Assembly Pipeline.*
