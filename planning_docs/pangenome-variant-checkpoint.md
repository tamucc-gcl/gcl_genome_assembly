# GREF catalog checkpoint and structural-variant roadmap

## First change ready for cluster validation

PANGENOME_VARIANT_AUDIT reads the exact `<taxid>.gref.vcf.gz` and its matching TBI, gated by the existing construction contract. The broad legacy gref output glob can include raw and standard VCFs, so it is deliberately not passed directly to consumers. Cactus and all accepted sharing/comparison computations are unchanged.

bcftools 1.21 supplies the sample list, index counts by coordinate sequence and full-file statistics. A small adapter verifies sample identity against the ledger and index record totals against the full scan. Published outputs under pangenome/<taxid>/variants include the original bcftools statistics, index coordinate counts, sample identity mapping, record summary, JSON audit and collapsed Markdown report. The VCF is neither rewritten nor duplicated here.

These are representation-level record statistics. Multiallelic records can contribute to several type categories; categories are not additive event counts. No population frequency or genotype-quality interpretation is enabled. GREF coordinate sequences cannot all be treated as ordinary assembly chromosomes. The audit does not establish REF/FASTA concordance, correct per-site ploidy, or traversal-to-allele alignment.

The six columns in the current export represent five biological individuals; this is a review expectation, not a sample-specific pipeline rule. The ledger includes ten haplotypes. Future compatible cohorts use their own ledger counts.

## Legacy inventory and disposition

- modules/pangenome_variants.nf: raw VCF -> vcfbub blocks -> default vcfwave -> normalization, with a separate parent view. Keep dormant. Mandatory fine decomposition conflicts with the agreed on-demand policy. Before enabling a consumer, choose its required representation and preserve provenance.
- modules/pangenome_classify.nf and its classifier: custom traversal-based categories and several bp summaries. Do not automatically reactivate. Inherited traversal annotations must not be interpreted as allele-specific after REF/ALT decomposition unless correspondence has been verified.
- modules/pangenome_rearrange.nf and rearrangement plotting: untangle-based candidates, orientation checks, inversion and duplication summaries. Preserve while assessing published replacements, but do not describe candidate flags as validated events.
- Legacy variant figures: retain useful record/allele counts and size spectra when their units are explicit. Reference-footprint summaries require correctly separated coordinate sequences; graph-derived coordinates must not be pooled into chromosome reference footprints. All replacement plots use tidyverse, ggplot2 and patchwork.
- Remove dataset-specific numerical assertions from active descriptions when migrating code. A block substitution does not by itself prove unrelated sequence; a representation change does not establish increased biological accuracy.

## Next bounded experiments

1. Accept the standard export inventory, identity mapping and record counts. Review header metadata and actual allele/genotype examples before deciding normalization/decomposition.
2. Define tiny known-answer cases: SNP/indel, nonreference insertion with internal variation, inversion, tandem duplication, inverted duplication, translocation and ambiguous repeat. Include fragmented-query cases and an explicit cannot-assess outcome.
3. Compare a limited published-tool shortlist against those cases and selected real loci. The agreed plan identifies SyRI for chromosome-scale assembly rearrangements and other graph/assembly tools for evaluation. Confirm pinned-version input requirements before choosing a benchmark command. No tool has been selected by this checkpoint.
4. Assess source-to-output event/allele correspondence, nonreference coordinates, false positives and total alignment-plus-analysis resources. Avoid adding overlapping callsets as if their counts were independent.
5. Choose a supported replacement or explicitly bounded retained custom component; then implement the standard structural-variant summaries and plots. INV/DUP/translocation remain required capabilities, not silently dropped scope.

The full original plan and candidate list remain in assembly-pangenome-agreed-plan.md, section 7. No full-cohort structural-variant benchmark is started in this checkpoint.

Toolmaker reference for the statistics interface: https://samtools.github.io/bcftools/bcftools.html#stats . The published raw statistics retain their original definitions; the Markdown summary does not relabel these as biological event counts.

## Run

After syncing, from the cluster project root:

```bash
export PYTHONDONTWRITEBYTECODE=1
python3 -m unittest discover -s gcl_genome_assembly/tests -p 'test_*.py'
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch \
  data/assembly_samplesheet_original.csv data/hic_readsets.none.csv validate
# After validation succeeds:
sbatch gcl_genome_assembly/run_assembly_checkpoint.sbatch \
  data/assembly_samplesheet_original.csv data/hic_readsets.none.csv graph_build
```

Expected new work: catalog audit plus the routine summary. Keep results/work; Cactus and the accepted analyses should remain cached. No pipeline commands, tests or analyses were run locally. Static diff checks passed. The audit requests 2 CPUs/8 GB/8h; measured resource needs remain to be checked.

## Return

```bash
mkdir -p comparisons
stamp=$(date +%Y%m%d-%H%M%S)
find genome_assembly/pangenome -type f \( -path '*/variants/*' \
  -o -name pangenome_identity.tsv -o -name pangenome_manifest.tsv \) \
  -print0 > "comparisons/pangenome-variants-${stamp}.files"
tar --null -czf "comparisons/pangenome-variants-${stamp}.tar.gz" \
  -T "comparisons/pangenome-variants-${stamp}.files"
ls -ltrh comparisons/pangenome-variants-*.tar.gz logs/nextflow_graph_build_*.log
```

Return the archive and corresponding Nextflow/Slurm logs. Do not include the multi-GB VCF or graphs. If a new check fails, return its task log and small tables rather than rebuilding Cactus.
