/*
========================================================================================
    PANGENOME VARIANT CATALOG MODULE
========================================================================================
    Repo location: modules/pangenome_variants.nf

    Turn the raw decomposed VCF (vg deconstruct output, with LV/PS/AT tags) into a report
    variant catalog: keep top-level, position-anchored sites (vcfbub --max-level 0 + a
    reference-allele-length cap), split multiallelics, then classify every allele relative
    to the reference path by REF/ALT length into SNP / indel / SV, with an SV size spectrum
    and an INS/DEL/COMPLEX type breakdown.

    INV/DUP typing is NOT attempted here -- deconstruct emits sequence-resolved alleles, so
    inversions/duplications appear as COMPLEX substitutions; proper typing needs realignment
    (a later refinement; minigraph's own inversion count is a separate signal).

    Input : tuple(taxid, raw_vcf)          (CACTUS_PANGENOME.out.raw_vcf)
    Output: filtered vcf (+tbi) / variant_summary / sv_sizes / versions
========================================================================================
*/

process PANGENOME_VARIANTS {
    tag "${taxid}"
    label 'pangenome_variants'

    publishDir "${params.outdir}/pangenome/${taxid}", mode: params.publish_dir_mode

    input:
    tuple val(taxid), path(raw_vcf)

    output:
    tuple val(taxid), path("${taxid}.variants.vcf.gz"),     emit: vcf
    tuple val(taxid), path("${taxid}.variants.vcf.gz.tbi"), emit: tbi
    tuple val(taxid), path("${taxid}.variant_summary.tsv"), emit: summary
    tuple val(taxid), path("${taxid}.sv_sizes.tsv"),        emit: sv_sizes
    path("versions.tsv"),                                   emit: versions

    script:
    def minsv = params.pangenome_sv_min_bp ?: 50
    def maxref = params.pangenome_vcfbub_max_ref ?: 100000
    """
    set -euo pipefail

    # top-level, position-anchored sites (drop nested + oversized bubbles), one allele per row
    vcfbub --input ${raw_vcf} --max-level 0 --max-ref-length ${maxref} \\
        | bcftools norm -m -any \\
        | bgzip > ${taxid}.variants.vcf.gz
    tabix -p vcf ${taxid}.variants.vcf.gz

    bcftools query -f '%REF\\t%ALT\\n' ${taxid}.variants.vcf.gz > refalt.tsv

    # counts: SNP (|REF|=|ALT|=1), INDEL (size < minsv), SV (size >= minsv)
    awk -v minsv=${minsv} 'BEGIN{FS=OFS="\\t"}
        { rl=length(\$1); al=length(\$2); d=(rl>al?rl-al:al-rl)
          if(rl==1 && al==1)      snp++
          else if(d<minsv)        indel++
          else { sv++; if(al>rl) ins++; else if(rl>al) del++; else cpx++ } }
        END{ print "class","count"
             print "SNP",snp+0; print "INDEL",indel+0; print "SV",sv+0
             print "SV_INS",ins+0; print "SV_DEL",del+0; print "SV_COMPLEX",cpx+0 }' \\
        refalt.tsv > ${taxid}.variant_summary.tsv

    # SV size spectrum (bp) + type, for the histogram
    awk -v minsv=${minsv} 'BEGIN{FS=OFS="\\t"; print "sv_size_bp","sv_type"}
        { rl=length(\$1); al=length(\$2); d=(rl>al?rl-al:al-rl)
          if(!(rl==1 && al==1) && d>=minsv){ t=(al>rl?"INS":(rl>al?"DEL":"COMPLEX")); print d,t } }' \\
        refalt.tsv > ${taxid}.sv_sizes.tsv

    {
      printf 'process\\ttool\\tversion\\n'
      printf '%s\\tbcftools\\t%s\\n' "${task.process}" "\$(bcftools --version 2>&1 | awk 'NR==1{print \$2}')"
      printf '%s\\tvcfbub\\t%s\\n'   "${task.process}" "\$(vcfbub --help 2>&1 | awk 'NR==1{print \$NF}')"
    } > versions.tsv
    """

    stub:
    """
    : > ${taxid}.variants.vcf.gz
    : > ${taxid}.variants.vcf.gz.tbi
    printf 'class\\tcount\\n' > ${taxid}.variant_summary.tsv
    printf 'sv_size_bp\\tsv_type\\n' > ${taxid}.sv_sizes.tsv
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
