/*
========================================================================================
    PANGENOME INVERSION RESCUE MODULE
========================================================================================
    Repo location: modules/pangenome_inversion_rescue.nf

    Recovers ALIGNMENT-RESCUED inversions from the SUBST bucket, and annotates the rest of
    SUBST by how much homology its alleles actually share.

    WHY IT IS A SEPARATE PROCESS FROM PANGENOME_CLASSIFY
    ---------------------------------------------------
    Classification is minutes; this is minimap2 over ~310,000 allele pairs. Splitting them
    is the same argument as splitting PANGENOME_VARIANTS from PANGENOME_CLASSIFY: the
    thresholds here will be tuned repeatedly and re-running the cheap step should not cost
    the expensive one. Both read the parent VCF independently and emit tables keyed on
    (chrom, pos, allele_idx), so there is no ordering constraint between them.

    WHAT IT FINDS, AND THE POSITIVE CONTROL THAT VALIDATES IT
    --------------------------------------------------------
    Romain et al. (bioRxiv 2025.03.14.643331) distinguish path-explicit inversions, visible
    to topology, from alignment-rescued ones whose alleles are disjoint paths and therefore
    land in SUBST. On Spratelloides: 411 rescued alleles, 2.0 Mb of minus-strand homology,
    every one under 100 kb, and NONE of the 288 SUBST alleles >=500 kb is an inversion.

    The control is what makes that trustworthy. --control routes the ~32 alleles that
    topology already calls path-explicit inversions through the same alignment test;
    30 of 32 were recovered. They are reported but excluded from the catalog, so there is no
    double counting. A low recovery rate would mean --min-frac is too strict and the SUBST
    calls are an underestimate -- the number is emitted every run for exactly that reason.

    THREE FAILURE MODES THE SCRIPT NOW GUARDS, ALL FOUND THE HARD WAY
    ----------------------------------------------------------------
    * minimap2 searches BOTH strands, so aligning revcomp(ALT) against REF is a no-op -- it
      finds the same homology and reports it on the opposite strand. Treating any alignment
      as evidence of inversion called 66% of SUBST and produced 1.74 Gb of "inverted"
      sequence in a 1 Gb genome. PAF column 5 is the signal: coverage is merged separately
      per strand and a call needs minus-strand dominance, --min-strand-margin.
    * Coverage must be read in the same frame as its denominator. Merging query intervals
      and dividing by min(ref_len, alt_len) produced fractions up to 2.8.
    * Summing PAF block lengths double-counts overlapping alignments. Intervals are merged.

    The k-mer prefilter is OFF by default. At this scale minimap2 on everything is
    affordable and cannot introduce a selection bias; the prefilter's sampling once passed
    only 0.08% of candidates, almost all of them large.

    Input : tuple(taxid, flavor, vcf), script
    Output: candidates / calls BED / summary / unresolved profile / versions
========================================================================================
*/

process PANGENOME_INVERSION_RESCUE {
    tag "${taxid}:${flavor}:${chrom}"
    label 'pangenome_inversion_rescue'

    publishDir "${params.outdir}/pangenome/${taxid}/variants", mode: params.publish_dir_mode

    input:
    tuple val(taxid), val(flavor), val(chrom), path(vcf)
    path(script)

    output:
    tuple val(taxid), val(flavor), path("*.rescue_candidates.tsv"),         emit: candidates
    tuple val(taxid), val(flavor), path("*.rescue_calls.bed"),              emit: calls
    tuple val(taxid), val(flavor), path("*.rescue_summary.tsv"),            emit: summary
    tuple val(taxid), val(flavor), path("*.rescue_unresolved_profile.tsv"), emit: profile
    path("versions.tsv"),                                                  emit: versions

    script:
    def minbp  = params.pangenome_inv_rescue_min_bp   ?: 1000
    def minfrac = params.pangenome_inv_rescue_min_frac ?: 0.8
    def margin = params.pangenome_inv_rescue_strand_margin ?: 2.0
    def minsv  = params.pangenome_sv_min_bp ?: 50
    def preset = params.pangenome_inv_rescue_preset ?: 'asm20'
    """
    set -euo pipefail

    python3 ${script} \\
        --vcf ${vcf} \\
        --outdir . \\
        --label ${taxid}.${flavor}.${chrom} \\
        --min-bp ${minbp} \\
        --minsv ${minsv} \\
        --min-frac ${minfrac} \\
        --min-strand-margin ${margin} \\
        --preset ${preset} \\
        --threads ${task.cpus} \\
        --control

    {
      printf 'process\\ttool\\tversion\\n'
      printf '%s\\tpython\\t%s\\n' "${task.process}" "\$(python3 --version 2>&1 | awk '{print \$2}')"
      printf '%s\\tminimap2\\t%s\\n' "${task.process}" "\$(minimap2 --version 2>&1 | head -n1)"
    } > versions.tsv
    """

    stub:
    """
    S=${taxid}.${flavor}.${chrom}
    for f in rescue_candidates rescue_summary rescue_unresolved_profile; do : > \$S.\$f.tsv; done
    : > \$S.rescue_calls.bed
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
