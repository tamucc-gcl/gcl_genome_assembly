/*
========================================================================================
    PANGENOME CLASSIFY MODULE
========================================================================================
    Repo location: modules/pangenome_classify.nf

    Replaces the length-based awk classifier that used to live inside PANGENOME_VARIANTS
    with a topological one derived from the graph's own allele traversals, plus an
    allele-frequency layer read from AC/AN.

    WHY THE OLD CLASSIFIER HAD TO GO
    --------------------------------
    It ran `bcftools query -f '%REF\\t%ALT\\n'` and classified on string lengths alone,
    which threw away AT (the signed allele traversal), LV/PS (nesting), AC/AN and every
    sample column one line before deciding. `SV_COMPLEX` and `SV_BLOCKSUB` were not variant
    classes at all -- they were the buckets the heuristic used when it could not decide.

    Topology reclassifies real content: on Spratelloides 233,015 length-called SV_INS and
    178,583 SV_DEL are actually SUBST, because they lose reference nodes as well as gaining
    novel ones, which length arithmetic cannot see.

    THREE INVARIANTS ENFORCED IN THE SCRIPT, NOT ASSUMED FROM THE WIRING
    -------------------------------------------------------------------
    1. AT IS ONLY VALID PRE-DECOMPOSITION. vcfwave and `bcftools norm -m` rewrite REF/ALT
       while AT is inherited from the parent, so allele i stops matching traversal i+1. The
       script detects decomposition markers in the header and refuses topology; with
       --tier parent it exits non-zero rather than producing confident nonsense.
    2. DUP is restricted to nodes the REFERENCE also visits. Testing
       alt_count > ref_count.get(n, 0) fires for any novel node, which made DUP a synonym
       for NOVEL_INS -- both returned an identical 493,623 during development.
    3. PER-ALLELE bp IS NOT A TOTAL. 41.9M records carry 55.3M alt alleles and max(REF,ALT)
       counts ALT length when ALT is longer, so the SV sum is 2.95 Gb against a ~1 Gb
       reference -- 28x inflated for SUBST. The summary therefore carries THREE measures:
       per-allele bp (size spectrum input only), merged reference footprint
       (reference-anchored, and reference-BIASED: INS scores 5.2 Mb against 247 Mb of novel
       sequence), and pangenome node footprint / novel node bp, which are graph-native and
       bounded by graph length rather than by one haplotype's length.

    --gfa is what enables the node measures, and its own S-line total is checked against
    odgi: 1,926,892,905 bp for the clip graph, matching `odgi stats -S` exactly.

    NOT SCATTERED BY CHROMOSOME. The whole 23,877,797-record parent view classifies in
    well under an hour single-threaded, so a tabix split would add 15x the task count, a
    merge step and a groupTuple barrier for no measurable gain. Same reasoning that keeps
    PANGENOME_HAP_COVERAGE whole-graph. If a future cohort makes this the bottleneck, adding
    the scatter is a channel change, not a module change.

    Input : tuple(taxid, flavor, tier, vcf), gfa, script
    Output: summary / sv_sizes / bed / spectrum / xtab / cooccurrence / af / private / audit
========================================================================================
*/

process PANGENOME_CLASSIFY {
    tag "${taxid}:${flavor}:${tier}"
    label 'pangenome_classify'

    publishDir "${params.outdir}/pangenome/${taxid}/variants", mode: params.publish_dir_mode

    input:
    tuple val(taxid), val(flavor), val(tier), path(vcf)
    tuple val(gtaxid), path(gfa)
    path(script)

    output:
    tuple val(taxid), val(flavor), val(tier),
          path("*.variant_summary.tsv"),       emit: summary
    tuple val(taxid), val(flavor), val(tier),
          path("*.sv_sizes.tsv"),              emit: sv_sizes
    tuple val(taxid), val(flavor), val(tier),
          path("*.variants.bed"),              emit: bed
    tuple val(taxid), val(flavor), val(tier),
          path("*.size_spectrum.tsv"),         emit: spectrum
    tuple val(taxid), val(flavor), val(tier),
          path("*.af_spectrum.tsv"),           emit: af
    tuple val(taxid), val(flavor), val(tier),
          path("*.representation_audit.tsv"),  emit: audit
    path("*.topology_xtab.tsv"),               emit: xtab,        optional: true
    path("*.label_cooccurrence.tsv"),          emit: cooccurrence, optional: true
    path("*.private_variants.tsv"),            emit: private_variants, optional: true
    path("versions.tsv"),                      emit: versions

    script:
    def minsv = params.pangenome_sv_min_bp ?: 50
    def bins  = params.pangenome_sv_bins ?: '50,100,500,1000,5000,10000,50000,100000,500000,1000000'
    def usegfa = gfa && gfa.name != 'NO_FILE' ? "--gfa ${gfa}" : ''
    """
    set -euo pipefail

    python3 ${script} \\
        --vcf ${vcf} \\
        --outdir . \\
        --label ${taxid} \\
        --flavor ${flavor} \\
        --tier ${tier} \\
        --chrom all \\
        --minsv ${minsv} \\
        --sv-bins '${bins}' \\
        ${usegfa}

    {
      printf 'process\\ttool\\tversion\\n'
      printf '%s\\tpython\\t%s\\n' "${task.process}" "\$(python3 --version 2>&1 | awk '{print \$2}')"
    } > versions.tsv
    """

    stub:
    """
    S=${taxid}.${flavor}.${tier}.all
    for f in variant_summary sv_sizes size_spectrum af_spectrum representation_audit \\
             topology_xtab label_cooccurrence private_variants; do : > \$S.\$f.tsv; done
    : > \$S.variants.bed
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
