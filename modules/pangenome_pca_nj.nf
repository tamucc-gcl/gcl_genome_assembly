/*
========================================================================================
    PANGENOME HAPLOTYPE SIMILARITY MODULE  (workstream D)
========================================================================================
    Repo location: modules/pangenome_pca_nj.nf

    Per-HAPLOTYPE distance matrix from the graph via `odgi similarity`, grouped by PanSN
    haplotype (sample#hap) so the graph's many contig-paths collapse to one unit per
    haploid genome -- INCLUDING the reference (it is just another path). This is the correct
    granularity for a pangenome ordination/tree: each haplotype is one point / one leaf.

    Distances come from shared graph content (node-presence), so the tree/ordination reflect
    overall sequence relatedness (SNPs + indels + SVs), not SNPs alone. Downstream is
    PANGENOME_POPSTRUCT (PCoA scatter + neighbour-joining tree).

    With >= 3 haplotypes this is non-degenerate (so the n=4 test plots for real). Uses odgi
    from the cactus image -- no plink.

    Input : tuple(taxid, og)     (CACTUS_PANGENOME.out.og -- the clip whole-genome graph)
    Output: similarity (long-format pairwise distances) / versions
========================================================================================
*/

process PANGENOME_PCA_NJ {
    tag "taxid_${taxid}"
    label 'cactus_tools'

    publishDir "${params.outdir}/pangenome/${taxid}/popstruct", mode: params.publish_dir_mode

    input:
    tuple val(taxid), path(og)

    output:
    tuple val(taxid), path("${taxid}.similarity.tsv"), emit: similarity
    path("versions.tsv"),                              emit: versions

    script:
    // group by PanSN haplotype: part of the path name before the 2nd '#' (sample#hap).
    // VERIFY on your odgi build with `odgi similarity --help` (delimiter-position flag name;
    // some builds also expose a --group-by-haplotype shortcut). Group names should be the 4
    // haplotypes, e.g. Sde-CBau_104#1, Sde-CBau_104#2, Sde-CMat_203_dip_hap1#0, ..._hap2#0.
    """
    set -euo pipefail
    export HOME="\$PWD"

    odgi similarity -i ${og} -D '#' -p 2 -d -t ${task.cpus} > ${taxid}.similarity.tsv

    {
      printf 'process\\ttool\\tversion\\n'
      printf '%s\\todgi\\t%s\\n' "${task.process}" "\$(odgi version 2>&1 | head -n1)"
    } > versions.tsv
    """

    stub:
    """
    printf 'group.a\\tgroup.b\\tjaccard\\teuclidean\\n' > ${taxid}.similarity.tsv
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
