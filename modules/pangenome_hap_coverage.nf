/*
========================================================================================
    PER-HAPLOTYPE COVERAGE DECOMPOSITION MODULE
========================================================================================
    Repo location: modules/pangenome_hap_coverage.nf

    panacus `hist` gives the MARGINAL coverage histogram h(i) = bp covered by exactly i
    haplotypes. It cannot say which haplotype owns the private (i == 1) sequence. This
    process computes the JOINT table

        M[haplotype, i] = bp this haplotype traverses that is covered by exactly i haplotypes

    directly from the clip GFA, via py_scripts/gfa_hap_coverage.py. Grouping is by the first
    two PanSN fields (sample, haplotype), i.e. identical to `panacus --groupby-haplotype`, so
    the marginal sum_h M[h,i]/i reproduces the panacus histogram exactly -- the emitted
    hap_coverage_check.tsv exists to be diffed against it.

    Single-threaded and I/O bound: three streaming passes over the GFA -- segment lengths +
    haplotype group scan, then one bitmask word per node, then per-contig attribution of the
    private column. Memory scales with NODE COUNT (~25 B/node), not genome size.

    Input : tuple(taxid, gfa)   (CACTUS_PANGENOME.out.gfa -- the clip GFA .gfa.gz)
            hapcov_script       (py_scripts/gfa_hap_coverage.py)
    Output: matrix / hap_private / by_contig / check / versions
========================================================================================
*/

process PANGENOME_HAP_COVERAGE {
    tag "${taxid}"
    label 'pangenome_hap_coverage'

    publishDir "${params.outdir}/pangenome/${taxid}", mode: params.publish_dir_mode

    input:
    tuple val(taxid), path(gfa)
    path(hapcov_script)

    output:
    tuple val(taxid), path("${taxid}.hap_coverage_matrix.tsv"), emit: matrix
    tuple val(taxid), path("${taxid}.hap_private.tsv"),         emit: hap_private
    tuple val(taxid), path("${taxid}.hap_coverage_check.tsv"),  emit: check
    tuple val(taxid), path("${taxid}.hap_private_by_contig.tsv"), emit: by_contig
    // ---- private SEGMENTS (pass 5) ----------------------------------------------------
    // Maximal runs of consecutive private (coverage == 1) nodes along each haplotype walk.
    // The four outputs above are TOTALS and cannot answer "what does the private-haplotype
    // size distribution look like"; these can. optional: true so the process still succeeds
    // if py_scripts/gfa_hap_coverage.py has not yet had its fifth pass added.
    tuple val(taxid), path("${taxid}.private_segments.tsv"),         emit: segments,  optional: true
    tuple val(taxid), path("${taxid}.private_segments.bed"),         emit: segments_bed, optional: true
    tuple val(taxid), path("${taxid}.private_segment_spectrum.tsv"), emit: spectrum,  optional: true
    path("versions.tsv"),                                       emit: versions

    script:
    // Floor for the per-segment TABLE and BED only. The spectrum and the totals are written
    // unfiltered regardless, so raising this bounds row count without hiding sequence --
    // which matters because 33.5M of the 33.7M segments are under 1 kb.
    def minpriv = params.pangenome_private_min_bp ?: 1000
    """
    set -euo pipefail

    # the script sniffs gzip magic bytes, so the .gfa.gz is read in place (no scratch copy)
    python3 ${hapcov_script} \\
        --gfa ${gfa} \\
        --label ${taxid} \\
        --min-private-bp ${minpriv} \\
        --outdir .

    {
      printf 'process\\ttool\\tversion\\n'
      printf '%s\\tpython\\t%s\\n' "${task.process}" "\$(python3 -c 'import sys; print("%d.%d.%d" % sys.version_info[:3])')"
      printf '%s\\tnumpy\\t%s\\n'  "${task.process}" "\$(python3 -c 'import numpy; print(numpy.__version__)')"
    } > versions.tsv
    """

    stub:
    """
    printf 'haplotype\\tcoverage_level\\tbp\\n' > ${taxid}.hap_coverage_matrix.tsv
    printf 'haplotype\\tsample\\thap\\tprivate_bp\\thap_bp\\tpct_of_private\\tpct_of_haplotype\\n' > ${taxid}.hap_private.tsv
    printf 'coverage_level\\tsum_over_haplotypes_bp\\tbp_reconstructed\\n' > ${taxid}.hap_coverage_check.tsv
    printf 'haplotype\\tcontig\\tprivate_bp\\tcontig_graph_bp\\tpct_of_contig\\tpct_of_hap_private\\tpct_of_pangenome_private\\n' > ${taxid}.hap_private_by_contig.tsv
    printf 'haplotype\\tcontig\\tseg_index\\tn_nodes\\tseg_bp\\tstart\\tend\\n' > ${taxid}.private_segments.tsv
    : > ${taxid}.private_segments.bed
    printf 'scope\\tsize_bin\\tn_segments\\tsegment_bp\\n' > ${taxid}.private_segment_spectrum.tsv
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
