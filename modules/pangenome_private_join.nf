/*
========================================================================================
    PANGENOME PRIVATE JOIN MODULE
========================================================================================
    Repo location: modules/pangenome_private_join.nf

    Collects every haplotype's PRIVATE_MAP and PRIVATE_KMER table for one graph flavour into
    a single per-segment table, and cross-tabulates the two verdicts.

    WHY THE CROSS-TABULATION IS THE DELIVERABLE
    -------------------------------------------
    PRIVATE_MAP and PRIVATE_KMER answer "is this really private?" by entirely independent
    routes: one aligns the segment against the other assemblies, the other counts its k-mers
    in the sample's own reads. Neither sees the other's evidence, so the four combinations
    each mean something specific:

      NOT_PRIVATE       + REPEAT_LIKE   present elsewhere AND high copy -> graph collapse
      PRIVATE_CONFIRMED + UNIQUE_LIKE   absent elsewhere AND single copy -> novel sequence
      NOT_PRIVATE       + UNIQUE_LIKE   present elsewhere but single copy -> the graph failed
                                        to merge homologous sequence; an alignment failure,
                                        not a repeat problem
      PRIVATE_CONFIRMED + REPEAT_LIKE   absent elsewhere but high copy WITHIN this sample ->
                                        haplotype-specific expansion

    The off-diagonals are the informative cases and they are invisible unless the two streams
    are tabulated together, which is the whole reason this process exists rather than
    publishing two independent tables.

    Broken out by chromosome, size bin and haplotype, because the live question is
    per-chromosome: chr8 and chr9 carry 23.6% and 22.9% private sequence against an ~11% floor
    across the other chromosomes, consistently across all ten haplotypes and both graph
    flavours. Consistency across independently assembled haplotypes argues for biology; this
    table is what decides it.

    ONE TABLE, BECAUSE THE TESTS ARE DEFERRED
    -----------------------------------------
    Statistical testing of these relationships is deliberately later work. The requirement is
    that when it happens it is a regression over columns that already exist rather than a
    fresh data-gathering exercise -- so segment size, alignment evidence and copy number all
    land in one row keyed on (haplotype, segment).

    Size bin edges default to the same values as the SV and private-segment spectra so all
    three are readable on one axis.

    Input : tuple(taxid, flavor, map_tables, kmer_tables), script
    Output: joined table / cross-tab / audit / versions
========================================================================================
*/

process PANGENOME_PRIVATE_JOIN {
    tag "${taxid}:${flavor}"
    label 'pangenome_private_join'

    publishDir "${params.outdir}/pangenome/${taxid}/private", mode: params.publish_dir_mode

    input:
    tuple val(taxid), val(flavor), path(map_tables), path(kmer_tables)
    path(script)

    output:
    tuple val(taxid), val(flavor),
          path("${taxid}.${flavor}.private_evidence.tsv"),       emit: evidence
    tuple val(taxid), val(flavor),
          path("${taxid}.${flavor}.private_evidence_xtab.tsv"),  emit: xtab
    tuple val(taxid), val(flavor),
          path("${taxid}.${flavor}.private_evidence_audit.tsv"), emit: audit
    path("versions.tsv"),                                        emit: versions

    script:
    def bins = params.pangenome_sv_bins ?: '50,100,500,1000,5000,10000,50000,100000,500000,1000000'
    """
    set -euo pipefail

    python3 ${script} \\
        --map ${map_tables} \\
        --kmer ${kmer_tables} \\
        --label ${taxid}.${flavor} \\
        --outdir . \\
        --bins '${bins}'

    # Both streams run on the SAME per-haplotype FASTAs, so every segment must appear in both.
    # A segment in one and not the other means a task failed or a haplotype was dropped from
    # one scatter but not the other -- which would quietly bias the cross-tabulation, since
    # the missing rows are not random with respect to verdict.
    A=${taxid}.${flavor}.private_evidence_audit.tsv
    mk=\$(awk -F'\\t' '\$1=="segments_without_kmer_row"{print \$2}' "\$A")
    mm=\$(awk -F'\\t' '\$1=="segments_without_map_row"{print \$2}'  "\$A")
    if [ "\${mk:-0}" -ne 0 ] || [ "\${mm:-0}" -ne 0 ]; then
        echo "[PRIVATE_JOIN ${taxid}:${flavor}] ERROR: \$mk segments lack k-mer data and" >&2
        echo "  \$mm lack mapping data. The two scatters processed different haplotype sets;" >&2
        echo "  the cross-tabulation would be biased. Check for a failed MAP or KMER task." >&2
        exit 1
    fi

    echo "[PRIVATE_JOIN ${taxid}:${flavor}] \$(awk -F'\\t' '\$1=="joined_segments"{print \$2}' "\$A") segments joined" >&2

    {
      printf 'process\\ttool\\tversion\\n'
      printf '%s\\tpython\\t%s\\n' "${task.process}" "\$(python3 --version 2>&1 | awk '{print \$2}')"
    } > versions.tsv
    """

    stub:
    """
    S=${taxid}.${flavor}
    printf 'haplotype\\tsample\\tsegment\\tcontig\\tstart\\tend\\tsegment_bp\\tn_other_assemblies\\tbest_identity\\taligned_frac_merged\\tmap_verdict\\tn_kmers\\tmean_copy\\tmedian_copy\\tfrac_absent\\tkmer_verdict\\tcombined\\n' > \$S.private_evidence.tsv
    printf 'scope\\tkey\\tcombined\\tn_segments\\tsegment_bp\\n' > \$S.private_evidence_xtab.tsv
    printf 'metric\\tvalue\\nsegments_without_kmer_row\\t0\\nsegments_without_map_row\\t0\\njoined_segments\\t0\\n' > \$S.private_evidence_audit.tsv
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
