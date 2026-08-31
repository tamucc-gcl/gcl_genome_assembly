/*
========================================================================================
    PANGENOME PRIVATE JOIN MODULE
========================================================================================
    Repo location: modules/pangenome_private_join.nf

    Collects every haplotype's PRIVATE_MAP and PRIVATE_KMER table for one graph flavour --
    each already carrying both the PRIVATE set and its matched NON-PRIVATE CONTROL -- into a
    single per-segment table, a cross-tabulation, and a flat CSV for modelling.

    WHY THE CROSS-TABULATION IS THE DELIVERABLE
    -------------------------------------------
    PRIVATE_MAP and PRIVATE_KMER answer "is this really private?" by independent routes: one
    aligns the segment against the other assemblies, the other counts its k-mers in the
    sample's own reads. Neither sees the other's evidence, so the four combinations each mean
    something specific:

      NOT_PRIVATE       + REPEAT_LIKE   present elsewhere AND high copy -> graph collapse
      PRIVATE_CONFIRMED + UNIQUE_LIKE   absent elsewhere AND single copy -> novel sequence
      NOT_PRIVATE       + UNIQUE_LIKE   present elsewhere but single copy -> the graph failed
                                        to merge homologous sequence; an alignment failure
                                        rather than a repeat problem
      PRIVATE_CONFIRMED + REPEAT_LIKE   absent elsewhere but high copy WITHIN this sample ->
                                        haplotype-specific expansion

    A fifth verdict, NO_ALIGNMENT, marks segments that aligned nowhere at all -- not even to
    their own assembly. Those are usually low-complexity, and minimap2 omits such a query
    entirely; enumerating them from the FASTA rather than the PAF is what keeps them from
    looking like a failed task.

    WHY THE CSV
    -----------
    The control set makes enrichment testable: is private sequence more repeat-like than
    size-matched non-private sequence from the same haplotype and contig? That is a model, not
    a threshold count, so the CSV carries one row per segment with is_private as 0/1,
    repeat_like as 0/1, and median_copy / copy_ratio for a continuous response -- enough for
    either a binomial or a lognormal fit with haplotype, chromosome and their interaction as
    random effects, and log(span_bp) as a fixed effect because multiplicity correlates with
    length.

    WHAT THE GUARD DOES AND DOES NOT CATCH
    --------------------------------------
    A WHOLE HAPLOTYPE absent from either stream is a failed task and is fatal. A handful of
    individual segments is not. The original check conflated the two and killed a run over 41
    unalignable segments while all 20 tasks had succeeded.

    Input : tuple(taxid, flavor, map_tables, kmer_tables), script
    Output: joined table / CSV / cross-tab / audit / versions
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
          path("${taxid}.${flavor}.private_evidence.csv"),       emit: csv
    tuple val(taxid), val(flavor),
          path("${taxid}.${flavor}.private_evidence_xtab.tsv"),   emit: xtab
    tuple val(taxid), val(flavor),
          path("${taxid}.${flavor}.private_evidence_audit.tsv"),  emit: audit
    path("versions.tsv"),                                         emit: versions

    script:
    def bins = params.pangenome_sv_bins ?: '50,100,500,1000,5000,10000,50000,100000,500000,1000000'
    """
    set -euo pipefail

    python3 ${script} \\
        --map ${map_tables} \\
        --kmer ${kmer_tables} \\
        --label ${taxid}.${flavor} \\
        --flavor ${flavor} \\
        --outdir . \\
        --bins '${bins}'

    A=${taxid}.${flavor}.private_evidence_audit.tsv

    # A WHOLE HAPLOTYPE absent from either stream means a failed MAP or KMER task, and the
    # cross-tabulation would be biased because the missing rows are not random with respect to
    # verdict. Fatal.
    miss_m=\$(awk -F'\\t' '\$1=="haplotypes_missing_from_map"{print \$2}'  "\$A")
    miss_k=\$(awk -F'\\t' '\$1=="haplotypes_missing_from_kmer"{print \$2}' "\$A")
    if [ "\${miss_m:--}" != "-" ] || [ "\${miss_k:--}" != "-" ]; then
        echo "[PRIVATE_JOIN ${taxid}:${flavor}] ERROR: whole haplotypes missing." >&2
        echo "  absent from map : \${miss_m}" >&2
        echo "  absent from kmer: \${miss_k}" >&2
        exit 1
    fi

    # Individual segments missing from the map stream are NOT an error: minimap2 omits a query
    # that aligns nowhere at all. Reported so the count is visible.
    nm=\$(awk -F'\\t' '\$1=="segments_without_map_row"{print \$2}' "\$A")
    if [ "\${nm:-0}" -gt 0 ]; then
        echo "[PRIVATE_JOIN ${taxid}:${flavor}] NOTE: \${nm} segments have no map row" >&2
        echo "  (aligned nowhere, incl. to their own assembly). Reported as NO_MAP_ROW." >&2
    fi

    joined=\$(awk -F'\\t' '\$1=="joined_segments"{print \$2}' "\$A")
    echo "[PRIVATE_JOIN ${taxid}:${flavor}] \${joined} segments joined" >&2

    {
      printf 'process\\ttool\\tversion\\n'
      printf '%s\\tpython\\t%s\\n' "${task.process}" "\$(python3 --version 2>&1 | awk '{print \$2}')"
    } > versions.tsv
    """

    stub:
    """
    S=${taxid}.${flavor}
    printf 'haplotype\\tsample\\tsegment\\tcontig\\tstart\\tend\\tsegment_bp\\tn_other_assemblies\\tbest_identity\\taligned_frac_merged\\tmap_verdict\\tn_kmers\\tmean_copy\\tmedian_copy\\tfrac_absent\\tkmer_verdict\\tcombined\\n' > \$S.private_evidence.tsv
    printf 'haplotype,individual,sample,flavor,chromosome,segment,start,end,span_bp,log_span,is_private,n_other_assemblies,best_identity,aligned_frac_merged,map_verdict,n_kmers_observed,n_kmers_expected,frac_absent,mean_copy,median_copy,max_copy,single_copy_ref,copy_ratio,kmer_verdict,repeat_like,combined\\n' > \$S.private_evidence.csv
    printf 'scope\\tkey\\tcombined\\tn_segments\\tsegment_bp\\n' > \$S.private_evidence_xtab.tsv
    printf 'metric\\tvalue\\nhaplotypes_missing_from_map\\t-\\nhaplotypes_missing_from_kmer\\t-\\nsegments_without_map_row\\t0\\njoined_segments\\t0\\n' > \$S.private_evidence_audit.tsv
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
