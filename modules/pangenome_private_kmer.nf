/*
========================================================================================
    PANGENOME PRIVATE KMER MODULE
========================================================================================
    Repo location: modules/pangenome_private_kmer.nf

    Per segment, k-mer copy number in that haplotype's SAMPLE'S OWN read database -- for the
    PRIVATE set and for a size-matched NON-PRIVATE CONTROL from the same haplotype.

    NO EXTERNAL DEPENDENCY, BY DESIGN
    ---------------------------------
    Needs no TE library, no reference database and no annotation pipeline, so it works for any
    species the pipeline is pointed at. It reuses the Merqury meryl database that
    BUILD_MERYL_DB already produces for assembly QC.

    READ-DERIVED, NOT ASSEMBLY-DERIVED
    ----------------------------------
    Read multiplicity reflects true genomic copy number and is unaffected by whether the
    assembler collapsed the repeat. An assembly-derived count would be circular, because
    collapse is one of the things that may be producing spurious private sequence.

    WHY BOTH SETS RUN IN ONE TASK, CONTROL FIRST
    -------------------------------------------
    Copy number is read multiplicity, so its scale is sequencing depth: single-copy k-mers
    measured ~14x on this cohort while a private satellite segment read ~360,000. Ratios are
    therefore only meaningful against a single-copy reference -- and that reference must NOT
    come from the private set, because private sequence is repeat-enriched and the reference
    would absorb the signal under test. Deriving it that way gave 233-436 where the true level
    is ~14, and called two thirds of private segments repeat-like.

    So the CONTROL run derives the reference and the PRIVATE run is normalised against it.
    Splitting them into separate tasks would require an ordering dependency between two tasks
    of the same process, which Nextflow cannot express -- hence one task, two calls, and the
    reference handed from the first to the second. It also guarantees both sets share a
    reference, which is what makes copy_ratio comparable between them.

    THE SAMPLE JOIN CANNOT BE DERIVED FROM THE HAPLOTYPE NAME
    --------------------------------------------------------
    meryl databases are per SAMPLE (five here); segments are per HAPLOTYPE (ten). The workflow
    builds PanSN names with
        sampleOf = id -> id.replaceFirst(/_hap[0-9]+$/,'').replaceAll(/\\./,'_')
    so dots become underscores irreversibly and the reference individual's haplotypes take a
    different shape from everyone else's. Parsing a sample back out would look correct on this
    cohort and break the moment a different individual became the reference. The sample id is
    therefore passed in from the workflow's own mapping.

    Input : tuple(taxid, flavor, haplotype, sample, private_fa, control_fa, meryl_db), script
    Output: one table carrying BOTH sets / audits / versions
========================================================================================
*/

process PANGENOME_PRIVATE_KMER {
    tag "${taxid}:${flavor}:${haplotype}"
    label 'pangenome_private_kmer'

    publishDir "${params.outdir}/pangenome/${taxid}/private", mode: params.publish_dir_mode

    input:
    tuple val(taxid), val(flavor), val(haplotype), val(sample),
          path(private_fa), path(control_fa), path(meryl_db)
    path(script)

    output:
    tuple val(taxid), val(flavor), val(haplotype),
          path("*.private_kmer.combined.tsv"),      emit: table
    tuple val(taxid), val(flavor),
          path("*.private_kmer_audit.tsv"),         emit: audit
    path("versions.tsv"),                           emit: versions

    script:
    // A MULTIPLE of the CONTROL's single-copy level, not a raw count -- see the header.
    def rratio  = params.pangenome_private_kmer_repeat_ratio ?: 3.0
    // k must match the meryl database: it sets the EXPECTED k-mer count per segment and hence
    // frac_absent, because -wig-count OMITS positions with no hit rather than zeroing them.
    def kmer    = params.kmer_size ?: 21
    def hapsafe = haplotype.replaceAll('#', '_').replaceAll('/', '_')
    def stem    = "${taxid}.${flavor}.${hapsafe}"
    """
    set -euo pipefail

    # meryl-lookup's report types are -bed, -bed-runs, -wig-count, -wig-depth, -existence,
    # -include and -exclude. There is NO -dump; an earlier version asked for one and every
    # task failed. -wig-count gives "the multiplicity of the kmer starting at each position",
    # which is the copy-number measure wanted. Capturing --help makes a future rename
    # diagnosable from the log rather than from guesswork.
    meryl-lookup -wig-count -help > meryl_lookup_help.txt 2>&1 || true

    # ---- CONTROL first: it derives the single-copy reference --------------------------
    meryl-lookup -wig-count -sequence ${control_fa} -mers ${meryl_db} \\
      | python3 ${script} \\
            --haplotype '${haplotype}' --sample ${sample} --set control \\
            --label ${taxid}.${flavor} --outdir . \\
            --kmer ${kmer} --repeat-ratio ${rratio}

    CA=${stem}.control.private_kmer_audit.tsv
    REF=\$(awk -F'\\t' '\$1=="single_copy_reference"{print \$2}' "\$CA")
    if [ -z "\${REF:-}" ] || awk -v r="\${REF:-0}" 'BEGIN{exit !(r+0 <= 0)}'; then
        echo "[PRIVATE_KMER ${taxid}:${flavor}:${haplotype}] ERROR: the control run produced" >&2
        echo "  no usable single-copy reference (got '\${REF:-}'). The private set cannot be" >&2
        echo "  normalised, and normalising it against itself is the bug this replaced." >&2
        exit 1
    fi
    echo "[PRIVATE_KMER ${taxid}:${flavor}:${haplotype}] single-copy reference \$REF (control)" >&2

    # ---- PRIVATE, normalised against the CONTROL reference ---------------------------
    meryl-lookup -wig-count -sequence ${private_fa} -mers ${meryl_db} \\
      | python3 ${script} \\
            --haplotype '${haplotype}' --sample ${sample} --set private \\
            --label ${taxid}.${flavor} --outdir . \\
            --kmer ${kmer} --repeat-ratio ${rratio} --single-copy "\$REF"

    # ---- one table carrying both sets ------------------------------------------------
    P=${stem}.private.private_kmer.tsv
    C=${stem}.control.private_kmer.tsv
    { grep -h '^#' "\$P"
      grep -v '^#' "\$P" | head -1
      grep -v '^#' "\$P" | tail -n +2
      grep -v '^#' "\$C" | tail -n +2
    } > ${stem}.private_kmer.combined.tsv

    np=\$(grep -vc '^#' "\$P" || true)
    nc=\$(grep -vc '^#' "\$C" || true)
    echo "[PRIVATE_KMER ${taxid}:${flavor}:${haplotype}] \$((np-1)) private + \$((nc-1)) control rows" >&2

    # High absence means the wrong meryl database was joined -- these are the sample's own
    # reads, so its own k-mers must be present. Measured against EXPECTED positions, since
    # -wig-count omits misses. Warn rather than fail: a genuinely low-coverage sample can
    # legitimately show some absence.
    for A in ${stem}.private.private_kmer_audit.tsv \$CA; do
        fa=\$(awk -F'\\t' '\$1=="overall_frac_absent"{print \$2}' "\$A")
        awk -v f="\${fa:-0}" -v a="\$A" 'BEGIN{
            if (f+0 > 0.5)
                printf("[PRIVATE_KMER] WARNING: %.1f%% of expected k-mers absent in %s -- check the sample join\\n",
                       100*f, a) > "/dev/stderr"
        }'
    done

    {
      printf 'process\\ttool\\tversion\\n'
      printf '%s\\tmeryl\\t%s\\n'  "${task.process}" "\$(meryl --version 2>&1 | sed 's/meryl //' | head -n1)"
      printf '%s\\tpython\\t%s\\n' "${task.process}" "\$(python3 --version 2>&1 | awk '{print \$2}')"
    } > versions.tsv
    """

    stub:
    """
    S=${taxid}.${flavor}.stub
    printf 'haplotype\\tsample\\tset\\tsegment\\tcontig\\tstart\\tend\\tspan_bp\\tn_kmers_observed\\tn_kmers_expected\\tfrac_absent\\tmean_copy\\tmedian_copy\\tmax_copy\\tsingle_copy_ref\\tcopy_ratio\\tverdict\\n' \\
      > \$S.private_kmer.combined.tsv
    printf 'metric\\tvalue\\nsingle_copy_reference\\t14\\noverall_frac_absent\\t0\\n' \\
      > \$S.private_kmer_audit.tsv
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
