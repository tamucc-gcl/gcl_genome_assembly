/*
========================================================================================
    PANGENOME PRIVATE KMER MODULE
========================================================================================
    Repo location: modules/pangenome_private_kmer.nf

    Per private segment, the k-mer copy number in that haplotype's SAMPLE'S OWN read
    database. The pipeline's self-contained answer to whether private sequence is
    repeat-derived.

    NO EXTERNAL DEPENDENCY, BY DESIGN
    ---------------------------------
    This needs no TE library, no reference database and no annotation pipeline, so it works
    for any species the pipeline is pointed at. It reuses the Merqury meryl database that
    BUILD_MERYL_DB already produces for assembly QC, so it costs one meryl-lookup pass and
    nothing else. A RepeatMasker/RepeatModeler analysis can be layered on afterwards from a
    separate annotation pipeline if wanted, but nothing here waits on that.

    READ-DERIVED, NOT ASSEMBLY-DERIVED
    ----------------------------------
    Copy number comes from the sample's reads. Read multiplicity reflects true genomic copy
    number and is unaffected by whether the assembler collapsed the repeat -- whereas an
    assembly-derived count would be circular, because collapse is one of the things that may
    be producing spurious private sequence in the first place.

    A segment that is BOTH private and high copy number is a candidate for graph collapse
    rather than biology. PRIVATE_MAP reaches the same conclusion by aligning the segment
    against the other assemblies, which is an entirely independent route, so the two provide a
    genuine cross-check rather than two views of one measurement.

    THE SAMPLE JOIN CANNOT BE DERIVED FROM THE HAPLOTYPE NAME
    --------------------------------------------------------
    meryl databases are per SAMPLE (five here); private segments are per HAPLOTYPE (ten). The
    workflow builds PanSN names with
        sampleOf = id -> id.replaceFirst(/_hap[0-9]+$/,'').replaceAll(/\\./,'_')
        nm       = (sampleOf(id) == refInd) ? flat(id) : "${sampleOf(id)}.${hapOf(id)}"
    so the reference individual's haplotypes come out as Sde-CMat_203_hap1#0 while everyone
    else comes out as Sde-CBau_104#1, and dots have become underscores irreversibly. Parsing
    a sample back out of a PanSN name would appear to work on this cohort and break silently
    the moment a different individual were chosen as the reference. The sample id is therefore
    passed in as a value from the workflow's own mapping.

    The audit's overall_frac_absent is the tripwire for a mis-join: these are the sample's own
    reads, so a large fraction of absent k-mers means the wrong database was paired.

    THE THRESHOLD IS SELF-CALIBRATING. Copy number is read multiplicity, so its scale is
    sequencing depth -- measured ~14x single-copy on this cohort, against ~360,000 for a
    satellite segment at chr10 position 0. An absolute cutoff is therefore meaningless, and a
    shared constant cannot work across samples of differing depth. The script takes the
    single-copy level from the data and pangenome_private_kmer_repeat_ratio is a MULTIPLE of
    it; the derived level is written into the audit so the calibration is checkable.

    Input : tuple(taxid, flavor, haplotype, sample, fasta, meryl_db), script
    Output: per-segment table / audit / versions
========================================================================================
*/

process PANGENOME_PRIVATE_KMER {
    tag "${taxid}:${flavor}:${haplotype}"
    label 'pangenome_private_kmer'

    publishDir "${params.outdir}/pangenome/${taxid}/private", mode: params.publish_dir_mode

    input:
    tuple val(taxid), val(flavor), val(haplotype), val(sample), path(fasta), path(meryl_db)
    path(script)

    output:
    tuple val(taxid), val(flavor), val(haplotype), path("*.private_kmer.tsv"), emit: table
    tuple val(taxid), val(flavor), path("*.private_kmer_audit.tsv"),           emit: audit
    path("versions.tsv"),                                                     emit: versions

    script:
    // A MULTIPLE of this haplotype's own single-copy coverage, not a raw count. Read
    // multiplicity is sequencing depth: single-copy measured ~14x on this cohort while a
    // satellite segment read ~360,000, so an absolute threshold is meaningless and cannot be
    // shared across samples of differing depth. The script derives the single-copy level from
    // the data (kmer-weighted median of per-segment medians) and reports it for auditing.
    def rratio = params.pangenome_private_kmer_repeat_ratio ?: 3.0
    // k must match the meryl database, because it sets the EXPECTED k-mer count per segment
    // and hence frac_absent -- -wig-count omits positions with no hit rather than zeroing them.
    def kmer   = params.analysis_kmer ?: 21
    """
    set -euo pipefail

    if [ ! -s ${fasta} ]; then
        echo "[PRIVATE_KMER ${taxid}:${flavor}:${haplotype}] empty FASTA; nothing to look up" >&2
    fi

    # Record the tool's own interface in the log: the report-type names have changed between
    # meryl releases, and this file is what makes a future change diagnosable from the log.
    meryl-lookup -wig-count -help > meryl_lookup_help.txt 2>&1 || true

    # -wig-count, NOT -dump. meryl-lookup's report types are -bed, -bed-runs, -wig-count,
    # -wig-depth, -existence, -include and -exclude; there is no -dump, and asking for one
    # failed every task. -wig-count gives "the multiplicity of the kmer starting at each
    # position", which is the copy-number measure wanted; -existence gives only present/absent
    # counts. Output goes to stdout when -output is omitted, so this streams: ~137 Mb of
    # private sequence per haplotype is minutes to pipe and tens of GB to store.
    meryl-lookup -wig-count \\
        -sequence ${fasta} \\
        -mers ${meryl_db} \\
      | python3 ${script} \\
            --haplotype '${haplotype}' \\
            --sample ${sample} \\
            --label ${taxid}.${flavor} \\
            --outdir . \\
            --kmer ${kmer} \\
            --repeat-ratio ${rratio}

    # A high absent fraction means the wrong meryl database was joined to this haplotype --
    # measured against EXPECTED k-mer positions, since -wig-count omits misses. --
    # these are the sample's own reads, so its own k-mers must be present. Warn loudly; do not
    # fail, because a genuinely low-coverage sample could legitimately produce some absence.
    A=\$(ls *.private_kmer_audit.tsv | head -1)
    fa=\$(awk -F'\\t' '\$1=="overall_frac_absent"{print \$2}' "\$A")
    awk -v f="\${fa:-0}" -v h="${haplotype}" -v s="${sample}" 'BEGIN{
        if (f+0 > 0.5)
            printf("[PRIVATE_KMER] WARNING: %.1f%% of %s k-mers absent from %s reads -- check the sample join\\n",
                   100*f, h, s) > "/dev/stderr"
    }'

    {
      printf 'process\\ttool\\tversion\\n'
      printf '%s\\tmeryl\\t%s\\n'  "${task.process}" "\$(meryl --version 2>&1 | sed 's/meryl //' | head -n1)"
      printf '%s\\tpython\\t%s\\n' "${task.process}" "\$(python3 --version 2>&1 | awk '{print \$2}')"
    } > versions.tsv
    """

    stub:
    """
    S=${taxid}.${flavor}.stub
    printf 'haplotype\\tsample\\tsegment\\tcontig\\tstart\\tend\\tn_kmers\\tmean_copy\\tmedian_copy\\tfrac_absent\\tverdict\\n' > \$S.private_kmer.tsv
    printf 'metric\\tvalue\\noverall_frac_absent\\t0\\n' > \$S.private_kmer_audit.tsv
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
