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
    def rthr = params.pangenome_private_kmer_repeat_threshold ?: 3.0
    def ncol = params.pangenome_private_kmer_name_column  ?: 0
    def vcol = params.pangenome_private_kmer_value_column ?: 0
    """
    set -euo pipefail

    if [ ! -s ${fasta} ]; then
        echo "[PRIVATE_KMER ${taxid}:${flavor}:${haplotype}] empty FASTA; nothing to look up" >&2
    fi

    # Record the tool's own interface in the log. meryl-lookup's -dump column layout has
    # varied between releases, and the summariser auto-detects it -- capturing --help here
    # means a future format change is diagnosable from the log rather than from guesswork.
    meryl-lookup --help > meryl_lookup_help.txt 2>&1 || true

    # Streamed, never stored: ~137 Mb of private sequence per haplotype is ~137M k-mer lines,
    # minutes to pipe but tens of GB to write out.
    meryl-lookup -dump \\
        -sequence ${fasta} \\
        -mers ${meryl_db} \\
      | python3 ${script} \\
            --haplotype '${haplotype}' \\
            --sample ${sample} \\
            --label ${taxid}.${flavor} \\
            --outdir . \\
            --repeat-threshold ${rthr} \\
            --name-column ${ncol} \\
            --value-column ${vcol}

    # A high absent fraction means the wrong meryl database was joined to this haplotype --
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
