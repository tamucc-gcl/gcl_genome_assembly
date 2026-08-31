/*
========================================================================================
    PANGENOME PRIVATE FASTA MODULE
========================================================================================
    Repo location: modules/pangenome_private_fasta.nf

    Extracts each haplotype's PRIVATE segments -- maximal runs of nodes covered by exactly
    one haplotype -- above a size floor, as one FASTA per haplotype.

    ONE TASK PER FLAVOUR, NOT PER HAPLOTYPE
    ---------------------------------------
    The GFA is 3.1 GB compressed and coverage must be computed over all 151.7M nodes before
    any haplotype's segments are known, so ten tasks would each redo the same whole-graph
    work. This reads it once and emits N FASTAs. PRIVATE_MAP and PRIVATE_KMER then scatter
    over those, 10 per flavour.

    RUN BOTH FLAVOURS
    -----------------
    Clipping removes 704,499,041 bp from the graph and 698,360,436 of that -- 99.1% -- is
    private sequence. Measured: clip yields 163,418 segments / 675 Mb at the 1 kb floor,
    full yields 170,311 / 1,371 Mb. The sub-1 kb population is identical between flavours, so
    clipping does not thin the distribution evenly; it removes a small number of very large
    private segments. It also inverts the reference's apparent rank -- highest of ten
    haplotypes on clip (15.08% of private bp), lowest of ten on full (8.08%) -- because the
    reference is the graph backbone and is never clipped.

    WHAT THE FASTAs ARE FOR
    -----------------------
    1. PRIVATE_MAP: a segment that aligns cleanly to several other assemblies is not private,
       it is graph collapse or an alignment failure. Nothing else in the pipeline can
       distinguish those two cases.
    2. PRIVATE_KMER: k-mer copy number against the sample's own Merqury database, as a
       repeat-content measure -- the pipeline's self-contained answer to whether private
       sequence is repeat-driven, independent of any external annotation.
    3. The live biological question: chr8 and chr9 carry 23.6% and 22.9% private sequence
       against an ~11% floor on other chromosomes, consistently across all ten haplotypes and
       both flavours. Real divergence, or collapsed repeat?

    Headers are <haplotype>|<contig>|<start>|<end>|seg<N> so downstream results join back to
    the segment table without a coordinate lookup.

    Input : tuple(taxid, flavor, gfa), script
    Output: per-haplotype PRIVATE and CONTROL FASTAs (flattened for scatter) / manifest /
            audit / versions

    The control windows are the other half of the deliverable. ~14% of a haplotype is private,
    so a control drawn without filtering would be ~14% contaminated toward the private value;
    candidate windows more than control_max_private_frac private are REJECTED rather than
    masked, because masking creates junction k-mers that exist nowhere in the genome and would
    be scored as absent. Sizes are drawn from the private distribution of the SAME
    (haplotype, contig), since multiplicity correlates with length and private sequence is not
    uniformly distributed across chromosomes.
========================================================================================
*/

process PANGENOME_PRIVATE_FASTA {
    tag "${taxid}:${flavor}"
    label 'pangenome_private_fasta'

    publishDir "${params.outdir}/pangenome/${taxid}/private", mode: params.publish_dir_mode

    input:
    tuple val(taxid), val(flavor), path(gfa)
    path(script)

    output:
    // flattened one-per-haplotype, so PRIVATE_MAP / PRIVATE_KMER scatter without a
    // groupTuple barrier on the consumer side
    tuple val(taxid), val(flavor), path("${taxid}.${flavor}.*.private.fa"),
        emit: fastas, optional: true
    // size-matched NON-PRIVATE windows from the same haplotype and contig. Without these the
    // private set has nothing to be compared against and "is private sequence repeat
    // enriched?" is not a testable question -- an earlier version calibrated against the
    // private segments themselves and produced a single-copy reference of 233-436 where the
    // true value is ~14, because the reference absorbed the signal under test.
    tuple val(taxid), val(flavor), path("${taxid}.${flavor}.*.control.fa"),
        emit: controls, optional: true
    tuple val(taxid), val(flavor), path("${taxid}.${flavor}.private_fasta_manifest.tsv"),
        emit: manifest
    tuple val(taxid), val(flavor), path("${taxid}.${flavor}.private_fasta_audit.tsv"),
        emit: audit
    path("versions.tsv"), emit: versions

    script:
    def minbp = params.pangenome_private_min_bp ?: 1000
    def ctrl  = (params.pangenome_private_control != false) ? '--control' : ''
    def cfrac = params.pangenome_private_control_max_private_frac ?: 0.05
    def cratio = params.pangenome_private_control_bp_ratio ?: 1.0
    def seed  = params.pangenome_private_control_seed ?: 42
    """
    set -euo pipefail

    python3 ${script} \\
        --gfa ${gfa} \\
        --label ${taxid}.${flavor} \\
        --outdir . \\
        --min-bp ${minbp} \\
        ${ctrl} \\
        --control-max-private-frac ${cfrac} \\
        --control-bp-ratio ${cratio} \\
        --seed ${seed}

    # The audit asserts internally that segment_bp == fasta_bp_total. It must ALSO match the
    # >=${minbp} bp rows of private_segment_spectrum.tsv from gfa_hap_coverage.py: the two
    # scripts compute coverage independently, so a disagreement means one of them is wrong
    # about what "private" means. Surfaced here rather than left for someone to notice.
    seg=\$(awk -F'\\t' '\$1=="segment_bp"{print \$2}'      ${taxid}.${flavor}.private_fasta_audit.tsv)
    fas=\$(awk -F'\\t' '\$1=="fasta_bp_total"{print \$2}'  ${taxid}.${flavor}.private_fasta_audit.tsv)
    if [ "\$seg" != "\$fas" ]; then
        echo "[PRIVATE_FASTA ${taxid}:${flavor}] ERROR: segment_bp=\$seg != fasta_bp_total=\$fas" >&2
        exit 1
    fi
    echo "[PRIVATE_FASTA ${taxid}:${flavor}] \$seg bp of private sequence extracted" >&2

    {
      printf 'process\\ttool\\tversion\\n'
      printf '%s\\tpython\\t%s\\n' "${task.process}" "\$(python3 --version 2>&1 | awk '{print \$2}')"
    } > versions.tsv
    """

    stub:
    """
    : > ${taxid}.${flavor}.stub_hap_1.private.fa
    : > ${taxid}.${flavor}.stub_hap_1.control.fa
    printf 'label\\thaplotype\\tfasta\\tn_segments\\tfasta_bp\\n' \\
      > ${taxid}.${flavor}.private_fasta_manifest.tsv
    printf 'metric\\tvalue\\nsegment_bp\\t0\\nfasta_bp_total\\t0\\n' \\
      > ${taxid}.${flavor}.private_fasta_audit.tsv
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
