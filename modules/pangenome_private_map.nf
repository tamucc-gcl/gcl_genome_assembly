/*
========================================================================================
    PANGENOME PRIVATE MAP MODULE
========================================================================================
    Repo location: modules/pangenome_private_map.nf

    Aligns ONE haplotype's private segments against the combined index of all assemblies and
    asks, per segment: does this sequence exist in any OTHER assembly?

    WHY THIS IS THE MOST IMPORTANT CHECK IN THE PRIVATE-SEQUENCE WORKSTREAM
    ----------------------------------------------------------------------
    "Private" is defined by GRAPH COVERAGE -- a node walked by exactly one haplotype. That is
    not the same claim as "this sequence is absent from the other assemblies". A segment can be
    private in the graph because minigraph-cactus collapsed a repeat, or because the aligner
    simply failed to merge it, while the sequence sits plainly in three other assemblies.
    Every private-sequence figure the pipeline reports rests on that distinction, and nothing
    else can test it.

    The live question: chr8 and chr9 carry 23.6% and 22.9% private sequence against an ~11%
    floor across the other chromosomes, consistently across all ten haplotypes and BOTH graph
    flavours. Consistency across independent assemblies argues for biology; if those segments
    align cleanly to other assemblies, it is collapse.

    SCATTERED PER HAPLOTYPE
    -----------------------
    10 haplotypes x 2 flavours = 20 tasks, each a single minimap2 query against the shared
    index from PANGENOME_PRIVATE_INDEX. Pairwise would be 90 alignments per flavour, each
    rebuilding a ~1 Gb index.

    SELF-EXCLUSION NEEDS AN EXPLICIT ASSEMBLY ID
    --------------------------------------------
    A private segment is present in its own assembly by construction and aligns there at
    ~100%. Those hits must be dropped -- but by ASSEMBLY TAG, not by identity, because a
    genuine near-identical copy in another assembly is exactly the collapsed-repeat case being
    looked for and an identity filter would hide it.

    The assembly id is therefore passed in as a value, not derived from the PanSN haplotype
    name. It cannot be derived: the workflow builds names with
        sampleOf = id -> id.replaceFirst(/_hap[0-9]+$/,'').replaceAll(/\\./,'_')
    so dots become underscores (not invertible), and the reference individual's haplotypes get
    a different naming shape from everyone else's. Parsing would appear to work and would
    break silently the moment a different individual became the reference.

    PRESET must match PANGENOME_PRIVATE_INDEX -- a minimap2 index bakes in k and w, and a
    mismatched query preset is an error. Both read pangenome_private_map_preset.

    Input : tuple(taxid, flavor, haplotype, assembly_id, fasta), index, script
    Output: per-segment table / audit / versions
========================================================================================
*/

process PANGENOME_PRIVATE_MAP {
    tag "${taxid}:${flavor}:${haplotype}"
    label 'pangenome_private_map'

    publishDir "${params.outdir}/pangenome/${taxid}/private", mode: params.publish_dir_mode

    input:
    tuple val(taxid), val(flavor), val(haplotype), val(assembly_id), path(fasta)
    tuple val(itaxid), path(index)
    path(script)

    output:
    tuple val(taxid), val(flavor), val(haplotype), path("*.private_map.tsv"), emit: table
    tuple val(taxid), val(flavor), path("*.private_map_audit.tsv"),           emit: audit
    path("versions.tsv"),                                                    emit: versions

    script:
    def preset  = params.pangenome_private_map_preset ?: 'asm10'
    def minid   = params.pangenome_private_map_min_identity ?: 0.90
    def minfrac = params.pangenome_private_map_min_frac ?: 0.5
    """
    set -euo pipefail

    if [ ! -s ${fasta} ]; then
        echo "[PRIVATE_MAP ${taxid}:${flavor}:${haplotype}] empty FASTA; nothing to map" >&2
    fi

    # --secondary=yes deliberately: a segment present in several places is the collapsed-repeat
    # signal, and suppressing secondaries would hide exactly that. -N raised for the same reason.
    minimap2 -x ${preset} -t ${task.cpus} --secondary=yes -N 50 \\
        ${index} ${fasta} 2> minimap2.log | gzip -c > hits.paf.gz

    python3 ${script} \\
        --paf hits.paf.gz \\
        --self-assembly ${assembly_id} \\
        --haplotype '${haplotype}' \\
        --label ${taxid}.${flavor} \\
        --outdir . \\
        --min-identity ${minid} \\
        --min-frac ${minfrac}

    # A run with zero self-hits means the --self-assembly value does not match the index's
    # tags, which would make every segment look NOT_PRIVATE. Fail rather than publish that.
    A=\$(ls *.private_map_audit.tsv | head -1)
    rows=\$(awk -F'\\t' '\$1=="paf_rows"{print \$2}' "\$A")
    self=\$(awk -F'\\t' '\$1=="rows_self_excluded"{print \$2}' "\$A")
    if [ "\${rows:-0}" -gt 0 ] && [ "\${self:-0}" -eq 0 ]; then
        echo "[PRIVATE_MAP ${taxid}:${flavor}:${haplotype}] ERROR: ${assembly_id} matched no" >&2
        echo "  index tags -- self-hits were not excluded, so results are meaningless." >&2
        exit 1
    fi

    rm -f hits.paf.gz

    {
      printf 'process\\ttool\\tversion\\n'
      printf '%s\\tminimap2\\t%s\\n' "${task.process}" "\$(minimap2 --version 2>&1 | head -n1)"
      printf '%s\\tpython\\t%s\\n'   "${task.process}" "\$(python3 --version 2>&1 | awk '{print \$2}')"
    } > versions.tsv
    """

    stub:
    """
    S=${taxid}.${flavor}.stub
    printf 'haplotype\\tsegment\\tcontig\\tstart\\tend\\tsegment_bp\\tn_other_assemblies\\tbest_identity\\taligned_frac_merged\\tmax_frac_single_assembly\\tverdict\\n' > \$S.private_map.tsv
    printf 'metric\\tvalue\\npaf_rows\\t0\\nrows_self_excluded\\t0\\n' > \$S.private_map_audit.tsv
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
