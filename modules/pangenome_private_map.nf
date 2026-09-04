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

    Both the PRIVATE set and its size-matched NON-PRIVATE CONTROL run in one task, against
    the same index with the same parameters. The control is only interpretable if it is
    measured identically, and one task guarantees that.

    Input : tuple(taxid, flavor, haplotype, assembly_id, private_fa, control_fa), index, script
    Output: per-segment table / audit / versions
========================================================================================
*/

process PANGENOME_PRIVATE_MAP {
    tag "${taxid}:${flavor}:${haplotype}"
    label 'pangenome_private_map'

    publishDir "${params.outdir}/pangenome/${taxid}/private", mode: params.publish_dir_mode

    input:
    tuple val(taxid), val(flavor), val(haplotype), val(assembly_id),
          path(private_fa), path(control_fa)
    tuple val(itaxid), path(index)
    path(script)

    output:
    tuple val(taxid), val(flavor), val(haplotype),
          path("*.private_map.combined.tsv"), emit: table
    tuple val(taxid), val(flavor), path("*.private_map_audit.tsv"),           emit: audit
    path("versions.tsv"),                                                    emit: versions

    script:
    def preset  = params.pangenome_private_map_preset ?: 'asm10'
    // 0 = DERIVE from the control (the default). A non-zero value overrides.
    def minid    = params.pangenome_private_map_min_identity ?: 0
    def dpct     = params.pangenome_private_map_derive_percentile ?: 5.0
    def maxnohit = params.pangenome_private_map_max_control_no_hit ?: 0.5
    def minfrac = params.pangenome_private_map_min_frac ?: 0
    def hapsafe = haplotype.replaceAll('#', '_').replaceAll('/', '_')
    def stem    = "${taxid}.${flavor}.${hapsafe}"
    """
    set -euo pipefail

    # CONTROL FIRST, then PRIVATE with the CONTROL's thresholds.
    #
    # min_identity and min_frac are DERIVED from the control unless overridden. Control
    # windows are non-private by construction, so the identity distribution of their hits is
    # the empirical cross-haplotype homology identity for this cohort -- 0.80 fits this
    # clupeid, but a plant or copepod cohort would differ and a hardcoded value would be
    # wrong invisibly. A hand-set 0.90 discarded nearly all real homology here, so 70% of
    # control windows reported no homologue anywhere.
    #
    # Both sets are then filtered IDENTICALLY, which is the only way the contrast means
    # anything. Same task, so the ordering is guaranteed.
    minimap2 -x ${preset} -t ${task.cpus} --secondary=yes -N 50 \
        ${index} ${control_fa} 2>> minimap2.log | gzip -c > hits.control.paf.gz

    python3 ${script} \
        --paf hits.control.paf.gz --fasta ${control_fa} \
        --self-assembly ${assembly_id} --haplotype '${haplotype}' --set control \
        --label ${taxid}.${flavor} --outdir . \
        --min-identity ${minid} --min-frac ${minfrac} \
        --derive-percentile ${dpct} --max-control-no-hit-frac ${maxnohit}

    CA=${stem}.control.private_map_audit.tsv
    MID=\$(awk -F'\\t' '\$1=="derived_min_identity"{print \$2}' "\$CA")
    MFR=\$(awk -F'\\t' '\$1=="derived_min_frac"{print \$2}' "\$CA")
    if [ -z "\${MID:-}" ] || [ -z "\${MFR:-}" ]; then
        echo "[PRIVATE_MAP ${taxid}:${flavor}:${haplotype}] ERROR: control run wrote no" >&2
        echo "  thresholds; the private set cannot be filtered comparably." >&2
        exit 1
    fi
    echo "[PRIVATE_MAP ${taxid}:${flavor}:${haplotype}] thresholds from control:" >&2
    echo "  min_identity=\$MID min_frac=\$MFR" >&2

    minimap2 -x ${preset} -t ${task.cpus} --secondary=yes -N 50 \
        ${index} ${private_fa} 2>> minimap2.log | gzip -c > hits.private.paf.gz

    python3 ${script} \
        --paf hits.private.paf.gz --fasta ${private_fa} \
        --self-assembly ${assembly_id} --haplotype '${haplotype}' --set private \
        --label ${taxid}.${flavor} --outdir . \
        --min-identity "\$MID" --min-frac "\$MFR"

    # Zero self-hits means the --self-assembly value does not match the index tags, and every
    # segment would come back NOT_PRIVATE -- a wrong answer that looks like a strong result.
    for A in *.private_map_audit.tsv; do
        rows=\$(awk -F'\\t' '\$1=="paf_rows"{print \$2}' "\$A")
        self=\$(awk -F'\\t' '\$1=="rows_self_excluded"{print \$2}' "\$A")
        if [ "\${rows:-0}" -gt 0 ] && [ "\${self:-0}" -eq 0 ]; then
            echo "[PRIVATE_MAP ${taxid}:${flavor}:${haplotype}] ERROR: ${assembly_id} matched no" >&2
            echo "  index tags in \$A -- self-hits were not excluded, results are meaningless." >&2
            exit 1
        fi
    done

    # one table carrying both sets. Single awk, NO PIPES: `grep -v '^#' f | head -1` sends
    # SIGPIPE to grep when head closes, and `set -o pipefail` turns that into exit 141. That
    # is how PANGENOME_PRIVATE_KMER failed, and the same trap is already documented in
    # pangenome_untangle.nf about piping `odgi paths -L` into head.
    P=${stem}.private.private_map.tsv
    C=${stem}.control.private_map.tsv
    awk '
      FNR==1          { nf++ }
      /^#/            { if (nf==1) print; next }
      \$1=="haplotype" { if (!hdr) { print; hdr=1 } next }
                      { print; rows[\$2]++ }
      END             { for (s in rows) printf("[PRIVATE_MAP] %s rows: %d\\n", s, rows[s]) > "/dev/stderr" }
    ' "\$P" "\$C" > ${stem}.private_map.combined.tsv

    rm -f hits.private.paf.gz hits.control.paf.gz

    {
      printf 'process\\ttool\\tversion\\n'
      printf '%s\\tminimap2\\t%s\\n' "${task.process}" "\$(minimap2 --version 2>&1 | head -n1)"
      printf '%s\\tpython\\t%s\\n'   "${task.process}" "\$(python3 --version 2>&1 | awk '{print \$2}')"
    } > versions.tsv
    """

    stub:
    """
    S=${taxid}.${flavor}.stub
    printf 'haplotype\\tset\\tsegment\\tcontig\\tstart\\tend\\tsegment_bp\\tn_other_assemblies\\tbest_identity\\taligned_frac_merged\\tmax_frac_single_assembly\\tverdict\\n' > \$S.private_map.combined.tsv
    printf 'metric\\tvalue\\npaf_rows\\t0\\nrows_self_excluded\\t0\\n' > \$S.private_map_audit.tsv
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
