/*
========================================================================================
    PANGENOME PRIVATE INDEX MODULE
========================================================================================
    Repo location: modules/pangenome_private_index.nf

    Builds ONE minimap2 index over all finalized assemblies, with sequence names tagged by
    assembly id, for PRIVATE_MAP to query.

    WHY AN INDEX RATHER THAN PAIRWISE ALIGNMENT
    -------------------------------------------
    The question is whether a private segment exists in any OTHER assembly. Done pairwise
    that is 10 x 9 = 90 alignments per flavour, each building its own index over a ~1 Gb
    target. Against a single combined index it is 20 queries (10 haplotypes x 2 flavours)
    with no index build at all, and one target hit tells us WHICH assembly matched.

    Flavour-independent: the assemblies do not change between clip and full, so this is one
    task feeding both arms.

    WHY SEQUENCE NAMES MUST BE TAGGED
    ---------------------------------
    Harmonization gives every assembly the SAME chromosome names -- chr1_1, chr2_1 and so on.
    That is the point of harmonization, and it means a naive concatenation produces ten
    sequences called chr1_1 and a target name that identifies nothing. Names become
    <assembly_id>::<contig>, so PRIVATE_MAP can both attribute a hit and exclude the query's
    own assembly.

    Self-exclusion is why the tag is the ASSEMBLY ID and not the PanSN haplotype key: the
    assembly id comes straight off meta.id with no derivation, whereas recovering an assembly
    from a PanSN name would require inverting
        sampleOf = id -> id.replaceFirst(/_hap[0-9]+$/,'').replaceAll(/\\./,'_')
    which is not invertible (dots become underscores) and changes shape depending on which
    individual is the reference.

    PRESET
    ------
    A minimap2 index is preset-specific -- the k and w baked in at build time must match the
    query preset, or minimap2 errors. pangenome_private_map_preset therefore governs both, and
    is used here as well as in PRIVATE_MAP.

    asm10 by default rather than asm5: the question is whether homologous sequence exists
    elsewhere at all, so sensitivity matters more than precision, and a segment that only
    matches at 5% divergence is still not private. asm20 would be more sensitive again at the
    cost of more repeat-driven noise.

    Input : all finalized assemblies as tuple(ids, fastas), flavour-independent
    Output: tagged multi-FASTA / minimap2 index / audit / versions
========================================================================================
*/

process PANGENOME_PRIVATE_INDEX {
    tag "${taxid}"
    label 'pangenome_private_index'

    input:
    tuple val(taxid), val(asm_ids), path(fastas)

    output:
    tuple val(taxid), path("${taxid}.all_assemblies.mmi"), emit: index
    tuple val(taxid), path("${taxid}.all_assemblies.tsv"), emit: contents
    path("versions.tsv"),                                  emit: versions

    script:
    def preset = params.pangenome_private_map_preset ?: 'asm10'
    def ids    = asm_ids instanceof List ? asm_ids : [asm_ids]
    def fas    = fastas  instanceof List ? fastas  : [fastas]
    if( ids.size() != fas.size() )
        error "PANGENOME_PRIVATE_INDEX: ${ids.size()} ids but ${fas.size()} fastas -- the " +
              "channel lost its pairing, which would tag sequences with the wrong assembly"
    def pairs = [ids, fas].transpose().collect { i, f -> "${i}\t${f}" }.join('\n')
    """
    set -euo pipefail

    # id<TAB>fasta, generated from the paired channel rather than from filenames: a filename
    # is not a reliable assembly id after staging, and a mis-pairing here would silently
    # attribute every hit to the wrong assembly.
    cat > pairs.tsv <<'PAIRS'
${pairs}
PAIRS

    printf 'assembly_id\\tn_seqs\\tbp\\n' > ${taxid}.all_assemblies.tsv

    : > all.fa
    while IFS=\$'\\t' read -r aid fa; do
        [ -n "\$aid" ] || continue
        # prefix every sequence name with the assembly id. Harmonized assemblies share
        # chromosome names, so without this the target name is meaningless.
        awk -v a="\$aid" '/^>/{ sub(/^>/,""); split(\$0,p," "); printf(">%s::%s\\n",a,p[1]); next } {print}' \\
            "\$fa" >> all.fa
        n=\$(grep -c '^>' "\$fa" || true)
        b=\$(awk '!/^>/{t+=length(\$0)} END{print t+0}' "\$fa")
        printf '%s\\t%s\\t%s\\n' "\$aid" "\$n" "\$b" >> ${taxid}.all_assemblies.tsv
    done < pairs.tsv

    # duplicate target names would make hit attribution ambiguous; fail rather than guess
    # awk rather than `| head -5`: head closing early sends SIGPIPE up the pipeline, which
    # `set -o pipefail` converts into a task failure. Harmless while there are no duplicates
    # and fatal exactly when there are -- i.e. precisely when the check needs to report.
    dup=\$(grep '^>' all.fa | sort | uniq -d | awk 'NR<=5')
    if [ -n "\$dup" ]; then
        echo "[PRIVATE_INDEX ${taxid}] ERROR: duplicate target names after tagging:" >&2
        echo "\$dup" >&2
        exit 1
    fi

    echo "[PRIVATE_INDEX ${taxid}] \$(grep -c '^>' all.fa) target sequences from \$(wc -l < pairs.tsv) assemblies" >&2

    # -I forces a SINGLE-part index. Default is 4G, and ~10 Gb of assembly splits into two
    # parts (measured: 8.0 Gb + 2.4 Gb). minimap2 then maps against each part in turn and
    # computes mapping quality and secondary selection PER PART, so results differ from a
    # single-part index. Sized from the actual input rather than hardcoded.
    TOT=\$(awk '!/^>/{n+=length(\$0)} END{printf "%d", (n/1000000000)+2}' all.fa)
    echo "[PRIVATE_INDEX ${taxid}] indexing \${TOT}G in one part" >&2
    minimap2 -x ${preset} -t ${task.cpus} -I \${TOT}G -d ${taxid}.all_assemblies.mmi all.fa
    rm -f all.fa

    {
      printf 'process\\ttool\\tversion\\n'
      printf '%s\\tminimap2\\t%s\\n' "${task.process}" "\$(minimap2 --version 2>&1 | head -n1)"
    } > versions.tsv
    """

    stub:
    """
    : > ${taxid}.all_assemblies.mmi
    printf 'assembly_id\\tn_seqs\\tbp\\n' > ${taxid}.all_assemblies.tsv
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
