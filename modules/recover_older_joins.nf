/* Direct source-flank recovery; never authorizes automatic cuts. */
process RECOVER_OLDER_JOINS {
    tag "${asm_id}"
    label 'older_join_recovery'
    publishDir "${params.outdir}/assembly/chimeras/older_joins", mode: params.publish_dir_mode,
        saveAs: { name -> name.endsWith('.source.json') ? null : name }
    input:
    tuple val(taxid), val(asm_id), path(original, stageAs: 'original/*'),
          path(original_agp, stageAs: 'original_agp/*'),
          path(current, stageAs: 'current/*'), path(native),
          path(ref_paf), path(candidates), path(ref_map)
    path resolver
    output:
    tuple val(taxid), val(asm_id), path("${asm_id}.review_joins.tsv"), emit: called
    tuple val(taxid), val(asm_id), path("${asm_id}.older_join_audit.tsv"), emit: audit
    tuple val(taxid), val(asm_id), path("${asm_id}.older_join_summary.json"), emit: summary
    path "${asm_id}.source.json", emit: source_manifest
    path 'versions.tsv', emit: versions
    script:
    """
    set -euo pipefail
    python3 ${resolver} prepare --prefix ${asm_id} \\
        --round1-agp ${original_agp} --round1-fasta ${original} \\
        --flank ${params.chimera_older_flank_bp} --minimum ${params.chimera_older_min_flank_bp}
    if [ -s ${asm_id}.flanks.fa ]; then
        minimap2 -x asm5 -c --eqx --secondary=yes -N 50 -p 0.5 \\
            --split-prefix older_split -t ${task.cpus} \\
            ${current} ${asm_id}.flanks.fa > ${asm_id}.flanks.paf
    else
        : > ${asm_id}.flanks.paf
    fi
    python3 ${resolver} assess --prefix ${asm_id} --assembly ${asm_id} \\
        --current ${current} --native ${native} --flank-paf ${asm_id}.flanks.paf \\
        --ref-paf ${ref_paf} --ref-map ${ref_map} --candidates ${candidates} \\
        --min-mapq ${params.chimera_older_min_mapq} --window ${params.chimera_older_window_bp} \\
        --min-span ${params.chimera_min_span} --min-bp ${params.chimera_component_min_bp} --margin ${params.chimera_component_margin}
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    printf 'older_join_recovery\\tminimap2\\t%s\\n' "\$(minimap2 --version)" >> versions.tsv
    """
}
