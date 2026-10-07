/* Apply source-bound actions before finishing, with sequence reconstruction verification. */

process BREAK_CHIMERAS {
    tag "${meta.id}"
    label 'break_chimeras'

    publishDir "${params.outdir}/assembly/chimeras", mode: params.publish_dir_mode,
               pattern: "*.{chimera_break_audit.tsv,coordinate_lift.tsv,verification.json}"

    input:
    // stageAs IS required here: the output is ${meta.id}.broken.fasta, so an input named
    // <id>.teloclip_extended.fasta could otherwise collide with it. But `.name` already
    // includes the staged prefix -- writing "input/${assembly_fasta.name}" produced
    // input/input/... Use the path object directly.
    tuple val(meta), path(assembly_fasta, stageAs: 'input/*'), path(name_map), path(actions)
    path(script)

    output:
    tuple val(meta), path("${meta.id}.broken.fasta"), path("${meta.id}.broken_name_map.tsv"),
        emit: assemblies
    tuple val(meta), path("${meta.id}.chimera_break_audit.tsv"), emit: audit
    path("versions.tsv"), emit: versions
    tuple val(meta), path("*.coordinate_lift.tsv"), emit: coordinate_lift, optional: true
    tuple val(meta), path("*.verification.json"), emit: verification, optional: true

    script:
    // Manual and automatic selections use the same action format and applicator.
    def mode = 'file'
    def minpiece = params.chimera_min_piece_bp ?: 1000000
    def entryScript = script instanceof List ? script.find { it.name == 'break_chimeras.py' } : script
    """
    set -euo pipefail

    python3 ${entryScript} \\
        --fasta ${assembly_fasta} \\
        --name-map ${name_map} \\
        --actions ${actions} \\
        --assembly ${meta.id} \\
        --out-fasta ${meta.id}.broken.fasta \\
        --out-name-map ${meta.id}.broken_name_map.tsv \\
        --audit ${meta.id}.chimera_break_audit.tsv \\
        --mode ${mode} \\
        --min-piece-bp ${minpiece}

    nb=\$(awk -F'\\t' '\$1=="scaffolds_broken"{print \$2}' ${meta.id}.chimera_break_audit.tsv)
    echo "[BREAK_CHIMERAS ${meta.id}] mode=${mode}, broke \${nb:-0} scaffold(s)" >&2

    {
      printf 'process\\ttool\\tversion\\n'
      printf '%s\\tpython\\t%s\\n' "${task.process}" "\$(python3 --version 2>&1 | awk '{print \$2}')"
    } > versions.tsv
    """

    stub:
    """
    cp ${assembly_fasta} ${meta.id}.broken.fasta
    printf 'old_name\\tnew_name\\torient\\torder\\tlength\\tclass\\tref_span\\tflags\\n' \\
      > ${meta.id}.broken_name_map.tsv
    printf 'metric\\tvalue\\nscaffolds_broken\\t0\\n' > ${meta.id}.chimera_break_audit.tsv
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
