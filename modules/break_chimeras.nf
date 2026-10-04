/* Reviewed, FASTA-bound gap cuts before finishing. Concordance votes are diagnostic; automatic cutting is unavailable pending independent evidence calibration. */

process BREAK_CHIMERAS {
    tag "${meta.id}"
    label 'break_chimeras'

    publishDir "${params.outdir}/assembly/chimeras", mode: params.publish_dir_mode,
               pattern: "*.chimera_break_audit.tsv"

    input:
    // stageAs IS required here: the output is ${meta.id}.broken.fasta, so an input named
    // <id>.teloclip_extended.fasta could otherwise collide with it. But `.name` already
    // includes the staged prefix -- writing "input/${assembly_fasta.name}" produced
    // input/input/... Use the path object directly.
    tuple val(meta), path(assembly_fasta, stageAs: 'input/*'), path(name_map), path(candidates)
    path(script)

    output:
    tuple val(meta), path("${meta.id}.broken.fasta"), path("${meta.id}.broken_name_map.tsv"),
        emit: assemblies
    tuple val(meta), path("${meta.id}.chimera_break_audit.tsv"), emit: audit
    path("versions.tsv"), emit: versions

    script:
    // File mode selects reviewed rows; Python also rejects auto if called directly.
    def mode    = (params.chimera_break?.toString() == 'auto') ? 'auto' : 'file'
    def minpiece = params.chimera_min_piece_bp ?: 1000000
    """
    set -euo pipefail

    # NO_HARMONIZE means this species has a single assembly, so there is no concordance vote
    # and nothing to detect against. Pass through: the chimera arm requires multiple
    # assemblies by construction, the same condition harmonization needs.
    if [ "${name_map.name}" = "NO_HARMONIZE" ] || [ ! -s "${candidates.name}" ]; then
        echo "[BREAK_CHIMERAS ${meta.id}] no name map or no candidates; passing through" >&2
        cp ${assembly_fasta} ${meta.id}.broken.fasta
        cp "${name_map.name}" ${meta.id}.broken_name_map.tsv 2>/dev/null \\
            || printf 'old_name\\tnew_name\\torient\\torder\\tlength\\tclass\\tref_span\\tflags\\n' \\
               > ${meta.id}.broken_name_map.tsv
        printf 'metric\\tvalue\\nassembly\\t${meta.id}\\nmode\\tpassthrough\\nscaffolds_broken\\t0\\n' \\
            > ${meta.id}.chimera_break_audit.tsv
        printf 'process\\ttool\\tversion\\n' > versions.tsv
        exit 0
    fi

    python3 ${script} \\
        --fasta ${assembly_fasta} \\
        --name-map ${name_map} \\
        --candidates ${candidates} \\
        --assembly ${meta.id} \\
        --out-fasta ${meta.id}.broken.fasta \\
        --out-name-map ${meta.id}.broken_name_map.tsv \\
        --audit ${meta.id}.chimera_break_audit.tsv \\
        --mode ${mode} \\
        --min-piece-bp ${minpiece}

    # SEQUENCE MUST BE CONSERVED. A split moves bases between records; it must never lose or
    # duplicate one. Compared on total non-header characters, because the record count and the
    # names both change by design.
    IN=\$(grep -v '^>' ${assembly_fasta} | tr -d '\\n' | wc -c)
    OUT=\$(grep -v '^>' ${meta.id}.broken.fasta | tr -d '\\n' | wc -c)
    if [ "\$IN" != "\$OUT" ]; then
        echo "[BREAK_CHIMERAS ${meta.id}] ERROR: sequence not conserved -- \$IN in, \$OUT out." >&2
        echo "  A split must only move bases between records, never lose or duplicate them." >&2
        exit 1
    fi
    echo "[BREAK_CHIMERAS ${meta.id}] sequence conserved: \$IN bp" >&2

    # Every name in the rewritten map must exist in the rewritten FASTA, or FINALIZE_ASSEMBLY
    # will silently drop records when it extracts by name.
    grep '^>' ${meta.id}.broken.fasta | sed 's/^>//; s/[[:space:]].*//' | sort -u > fa_names.txt
    awk -F'\\t' 'NR>1 && \$1!~/^#/{print \$1}' ${meta.id}.broken_name_map.tsv | sort -u > map_names.txt
    if ! miss=\$(comm -13 fa_names.txt map_names.txt) || [ -n "\${miss:-}" ]; then
        echo "[BREAK_CHIMERAS ${meta.id}] ERROR: name map references records absent from the" >&2
        echo "  FASTA; FINALIZE_ASSEMBLY would drop them:" >&2
        echo "\${miss}" | sed 's/^/    /' >&2
        exit 1
    fi

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
