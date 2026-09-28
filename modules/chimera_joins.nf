/* Validate and classify joins from the last scaffolding AGP against pre-finishing reference PAFs. Earlier joins are unresolved; no AGP chain crosses correction. The output table is bound to its assessment FASTA by SHA-256. */

process CHIMERA_JOINS {
    tag "${asm_id}"
    label 'chimera_joins'

    publishDir "${params.outdir}/assembly/chimeras", mode: params.publish_dir_mode

    input:
    tuple val(taxid), val(asm_id), path(round1_agp), path(round2_agp), path(ref_paf),
          path(candidates), path(ref_name_map), path(assessment_fasta, stageAs: 'assessment/*')
    path(agp_script)
    path(joins_script)
    path(coordinate_guard)

    output:
    tuple val(taxid), val(asm_id), path("${asm_id}.agp_joins.tsv"),      emit: agp_joins
    tuple val(taxid), val(asm_id), path("${asm_id}.chimeric_joins.tsv"), emit: called
    tuple val(taxid), val(asm_id), path("${asm_id}.coordinate_audit.tsv"), emit: coordinate_audit
    path("versions.tsv"),                                                emit: versions

    script:
    // Read from the AGP rather than assumed: every round-2 gap on this cohort is
    // `scaffold proximity_ligation`, 100 bp. A different scaffolder simply yields no matches
    // and every candidate goes to REVIEW rather than being cut at a guessed position.
    def gap_ev  = params.chimera_gap_evidence ?: 'proximity_ligation'
    def gap_len = params.chimera_gap_len ?: 100
    def minblk  = params.chimera_paf_min_block ?: 2000
    def cmin    = params.chimera_component_min_bp ?: 100000
    def cmarg   = params.chimera_component_margin ?: 2.0
    def maxdist = params.chimera_max_join_distance ?: 250000
    """
    set -euo pipefail

    R2=""
    if [ "${round2_agp.name}" != "NO_ROUND2" ] && [ -s "${round2_agp.name}" ]; then
        R2="--round2 ${round2_agp}"
    else
        echo "[CHIMERA_JOINS ${asm_id}] no round-2 AGP: round-1 objects are final" >&2
    fi

    python3 ${agp_script} \\
        --round1 ${round1_agp} \\
        \$R2 \\
        --assembly ${asm_id} \\
        --out ${asm_id}.agp_joins.tsv \\
        --require-evidence '${gap_ev}' \\
        --require-gap-len ${gap_len}

    # The reference has no PAF -- it is not aligned against itself -- and cannot be chimeric
    # with respect to its own coordinate system. A composite IN the reference is caught by
    # harmonization's consensus vote instead. Emit an empty table so the join downstream
    # still has a row for this assembly rather than dropping it.
    if [ "${ref_paf.name}" = "NO_PAF" ] || [ ! -s "${ref_paf.name}" ]; then
        echo "[CHIMERA_JOINS ${asm_id}] no reference PAF (this is the reference itself);" >&2
        echo "  emitting an empty called table" >&2
        printf '# no reference PAF: this assembly IS the reference\\n' \\
            > ${asm_id}.chimeric_joins.tsv
        printf 'assembly\\tscaffold\\tname\\tcut_bp\\tleft_chrom\\tright_chrom\\tleft_component\\tright_component\\tn_components\\tn_transitions\\tagp_join_bp\\tagp_join_distance\\tagp_source\\tgap_len\\tcallable\\treason\\tspan_bp\\tvote\\tcandidate_verdict\\n' \\
            >> ${asm_id}.chimeric_joins.tsv
    else
        python3 ${joins_script} \\
            --joins ${asm_id}.agp_joins.tsv \\
            --round1 ${round1_agp} \\
            \$R2 \\
            --paf ${ref_paf} \\
            --assembly ${asm_id} \\
            --scaffold-map ${candidates} \\
            --ref-name-map ${ref_name_map} \\
            --out ${asm_id}.chimeric_joins.tsv \\
            --min-block ${minblk} \\
            --component-min-bp ${cmin} \\
            --component-margin ${cmarg} \\
            --max-join-distance ${maxdist}
    fi

    python3 ${coordinate_guard} --fasta ${assessment_fasta} --agp ${round1_agp} \\
        --table ${asm_id}.chimeric_joins.tsv --audit ${asm_id}.coordinate_audit.tsv \\
        --round ${params.run_scaffold_round2 ? 'round2' : 'round1'}
    # Surface what was called, and what was NOT tested. A scaffold whose components could
    # none be assigned a chromosome is reported by the script as a warning -- it must not be
    # confused with a scaffold that was tested and came back clean.
    nc=\$(awk -F'\\t' 'NR>1 && \$1!~/^#/ && \$15=="yes"' ${asm_id}.chimeric_joins.tsv | wc -l)
    nb=\$(awk -F'\\t' 'NR>1 && \$1!~/^#/ && \$15=="yes" && \$19=="BREAK_CANDIDATE"' \\
          ${asm_id}.chimeric_joins.tsv | wc -l)
    echo "[CHIMERA_JOINS ${asm_id}] \${nc:-0} callable chimeric join(s), \${nb:-0} on a" >&2
    echo "  BREAK_CANDIDATE scaffold" >&2
    if [ "\${nb:-0}" -gt 0 ]; then
        awk -F'\\t' 'NR>1 && \$1!~/^#/ && \$15=="yes" && \$19=="BREAK_CANDIDATE" {
            printf "    %-16s cut %-12s %s -> %s  (%s)\\n", \$2, \$4, \$5, \$6, \$13 }' \\
            ${asm_id}.chimeric_joins.tsv >&2
    fi

    {
      printf 'process\\ttool\\tversion\\n'
      printf '%s\\tpython\\t%s\\n' "${task.process}" "\$(python3 --version 2>&1 | awk '{print \$2}')"
    } > versions.tsv
    """

    stub:
    """
    printf 'metric\\tvalue\\nstatus\\tstub_unvalidated\\n' > ${asm_id}.coordinate_audit.tsv
    printf 'assembly\\tfinal_object\\tfinal_cut\\tsource\\tlift\\tgap_len\\n' \\
      > ${asm_id}.agp_joins.tsv
    printf 'assembly\\tscaffold\\tname\\tcut_bp\\tleft_chrom\\tright_chrom\\tleft_component\\tright_component\\tn_components\\tn_transitions\\tagp_join_bp\\tagp_join_distance\\tagp_source\\tgap_len\\tcallable\\treason\\tspan_bp\\tvote\\tcandidate_verdict\\n' \\
      > ${asm_id}.chimeric_joins.tsv
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
