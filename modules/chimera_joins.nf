/*
========================================================================================
    CHIMERA JOINS MODULE
========================================================================================
    Repo location: modules/chimera_joins.nf

    Which scaffolding joins separate two different CONSENSUS chromosomes -- the chimeric
    ones -- and exactly where each sits.

    Runs both scripts, because they are never useful apart:
      agp_joins.py      both AGPs chained -> every join, exact, in final coordinates
      chimera_joins.py  AGP components + the reference PAF -> which joins are chimeric

    WHY BOTH, AND WHY NEITHER ALONE
    -------------------------------
    The AGP says where joins ARE but not which matter: Sde-CTlk_104_hap2 scaffold_3 has 173
    joins and one is chimeric; Sde-CPla_115_hap1 has 919 across the assembly. Breaking at
    every join would undo scaffolding entirely.

    Inference alone cannot give a cut point. Estimating the junction by binning the PAF put
    it within ~0.7 Mb on both known candidates -- close enough to look right, and landing in
    sequence rather than in the 100 bp scaffolding gap that the join actually is.

    CHAINING BOTH ROUNDS IS REQUIRED
    --------------------------------
    Round 1 makes 1,830-1,891 joins per haplotype here; round 2 makes 32-54. BOTH known
    junctions are round-1 joins, so without lifting them ~97% of the search space -- including
    every real answer -- is invisible. --round2 is optional: when round 2 did not run, the
    round-1 objects are final.

    THE REFERENCE NAME MAP IS NOT OPTIONAL
    --------------------------------------
    Harmonization keeps two namespaces: refN is a reference-frame PIECE, chrN a CONSENSUS
    chromosome from the join graph's connected components across voters. A consensus
    chromosome can span several reference pieces, so a scaffold joining ref5 and ref12 into
    chr5 is CORRECTLY JOINED, not chimeric -- and calling it chimeric would cut a good
    scaffold. Transitions must be measured in consensus chromosomes, which is what the
    reference's own name map supplies.

    It also fixes a namespace mismatch: harmonization's PAFs align the INPUT fastas so their
    targets are `scaffold_N`, while PAIRWISE_ALIGNMENT's align finalized assemblies so theirs
    are `chr5_1`. The map handles either, and an untranslatable target aborts rather than
    yielding zero transitions -- which would read as a clean result.

    DETECTION ONLY. Nothing is cut here. The concordance vote in chimera_candidates.tsv is
    the gate and BREAK_CHIMERAS applies it.

    Input : tuple(taxid, meta.id, round1_agp, round2_agp|NO_ROUND2, ref_paf|NO_PAF,
                  candidates, ref_name_map), agp_script, joins_script
    Output: per-assembly joins table, called-chimera table, versions
========================================================================================
*/

process CHIMERA_JOINS {
    tag "${asm_id}"
    label 'chimera_joins'

    publishDir "${params.outdir}/assembly/chimeras", mode: params.publish_dir_mode

    input:
    tuple val(taxid), val(asm_id), path(round1_agp), path(round2_agp), path(ref_paf),
          path(candidates), path(ref_name_map)
    path(agp_script)
    path(joins_script)

    output:
    tuple val(taxid), val(asm_id), path("${asm_id}.agp_joins.tsv"),      emit: agp_joins
    tuple val(taxid), val(asm_id), path("${asm_id}.chimeric_joins.tsv"), emit: called
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
    printf 'assembly\\tfinal_object\\tfinal_cut\\tsource\\tlift\\tgap_len\\n' \\
      > ${asm_id}.agp_joins.tsv
    printf 'assembly\\tscaffold\\tname\\tcut_bp\\tleft_chrom\\tright_chrom\\tleft_component\\tright_component\\tn_components\\tn_transitions\\tagp_join_bp\\tagp_join_distance\\tagp_source\\tgap_len\\tcallable\\treason\\tspan_bp\\tvote\\tcandidate_verdict\\n' \\
      > ${asm_id}.chimeric_joins.tsv
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
