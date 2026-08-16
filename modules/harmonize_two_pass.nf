/*
========================================================================================
    TWO-PASS REFERENCE SELECTION
========================================================================================
    Repo location: modules/harmonize_two_pass.nf

    Reference selection used to guess from scaffold counts in bash: count scaffolds above
    a size threshold, take the assembly nearest the batch median, break ties on N50. Four
    problems -- it optimised a number the harmonizer never uses (bash nchrom() and the
    drop-off caller disagree by 0-2 on good assemblies and ~240 on failed ones); nearest-
    median assumes the batch is centred on truth, and two failed assemblies moved the
    median from 17.5 to 19.5 and chose a reference split at five junctions; passengers
    voted; and no quality dimension was visible at all.

    A candidate can only really be judged by running the alignment and the join graph
    against it. So: score every candidate that way in parallel, then harmonize once with
    the winner.

        HARMONIZE_CANDIDATES   fais only, one task    -> ranked candidate list
        HARMONIZE_SCORE        one task per candidate -> one metrics row each
        HARMONIZE_SELECT       collects the rows      -> reference_id + ranked table

    THE METRIC (validated on four references spanning 15-26 pieces before being built)

        multiplicity = placed_scaffolds(other voters) / (n_consensus x n_other_voters)

        1.0 is a clean 1:1 chromosome mapping. Above 1, assemblies contribute several
        scaffolds per chromosome. Below 1, chromosomes are missing. The reference's OWN
        scaffolds are excluded -- they define the frame and are trivially placed, so
        counting them measures which assembly was picked rather than how good the frame is.

        Rule: voter -> min |multiplicity - 1| -> fewest graph edges -> highest
        genome_fraction. Edges is the tie-break and not the primary because a reference
        that wrongly FUSES two chromosomes needs no edges to correct and would score best
        on edges alone.

    HARMONIZE_SCORE runs the identical alignment and harmonizer invocation a real run
    would -- same preset, same arguments, every assembly aligned -- so a score is a
    faithful preview of that frame rather than an approximation of it. The arguments are
    built once in workflows/harmonize_scaffolds.nf and passed to both this and
    HARMONIZE_SPECIES, so the two cannot drift apart.

    Cost: n_candidates x n_assemblies alignments, fully parallel across candidates.
========================================================================================
*/

nextflow.enable.dsl = 2

process HARMONIZE_CANDIDATES {
    tag "taxid_${taxid}"
    label 'harmonize_candidates'

    publishDir "${params.outdir}/assembly/harmonization", mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.tsv' ? null : fn }

    input:
    tuple val(taxid), val(ids), path(assemblies)
    path(resolver)
    val(hargs)

    output:
    tuple val(taxid), path("${taxid}.reference_candidates.tsv"), emit: candidates
    path("versions.tsv"),                                        emit: versions

    script:
    def ncand = (params.harmonize_reference_candidates != null) ? params.harmonize_reference_candidates : 8
    """
    set -euo pipefail

    ids=(${ids.join(' ')})
    fastas=(${assemblies})
    for i in "\${!ids[@]}"; do
        ln -sf "\${fastas[\$i]}" "\${ids[\$i]}.input.fasta"
        samtools faidx "\${ids[\$i]}.input.fasta"
    done

    # --list-candidates reads ONLY the fais, so no alignment is needed to decide which
    # assemblies are worth aligning against. --reference-id is required by the manifest
    # reader but is not used for role assignment in this mode, so no candidate can be
    # promoted into the pool by being named here.
    printf 'id\\trole\\tfai\\tpaf\\n' > manifest.tsv
    for i in "\${!ids[@]}"; do
        R=query; [ "\$i" = "0" ] && R=ref
        printf '%s\\t%s\\t%s\\t-\\n' "\${ids[\$i]}" "\$R" "\${ids[\$i]}.input.fasta.fai" \\
            >> manifest.tsv
    done

    python3 ${resolver} \\
        --manifest manifest.tsv \\
        --reference-id "\${ids[0]}" \\
        --species "${taxid}" \\
        --list-candidates ${ncand} \\
        ${hargs} \\
        --outdir .

    N=\$(grep -vc '^#' "${taxid}.reference_candidates.tsv" || true)
    if [ "\${N}" -lt 2 ]; then
        echo "[HARMONIZE taxid ${taxid}] ERROR: only \$((N-1)) reference candidate(s)." >&2
        exit 1
    fi
    echo "[HARMONIZE taxid ${taxid}] \$((N-1)) reference candidate(s) to score"

    ST=\$(samtools --version 2>&1 | awk 'NR==1{print \$2}')
    PY=\$(python3 --version 2>&1 | awk '{print \$2}')
    {
        printf 'process\\ttool\\tversion\\n'
        printf '%s\\tsamtools\\t%s\\n' "${task.process}" "\${ST}"
        printf '%s\\tpython\\t%s\\n'   "${task.process}" "\${PY}"
    } > versions.tsv
    """

    stub:
    """
    printf 'candidate\\trole\\tn_chrom_set\\tgenome_fraction\\tcut_ratio\\tscaffold_n50\\tchrom_set_method\\treason\\n' \\
        > "${taxid}.reference_candidates.tsv"
    printf '%s\\tvoter\\t0\\t0\\t0\\t0\\tdropoff\\t-\\n' "${ids[0]}" >> "${taxid}.reference_candidates.tsv"
    printf '%s\\tvoter\\t0\\t0\\t0\\t0\\tdropoff\\t-\\n' "${ids[1]}" >> "${taxid}.reference_candidates.tsv"
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}


process HARMONIZE_SCORE {
    tag "taxid_${taxid}:${cand}"
    label 'harmonize_scaffolds'

    publishDir "${params.outdir}/assembly/harmonization/reference_scores",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.tsv' ? null : fn }

    input:
    tuple val(taxid), val(cand), val(ids), path(assemblies)
    path(resolver)
    val(hargs)

    output:
    tuple val(taxid), path("${taxid}.${cand}.score.tsv"), emit: score
    path("versions.tsv"),                                 emit: versions

    script:
    def preset = params.harmonize_minimap2_preset ?: 'asm5'
    """
    set -euo pipefail

    ids=(${ids.join(' ')})
    fastas=(${assemblies})
    for i in "\${!ids[@]}"; do
        ln -sf "\${fastas[\$i]}" "\${ids[\$i]}.input.fasta"
        samtools faidx "\${ids[\$i]}.input.fasta"
    done

    # identical to a real run against this reference: same preset, same flags, every
    # assembly aligned. A score is a preview of the frame, not an approximation of it.
    printf 'id\\trole\\tfai\\tpaf\\n' > manifest.tsv
    printf '%s\\tref\\t%s\\t-\\n' "${cand}" "${cand}.input.fasta.fai" >> manifest.tsv
    for id in "\${ids[@]}"; do
        [ "\$id" = "${cand}" ] && continue
        minimap2 -t ${task.cpus} -x ${preset} --eqx --secondary=no \\
            "${cand}.input.fasta" "\${id}.input.fasta" \\
            2> "\${id}.minimap2.log" > "\${id}.paf"
        printf '%s\\tquery\\t%s\\t%s\\n' "\$id" "\${id}.input.fasta.fai" "\${id}.paf" \\
            >> manifest.tsv
    done

    python3 ${resolver} \\
        --manifest manifest.tsv \\
        --reference-id "${cand}" \\
        --species "${taxid}" \\
        --score-only \\
        ${hargs} \\
        --outdir .

    MM=\$(minimap2 --version 2>&1 | awk 'NR==1{print}')
    ST=\$(samtools --version 2>&1 | awk 'NR==1{print \$2}')
    PY=\$(python3 --version 2>&1 | awk '{print \$2}')
    {
        printf 'process\\ttool\\tversion\\n'
        printf '%s\\tminimap2\\t%s\\n' "${task.process}" "\${MM}"
        printf '%s\\tsamtools\\t%s\\n' "${task.process}" "\${ST}"
        printf '%s\\tpython\\t%s\\n'   "${task.process}" "\${PY}"
    } > versions.tsv
    """

    stub:
    """
    printf 'reference_id\\trole\\tn_reference_pieces\\tn_consensus_chromosomes\\tgraph_rule\\tn_edges\\tn_voters\\tn_other_voters\\tplaced_other_voters\\tmultiplicity\\tabs_multiplicity_dev\\tgenome_fraction\\tcut_ratio\\tchrom_set_method\\tchrom_set_flags\\tscaffold_n50\\n' \\
        > "${taxid}.${cand}.score.tsv"
    printf '${cand}\\tvoter\\t0\\t0\\tplurality\\t0\\t0\\t0\\t0\\t1.0\\t0.0\\t0.0\\t0\\tdropoff\\t-\\t0\\n' \\
        >> "${taxid}.${cand}.score.tsv"
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}


process HARMONIZE_SELECT {
    tag "taxid_${taxid}"
    label 'harmonize_candidates'

    publishDir "${params.outdir}/assembly/harmonization", mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.tsv' ? null : fn }

    input:
    tuple val(taxid), path(scores)
    path(selector)

    output:
    tuple val(taxid), path("${taxid}.chosen_reference.txt"), emit: reference_id
    path("${taxid}.reference_scores.tsv"),                   emit: scores
    path("versions.tsv"),                                    emit: versions

    script:
    """
    set -euo pipefail

    python3 ${selector} \\
        --scores ${scores} \\
        --out-id "${taxid}.chosen_reference.txt" \\
        --out-table "${taxid}.reference_scores.tsv"

    if [ ! -s "${taxid}.chosen_reference.txt" ]; then
        echo "[HARMONIZE taxid ${taxid}] ERROR: no reference selected" >&2
        exit 1
    fi
    echo "[HARMONIZE taxid ${taxid}] two-pass reference = \$(cat ${taxid}.chosen_reference.txt)"

    PY=\$(python3 --version 2>&1 | awk '{print \$2}')
    {
        printf 'process\\ttool\\tversion\\n'
        printf '%s\\tpython\\t%s\\n' "${task.process}" "\${PY}"
    } > versions.tsv
    """

    stub:
    """
    printf 'STUB_REFERENCE\\n' > "${taxid}.chosen_reference.txt"
    printf 'reference_id\\trole\\tselected\\trank\\n' > "${taxid}.reference_scores.tsv"
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
