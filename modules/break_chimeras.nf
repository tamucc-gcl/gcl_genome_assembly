/*
========================================================================================
    BREAK CHIMERAS MODULE
========================================================================================
    Repo location: modules/break_chimeras.nf

    Splits chimeric scaffolds at the detected junction and rewrites the harmonization name
    map so each half gets the chromosome name it deserves.

    WHERE IT SITS
    -------------
    Between HARMONIZE_SCAFFOLDS and FINALIZE_ASSEMBLY, and that is the only place it can go.
    Harmonization emits a NAME MAP, not a renamed FASTA -- FINALIZE_ASSEMBLY applies the map.
    So this is the last point at which the FASTA is still in original coordinates and the
    composite is still one record, and the first at which the cross-haplotype concordance
    vote exists to justify cutting it.

    It rewrites both the FASTA (one record becomes two, `<scaffold>_sub_<start>_<end>`) and
    the name map (one row becomes two, mapping those to chr5_1 and chr9_1), so
    FINALIZE_ASSEMBLY works unchanged.

    THREE MODES
    -----------
      params.chimera_break = false      pass through untouched; candidates already written
                                        by harmonization are the reviewable artifact
      params.chimera_break = '<path>'   break every row in that file -- the file IS the
                                        instruction, so editing it is how you choose
      params.chimera_break = 'auto'     break rows harmonization marked BREAK_CANDIDATE

    `auto` gates on the CONCORDANCE VOTE alone, and that is deliberate rather than a
    limitation accepted reluctantly. The vote is the only signal independent of how the
    scaffold was BUILT. Hi-C cannot be used: measured on the chr5+chr9 fusion, cross-junction
    contact came out at 1.209x matched distance -- ELEVATED, because the junction sits in
    subtelomeric repeat that attracts spurious contacts. That is how YaHS made the error.
    Testing a Hi-C-made join with Hi-C is circular and here it points the wrong way.

    Interstitial telomere signal is a strong independent corroborator -- arrays in BOTH
    orientations at 37.6-43.7 Mb on that scaffold while its own termini have almost none --
    but tidk runs in FINAL_VIZ, downstream of here. CHIMERA_EVIDENCE adds it to the candidates
    file afterwards for the reviewed path. Telomere ABSENCE must never veto a break: a mid-arm
    fusion leaves none.

    WHAT IT WILL DO ON THIS COHORT
    ------------------------------
    Two scaffolds, both in Sde-CTlk_104, both voted 1f/7s:
      hap1  chr5_1+chr9_1   111.6 Mb   junction ~37 Mb
      hap2  chr6_3+chr12_1   76.7 Mb
    Sde-CPla_115's 135 composites are 1-14 Mb and fail the 20 Mb span floor -- correctly, a
    2 Mb chimeric fragment is a fragmented assembly rather than a mis-joined chromosome, and
    breaking it leaves two unplaced halves.

    Input : tuple(meta, fasta, name_map), candidates, script
    Output: tuple(meta, fasta, name_map) -- same shape in as out, so the DAG is unchanged
========================================================================================
*/

process BREAK_CHIMERAS {
    tag "${meta.id}"
    label 'break_chimeras'

    publishDir "${params.outdir}/assembly/chimeras", mode: params.publish_dir_mode,
               pattern: "*.chimera_break_audit.tsv"

    input:
    // `candidates` is now CHIMERA_JOINS's called-joins table, not harmonization's candidates
    // file: it carries cut_bp taken from the AGP and ONE ROW PER CHIMERIC JOIN, so a scaffold
    // with two of them produces three pieces. It also carries candidate_verdict, so auto mode
    // gates without re-deriving the concordance vote.
    tuple val(meta), path(assembly_fasta, stageAs: 'input/*'), path(name_map), path(candidates)
    path(script)

    output:
    tuple val(meta), path("${meta.id}.broken.fasta"), path("${meta.id}.broken_name_map.tsv"),
        emit: assemblies
    tuple val(meta), path("${meta.id}.chimera_break_audit.tsv"), emit: audit
    path("versions.tsv"), emit: versions

    script:
    // 'auto' gates on the verdict harmonization wrote; a path means break every row present
    def mode    = (params.chimera_break?.toString() == 'auto') ? 'auto' : 'file'
    def minpiece = params.chimera_min_piece_bp ?: 1000000
    """
    set -euo pipefail

    # NO_HARMONIZE means this species has a single assembly, so there is no concordance vote
    # and nothing to detect against. Pass through: the chimera arm requires multiple
    # assemblies by construction, the same condition harmonization needs.
    if [ "${name_map.name}" = "NO_HARMONIZE" ] || [ ! -s "${candidates.name}" ]; then
        echo "[BREAK_CHIMERAS ${meta.id}] no name map or no candidates; passing through" >&2
        cp "input/${assembly_fasta.name}" ${meta.id}.broken.fasta
        cp "${name_map.name}" ${meta.id}.broken_name_map.tsv 2>/dev/null \\
            || printf 'old_name\\tnew_name\\torient\\torder\\tlength\\tclass\\tref_span\\tflags\\n' \\
               > ${meta.id}.broken_name_map.tsv
        printf 'metric\\tvalue\\nassembly\\t${meta.id}\\nmode\\tpassthrough\\nscaffolds_broken\\t0\\n' \\
            > ${meta.id}.chimera_break_audit.tsv
        printf 'process\\ttool\\tversion\\n' > versions.tsv
        exit 0
    fi

    python3 ${script} \\
        --fasta "input/${assembly_fasta.name}" \\
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
    IN=\$(grep -v '^>' "input/${assembly_fasta.name}" | tr -d '\\n' | wc -c)
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
    cp "input/${assembly_fasta.name}" ${meta.id}.broken.fasta
    printf 'old_name\\tnew_name\\torient\\torder\\tlength\\tclass\\tref_span\\tflags\\n' \\
      > ${meta.id}.broken_name_map.tsv
    printf 'metric\\tvalue\\nscaffolds_broken\\t0\\n' > ${meta.id}.chimera_break_audit.tsv
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
