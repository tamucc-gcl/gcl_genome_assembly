/*
========================================================================================
    HI-C PAIRS LIFT — HARMONIZED RENAME + REVCOMP (NO REMAPPING)
========================================================================================
    Repo location: modules/hic_lift_harmonized_pairs.nf

    Converts Hi-C pairs.gz from PRE-harmonization scaffold coordinates into FINALIZED
    (harmonized) coordinates, using the harmonized name map. Pure arithmetic; no aligner.

    WHY THIS IS EXACT
    -----------------
    FINALIZE_ASSEMBLY on the harmonized path performs exactly three operations: rename
    old_name -> new_name, reverse-complement where orient == 'rev', and reorder. It does
    NOT filter (the min_scaffold_bp threshold applies only on the NO_HARMONIZE size-rank
    path) and does NOT alter sequence. Scaffold set and lengths are invariant. So:

        chrom  : old_name -> new_name
        pos    : if rev,  pos' = L - pos + 1        (1-based; no off-by-one)
        strand : if rev,  '+' <-> '-'
        order  : restore the upper-triangular convention in the NEW chromosome order

    Contrast with modules/hic_liftover_pairs.nf, which lifts contig -> scaffold coordinates
    through an AGP. That is a harder problem (component offsets, gaps) and needs agptools.
    This is the degenerate case and needs nothing but the name map.

    WHAT THIS DOES NOT HANDLE
    -------------------------
    Anything that changes sequence content or coordinates: gap filling, telomere extension,
    Inspector breaks, decontamination removal. All of those happen UPSTREAM of the assembly
    this lift starts from. Pointing it at pairs mapped to contigs, or to a pre-gap-fill
    assembly, would produce silently wrong coordinates -- so the lifter cross-checks the
    pairs header chromsizes against the map lengths and REFUSES on a mismatch. There is no
    flag to bypass that.

    WIRING (pipeline reorder — see the plan, step 1)
    -----------------------------------------------
    Today MAP_HIC_TO_FINAL aligns to the FINALIZED assembly, which is one step too late for
    HARMONIZE to consume the contacts (it would be a dependency cycle). Retarget it at the
    post-teloclip assembly and lift afterwards:

        MAP_HIC_TO_ASSEMBLY(ch_final_assembly)     // post-teloclip, pre-harmonize
          -> FILTER_HIC_BAM  -> pairs  (post-teloclip coordinates)
             |
             +--> HARMONIZE_SCAFFOLDS  ..... junction evidence, SAME coordinate system as
             |                               the assemblies being harmonized
             +--> HIC_LIFT_HARMONIZED_PAIRS(pairs, name_map)
                    -> CONTACT_MAP_FINAL(lifted pairs, finalized fasta)

    Net Hi-C mapping count is unchanged at one. The final contact maps are equivalent,
    because a rename+revcomp lift is exact. The reordering is also what makes the Hi-C
    junction evidence channel possible at all.

    VALIDATION
    ----------
    Round trip, on already-published files: invert the map, apply it forward, require the
    input back exactly. Reverse-complement is its own inverse, so this exercises the whole
    arithmetic on real data. See py_scripts/validate_pairs_lift.py --roundtrip.

    Once pre-finalize pairs exist, that script also compares them against the lifted output
    directly (Test A, exact) and against a cool built the old way (Test B, tolerance -- the
    two paths align to different references so multi-mapping placement can diverge).

    NAME MAP
    --------
    FINALIZE_ASSEMBLY.out.applied_lift -- <id>.applied_lift.tsv, columns
    old_name new_name orient order length. FINALIZE writes it on BOTH naming paths
    (harmonized, where orient may be 'rev'; and size-rank, where orient is always 'fwd'),
    so there is one schema and this module has no parameters. FINALIZE's 4-column
    name_map.tsv is NOT usable here -- it has no orient column.

    NO PARAMETERS BY DESIGN
    -----------------------
    A length mismatch between the pairs header and the map means the pairs were mapped to a
    different assembly; a scaffold absent from the map means the same. There is no setting
    where continuing is correct, so both are unconditional hard errors.

    Input : tuple(meta, stage, pairs_gz, applied_lift) ; path(lifter)
    Output: pairs / chrom_sizes / stats / versions
========================================================================================
*/

nextflow.enable.dsl = 2

process HIC_LIFT_HARMONIZED_PAIRS {
    tag "${meta.id}_${stage}"
    label 'hic_liftover_pairs'

    publishDir "${params.outdir}/qc/hic_mapping/${stage}/harmonized_lift",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'versions.tsv' ? null : fn }

    input:
    tuple val(meta), val(stage), path(pairs_gz), path(name_map)
    path(lifter)

    output:
    tuple val(meta), val(stage), path("${meta.id}_${stage}.harmonized.pairs.gz"),
          emit: pairs
    tuple val(meta), val(stage), path("${meta.id}_${stage}.harmonized.chrom.sizes"),
          emit: chrom_sizes
    tuple val(meta), val(stage), path("${meta.id}_${stage}.lift_stats.tsv"),
          emit: stats
    path("versions.tsv"), emit: versions

    script:
    """
    set -euo pipefail
    export LC_ALL=C

    # pairtools sort reads TMPDIR; keep it in the task dir (already on /scratch)
    TMPDIR="\${TMPDIR:-\$PWD}"

    # ------------------------------------------------------------------
    # 1) translate coordinates (streaming; O(scaffolds) memory)
    # ------------------------------------------------------------------
    python3 ${lifter} \\
        --pairs ${pairs_gz} \\
        --name-map ${name_map} \\
        --out lifted.unsorted.pairs \\
        --chrom-sizes ${meta.id}_${stage}.harmonized.chrom.sizes \\
        --stats ${meta.id}_${stage}.lift_stats.tsv \\
        --sample-id ${meta.id}

    # ------------------------------------------------------------------
    # 2) re-sort in the new coordinate space
    #    pairtools sort takes the chromosome order from the input pairs HEADER (it does
    #    not accept --chroms-path). The lifter rewrote #chromsize: lines in finalized
    #    order and dropped the stale @SQ samheader records, so the header the sorter
    #    sees describes the body it is sorting. NB modules/hic_liftover_pairs.nf reuses
    #    the ORIGINAL header after renaming chromosomes, which leaves it describing
    #    sequences absent from the body -- worth revisiting there.
    # ------------------------------------------------------------------
    pairtools sort --nproc ${task.cpus} --tmpdir "\${TMPDIR}" lifted.unsorted.pairs \\
      | bgzip -c > ${meta.id}_${stage}.harmonized.pairs.gz

    rm -f lifted.unsorted.pairs

    # ------------------------------------------------------------------
    # 3) sanity gate: every chromosome in the sizes file must appear in the header
    # ------------------------------------------------------------------
    N_SIZES=\$(wc -l < ${meta.id}_${stage}.harmonized.chrom.sizes)
    N_HDR=\$(zcat ${meta.id}_${stage}.harmonized.pairs.gz \\
             | awk '/^#chromsize:/{c++} !/^#/{exit} END{print c+0}')
    if [ "\${N_SIZES}" -ne "\${N_HDR}" ]; then
        echo "[LIFT ${meta.id}] ERROR: chrom.sizes has \${N_SIZES} entries but the sorted" \\
             "pairs header has \${N_HDR} #chromsize lines" >&2
        exit 1
    fi
    echo "[LIFT ${meta.id}] \${N_SIZES} scaffolds; see ${meta.id}_${stage}.lift_stats.tsv"

    # ------------------------------------------------------------------
    # 4) versions (awk drains input -> no SIGPIPE under pipefail; no `head`)
    # ------------------------------------------------------------------
    PT=\$(pairtools --version 2>&1 | awk 'NR==1{print \$NF}')
    BG=\$(bgzip --version 2>&1 | awk 'NR==1{print \$NF}')
    PY=\$(python3 --version 2>&1 | awk '{print \$2}')
    {
        printf 'process\\ttool\\tversion\\n'
        printf '%s\\tpairtools\\t%s\\n' "${task.process}" "\${PT}"
        printf '%s\\tbgzip\\t%s\\n'     "${task.process}" "\${BG}"
        printf '%s\\tpython\\t%s\\n'    "${task.process}" "\${PY}"
    } > versions.tsv
    """

    stub:
    """
    printf '' | bgzip -c > ${meta.id}_${stage}.harmonized.pairs.gz 2>/dev/null \\
        || touch ${meta.id}_${stage}.harmonized.pairs.gz
    touch ${meta.id}_${stage}.harmonized.chrom.sizes
    printf 'metric\\tvalue\\n' > ${meta.id}_${stage}.lift_stats.tsv
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
