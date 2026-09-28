/* Assess current pre-finishing joins using telomeres, N gaps, and existing Hi-C mappings lifted through the last-round AGP only. Hi-C is supplementary evidence, not independent proof of a scaffolding error. */

process CHIMERA_EVIDENCE {
    tag "${asm_id}"
    label 'chimera_evidence'

    publishDir "${params.outdir}/assembly/chimeras/evidence", mode: params.publish_dir_mode

    input:
    // No stageAs: every output here is ${asm_id}.*.chimera_evidence.*, so the input FASTA
    // cannot collide with one. (break_chimeras.nf needs it because its output IS a FASTA.)
    // With `stageAs: 'asm/*'` the .name property already carries the prefix, so prefixing it
    // again by hand produced a doubled path and samtools faidx failed.
    tuple val(taxid), val(asm_id), path(assembly_fasta), path(called_joins),
          path(round1_agp), path(round2_agp), path(contig_pairs), val(telomere_motif)
    path(hic_script)
    path(evidence_script)

    output:
    tuple val(taxid), val(asm_id), path("${asm_id}.*.chimera_evidence.tsv"),
        emit: evidence, optional: true
    tuple val(taxid), val(asm_id), path("${asm_id}.*.chimera_evidence.png"),
        emit: figures, optional: true
    path("${asm_id}.hic_pairs_audit.tsv"), emit: pairs_audit, optional: true
    path("versions.tsv"),                  emit: versions

    script:
    def res    = params.chimera_hic_resolution ?: 100000
    def winbin = params.chimera_hic_window_bins ?: 20
    def gapwin = params.chimera_gap_snap_window ?: 500000
    def motif  = telomere_motif ?: 'CCCTAA'
    """
    set -euo pipefail

    # Only CALLABLE joins on a BREAK_CANDIDATE scaffold are worth confirming. A REVIEW row is
    # by definition one the vote could not resolve, and a NOT_A_CANDIDATE row was excluded on
    # span or on the vote -- neither is about to be cut, so neither needs a contact map built
    # for it. Sde-CPla_115_hap1 has ~79 callable joins and none is a break candidate.
    awk -F'\\t' 'NR>1 && \$1!~/^#/ && \$15=="yes" && \$19=="BREAK_CANDIDATE" {print \$2}' \\
        ${called_joins} | sort -u > scaffolds.txt

    if [ ! -s scaffolds.txt ]; then
        echo "[CHIMERA_EVIDENCE ${asm_id}] no callable break candidates; nothing to confirm" >&2
        printf 'process\\ttool\\tversion\\n' > versions.tsv
        exit 0
    fi
    echo "[CHIMERA_EVIDENCE ${asm_id}] confirming \$(wc -l < scaffolds.txt) scaffold(s)" >&2

    R2=""
    if [ "${round2_agp.name}" != "NO_ROUND2" ] && [ -s "${round2_agp.name}" ]; then
        R2="--round2 ${round2_agp}"
    fi

    # ---- one pairs pass for ALL this assembly's candidates -------------------------
    # The input is ~1.9 GB gzipped, so filtering per candidate would decompress it once per
    # scaffold. chimera_hic_pairs.py takes the whole list and writes one file each.
    HAVE_HIC=0
    if [ "${contig_pairs.name}" != "NO_PAIRS" ] && [ -s "${contig_pairs.name}" ]; then
        python3 ${hic_script} \\
            --pairs ${contig_pairs} \\
            --round1 ${round1_agp} \\
            \$R2 \\
            --scaffolds "\$(paste -sd, scaffolds.txt)" \\
            --outdir . \\
            --label ${asm_id}
        HAVE_HIC=1
    else
        echo "[CHIMERA_EVIDENCE ${asm_id}] no contig-space pairs: telomere and N-gap" >&2
        echo "  evidence only. Those are independent of Hi-C, which confirms rather than" >&2
        echo "  justifies, so their absence weakens nothing about the decision." >&2
    fi

    samtools faidx ${assembly_fasta}

    while read -r SC; do
        [ -n "\$SC" ] || continue
        SAFE=\$(echo "\$SC" | tr '+' '_')

        # ---- mini reference: this scaffold alone ----------------------------------
        samtools faidx ${assembly_fasta} "\$SC" > "\${SAFE}.mini.fa"
        if [ ! -s "\${SAFE}.mini.fa" ]; then
            echo "[CHIMERA_EVIDENCE ${asm_id}] WARNING \$SC not in the assembly; skipping" >&2
            continue
        fi

        # ---- telomeres: tidk is the source of record -------------------------------
        # Not a hand-rolled motif count. tidk is already a pipeline dependency and its output
        # is what FINAL_VIZ reports, so computing the same quantity a second way would leave
        # two numbers and no way to adjudicate a disagreement. It also normalises the
        # canonical repeat. Note tidk reports the window END; the reader converts to START.
        TW=""
        if tidk search --string '${motif}' --output "\${SAFE}" --dir . --extension tsv \\
               "\${SAFE}.mini.fa" > "\${SAFE}.tidk.log" 2>&1; then
            for c in "\${SAFE}_telomeric_repeat_windows.tsv" "\${SAFE}.tsv"; do
                [ -s "\$c" ] && TW="--telomere-windows \$c" && break
            done
        fi
        if [ -z "\$TW" ]; then
            echo "[CHIMERA_EVIDENCE ${asm_id}] tidk produced no windows for \$SC;" >&2
            echo "  falling back to a sequence-derived motif count (recorded in the audit)" >&2
        fi

        # ---- Hi-C: translate, RE-SORT, load ---------------------------------------
        COOL=""
        if [ "\$HAVE_HIC" = "1" ] && [ -s "${asm_id}.\${SC}.pairs" ]; then
            # translation destroys the input's contig-order sort, and cooler cload requires
            # sorted input. Small file by now: one scaffold's pairs only.
            sort -k2,2 -k4,4 -k3,3n -k5,5n "${asm_id}.\${SC}.pairs" > "\${SAFE}.sorted.pairs"
            if cooler cload pairs \\
                   -c1 2 -p1 3 -c2 4 -p2 5 \\
                   "${asm_id}.\${SC}.chromsizes:${res}" \\
                   "\${SAFE}.sorted.pairs" "\${SAFE}.cool" >"\${SAFE}.cload.log" 2>&1; then
                COOL="--cool \${SAFE}.cool"
            else
                echo "[CHIMERA_EVIDENCE ${asm_id}] cooler cload failed for \$SC:" >&2
                tail -5 "\${SAFE}.cload.log" | sed 's/^/    /' >&2
            fi
        fi

        # ---- one evidence record per CUT on this scaffold ---------------------------
        awk -F'\\t' -v S="\$SC" \\
            'NR>1 && \$1!~/^#/ && \$2==S && \$15=="yes" && \$19=="BREAK_CANDIDATE" {print \$4}' \\
            ${called_joins} | sort -un > "\${SAFE}.cuts.txt"

        while read -r CUT; do
            [ -n "\$CUT" ] || continue
            python3 ${evidence_script} \\
                --fasta "\${SAFE}.mini.fa" \\
                --candidates ${called_joins} \\
                --assembly ${asm_id} \\
                --scaffold "\$SC" \\
                --cut-bp "\$CUT" \\
                --telomere-motif '${motif}' \\
                \$TW \$COOL \\
                --outdir . \\
                --label "${asm_id}.\${SAFE}_\${CUT}" \\
                --win-bins ${winbin} \\
                --gap-snap-window ${gapwin}
        done < "\${SAFE}.cuts.txt"
    done < scaffolds.txt

    # What the confirmation actually said, in the task log: a reviewer should not have to open
    # a file to learn whether the evidence supported the cut.
    for f in ${asm_id}.*.chimera_evidence.tsv; do
        [ -s "\$f" ] || continue
        awk -F'\\t' -v F="\$(basename \$f)" '
            \$1=="cut_bp"{c=\$2} \$1=="verdict"{v=\$2} \$1=="notes"{n=\$2}
            END{printf "[CHIMERA_EVIDENCE] %s cut=%s %s [%s]\\n", F, c, v, n}' "\$f" >&2
    done

    {
      printf 'process\\ttool\\tversion\\n'
      printf '%s\\tsamtools\\t%s\\n' "${task.process}" "\$(samtools --version 2>&1 | awk 'NR==1{print \$2}')"
      printf '%s\\ttidk\\t%s\\n'     "${task.process}" "\$(tidk --version 2>&1 | awk '{print \$NF}')"
      printf '%s\\tcooler\\t%s\\n'   "${task.process}" "\$(cooler --version 2>&1 | awk '{print \$NF}')"
      printf '%s\\tpython\\t%s\\n'   "${task.process}" "\$(python3 --version 2>&1 | awk '{print \$2}')"
    } > versions.tsv
    """

    stub:
    """
    printf 'metric\\tvalue\\ncut_bp\\t0\\nverdict\\tSTUB\\n' \\
      > ${asm_id}.stub.chimera_evidence.tsv
    printf 'metric\\tvalue\\npairs_written\\t0\\n' > ${asm_id}.hic_pairs_audit.tsv
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
