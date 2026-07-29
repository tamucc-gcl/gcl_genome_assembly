/*
========================================================================================
    GETORGANELLE MODULE
========================================================================================
    Repo location: modules/getorganelle.nf

    Short-read organelle assembly (the short-read analog of MITOHIFI). One task per
    (sample x organelle type). Implements the contract we validated on real data:

      1. ATTEMPT   get_organelle_from_reads.py, coverage-capped (--reduce-reads-for-coverage).
                   Capping organelle coverage is the decisive lever: deep WGS over-covers
                   organelles (30,000x observed) and yields multi-Mb error hairballs; ~50-100x
                   gives a clean plastome/mito-sized graph in minutes.
      2. PRUNE     if no circle was closed, re-disentangle the graph with
                   get_organelle_from_assembly.py --min-depth to drop low-depth contamination
                   (rescues contamination-limited samples; harmless otherwise). Toggle via
                   params.getorganelle_from_assembly.
      3. SELECT    prefer a certified-circular config; else take the largest scaffolded
                   (linear) path. Always emit a usable FASTA (unless nothing assembled) plus
                   status in {circular, linear, failed}. Plant mitogenomes and IR/MTPT-tangled
                   plastomes legitimately end at 'linear' — a valid result, not a failure —
                   and it feeds annotation + organelle-removal unchanged.

    Input:
      tuple(meta, sr1, sr2, org)   org = [type, recursion, kmers, coverage]  (one target)
      val config_dir               GetOrganelle DB dir (DOWNLOAD_GETORGANELLE_DB)
    Output:
      assembly  tuple(meta, org.type, fasta)   canonical FASTA, sample/status-tagged headers
                                               (absent if status=failed -> sample drops out)
      graph     tuple(meta, org.type, gfa)     largest selected graph (Bandage/QC)
      stats     tuple(meta, org.type, tsv)     status/length/params (ALWAYS, incl. failed)
      log       the GetOrganelle reads-run log
========================================================================================
*/

process GETORGANELLE {
    tag "${meta.sample}:${org.type}"
    label 'getorganelle'

    publishDir "${params.outdir}/organelle", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(sr1), path(sr2), val(org)
    val config_dir

    output:
    tuple val(meta), val(org.type), path("${meta.sample}.${org.type}.fasta"),           emit: assembly, optional: true
    tuple val(meta), val(org.type), path("${meta.sample}.${org.type}.gfa"),             emit: graph,    optional: true
    tuple val(meta), val(org.type), path("${meta.sample}.${org.type}.stats.tsv"),       emit: stats
    tuple val(meta), val(org.type), path("${meta.sample}.${org.type}.get_org.log.txt"), emit: log,      optional: true
    path "versions.tsv", emit: versions

    script:
    def cov       = org.coverage
    def rounds    = org.recursion
    def kmers     = org.kmers
    def wsize     = org.word_size ?: params.getorganelle_word_size
    def wflag     = wsize ? "-w ${wsize}" : ''
    def extra     = params.getorganelle_extra_args ?: ''
    def do_prune  = (params.getorganelle_from_assembly == false) ? 'false' : 'true'
    def min_depth = params.getorganelle_min_depth ?: 10
    """
    set -eu
    shopt -s nullglob

    SAMPLE="${meta.sample}"
    ORG="${org.type}"
    echo "[GETORGANELLE] sample=\${SAMPLE} organelle=\${ORG} cov=${cov} -R=${rounds} -k=${kmers} ${wflag}"

    # ---------- 1) primary attempt: recruit + assemble from reads (coverage-capped) ----------
    get_organelle_from_reads.py \\
        --config-dir "${config_dir}" \\
        -1 ${sr1} -2 ${sr2} \\
        -F \${ORG} \\
        --reduce-reads-for-coverage ${cov} \\
        ${wflag} \\
        -R ${rounds} \\
        -k ${kmers} \\
        -o reads_out \\
        -t ${task.cpus} \\
        ${extra} || echo "[GETORGANELLE] reads run exited non-zero; will inspect the graph"

    [ -f reads_out/get_org.log.txt ] && cp reads_out/get_org.log.txt "\${SAMPLE}.\${ORG}.get_org.log.txt" || true

    # labelled graph from the reads run (for the prune-retry and QC)
    rg=(reads_out/extended_spades/K*/*extend-embplant*.fastg)
    READS_GRAPH=""
    [ \${#rg[@]} -gt 0 ] && READS_GRAPH=\$(ls -S "\${rg[@]}" | head -n1)

    # ---------- 2) optional prune + retry if no circle was closed ----------
    comp0=(reads_out/*complete*.path_sequence.fasta)
    if [ \${#comp0[@]} -eq 0 ] && [ "${do_prune}" = "true" ] && [ -n "\${READS_GRAPH}" ]; then
        echo "[GETORGANELLE] no circle from reads; pruning graph (--min-depth ${min_depth}) and retrying"
        get_organelle_from_assembly.py \\
            --config-dir "${config_dir}" \\
            -F \${ORG} \\
            -g "\${READS_GRAPH}" \\
            --min-depth ${min_depth} \\
            -o fromasm_out \\
            -t ${task.cpus} || echo "[GETORGANELLE] from-assembly retry exited non-zero; inspecting"
    fi

    # ---------- 3) select best: circle > linear; prefer reads_out over fromasm_out ----------
    STATUS=failed; PICK=""; GRAPHDIR=""
    for d in reads_out fromasm_out; do
        [ -d "\$d" ] || continue
        comp=("\$d"/*complete*.path_sequence.fasta)
        if [ \${#comp[@]} -gt 0 ]; then
            PICK=\$(ls -S "\${comp[@]}" | head -n1); STATUS=circular; GRAPHDIR="\$d"; break
        fi
    done
    if [ -z "\${PICK}" ]; then
        for d in reads_out fromasm_out; do
            [ -d "\$d" ] || continue
            scaf=("\$d"/*scaffolds*.path_sequence.fasta)
            anyp=("\$d"/*.path_sequence.fasta)
            if [ \${#scaf[@]} -gt 0 ]; then
                PICK=\$(ls -S "\${scaf[@]}" | head -n1); STATUS=linear; GRAPHDIR="\$d"; break
            elif [ \${#anyp[@]} -gt 0 ]; then
                PICK=\$(ls -S "\${anyp[@]}" | head -n1); STATUS=linear; GRAPHDIR="\$d"; break
            fi
        done
    fi

    # ---------- 4) write canonical outputs ----------
    CANON="\${SAMPLE}.\${ORG}.fasta"
    n=0
    if [ -n "\${PICK}" ]; then
        while IFS= read -r line; do
            case "\${line}" in
                ">"*) n=\$((n+1)); printf '>%s.%s.%s.%d %s\\n' "\${SAMPLE}" "\${ORG}" "\${STATUS}" "\${n}" "\${line#>}" ;;
                *)    printf '%s\\n' "\${line}" ;;
            esac
        done < "\${PICK}" > "\${CANON}"
    fi

    if [ -n "\${GRAPHDIR}" ]; then
        gfas=("\${GRAPHDIR}"/*.selected_graph.gfa)
        [ \${#gfas[@]} -gt 0 ] && cp "\$(ls -S "\${gfas[@]}" | head -n1)" "\${SAMPLE}.\${ORG}.gfa"
    fi

    LEN=0
    [ -s "\${CANON}" ] && LEN=\$(grep -v '^>' "\${CANON}" | tr -d '\\n' | wc -c)

    printf 'sample\\torganelle\\tstatus\\tn_seqs\\ttotal_len\\tcoverage_target\\trecursion\\tsource\\n' > "\${SAMPLE}.\${ORG}.stats.tsv"
    printf '%s\\t%s\\t%s\\t%d\\t%d\\t%s\\t%s\\t%s\\n' "\${SAMPLE}" "\${ORG}" "\${STATUS}" "\${n}" "\${LEN}" "${cov}" "${rounds}" "\${PICK:-none}" >> "\${SAMPLE}.\${ORG}.stats.tsv"
    echo "[GETORGANELLE] \${SAMPLE} \${ORG}: status=\${STATUS} n_seqs=\${n} len=\${LEN}"

    # drop an empty canonical so 'assembly' is simply absent on failure (sample drops out)
    [ -s "\${CANON}" ] || rm -f "\${CANON}"

    printf 'GetOrganelle\\t%s\\n' "\$(get_organelle_from_reads.py --version 2>&1 | sed -n 's/^GetOrganelle //p' | head -n1)" > versions.tsv
    """

    stub:
    """
    echo ">${meta.sample}.${org.type}.linear.1 stub" > ${meta.sample}.${org.type}.fasta
    echo "ACGTACGTACGTACGT" >> ${meta.sample}.${org.type}.fasta
    touch ${meta.sample}.${org.type}.gfa
    printf 'sample\\torganelle\\tstatus\\tn_seqs\\ttotal_len\\tcoverage_target\\trecursion\\tsource\\n%s\\t%s\\tlinear\\t1\\t16\\t${org.coverage}\\t${org.recursion}\\tstub\\n' "${meta.sample}" "${org.type}" > ${meta.sample}.${org.type}.stats.tsv
    echo stub > ${meta.sample}.${org.type}.get_org.log.txt
    touch versions.tsv
    """
}
