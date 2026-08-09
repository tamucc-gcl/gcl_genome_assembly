/*
========================================================================================
    CACTUS PANGENOME MODULE
========================================================================================
    Repo location: modules/cactus_pangenome.nf

    One minigraph-cactus graph per species. Inputs are already PanSN-named by the
    PANGENOME subworkflow: `names` are the seqfile sample names (reference has NO
    haplotype suffix; other haplotypes are <sample>.<hap>), aligned by position with the
    staged `fastas`. This module builds the two-column seqfile, derives --refContigs from
    the reference's harmonized chromosome scaffolds (cactus auto-detection only matches
    'chr'+<=3 chars, so it misses chr10_1..; we set them explicitly), and runs cactus.

    Toil runs single-machine inside the container. The jobstore must not pre-exist and the
    workDir must be on fast scratch -- both live in the task dir, which this pipeline
    already places on /scratch (process.scratch), and are discarded after the task.

    Input : tuple(taxid, ref_name, names, fastas)
    Output: gbz / gfa / vcf / odgi(.og) / full graph dir  + versions
========================================================================================
*/

process CACTUS_PANGENOME {
    tag "taxid_${taxid}"
    label 'cactus_pangenome'

    // publish flat: cactus nests everything under out/; strip that prefix at publish time
    // (leaves the workdir untouched, so the graph itself is never moved/rewritten).
    publishDir "${params.outdir}/pangenome/${taxid}", mode: params.publish_dir_mode,
        saveAs: { fn -> fn.startsWith('out/') ? fn.substring(4) : fn }

    input:
    tuple val(taxid), val(ref_name), val(names), path(fastas)

    output:
    // primary clip-graph handles (exact names -> single files, never the .full.* variants)
    tuple val(taxid), path("out/${taxid}.gbz"),        emit: gbz
    tuple val(taxid), path("out/${taxid}.og"),         emit: og
    tuple val(taxid), path("out/${taxid}.hapl"),       emit: hapl
    tuple val(taxid), path("out/${taxid}.snarls"),     emit: snarls
    tuple val(taxid), path("out/${taxid}.gfa.gz"),     emit: gfa
    tuple val(taxid), path("out/${taxid}.vcf.gz"),     emit: vcf
    tuple val(taxid), path("out/${taxid}.vcf.gz.tbi"), emit: vcf_tbi
    tuple val(taxid), path("out/${taxid}.raw.vcf.gz"), emit: raw_vcf
    tuple val(taxid), path("out/${taxid}.chroms/*"),   emit: chrom_og
    tuple val(taxid), path("out/${taxid}.gaf.gz"),      emit: gaf, optional: true
    tuple val(taxid), path("out/${taxid}.viz/*"),      emit: viz
    // catch-all: keeps + publishes everything cactus produced EXCEPT the construction scratch
    // removed in-script (full graphs, HAL, GAF/PAF, SV graph, raw VCF, stats bundle, ...).
    // '**' so future cactus outputs are retained automatically (general-purpose).
    tuple val(taxid), path("out/**"),                  emit: all
    tuple val(taxid), path("seqfile.txt"),             emit: seqfile
    path("versions.tsv"),                              emit: versions

    script:
    def extra = params.pangenome_cactus_extra ?: ''
    // -gpu image runs KegAlign automatically; --gpu 1 pins it to the single requested GPU
    // and --lastzMemory is the recommended cluster safeguard for the alignment jobs.
    def gpu   = params.pangenome_use_gpu ? '--gpu 1 --lastzMemory 100G' : ''
    """
    set -euo pipefail

    # cactus/Toil write config under ~/.toil and heavy scratch under TMPDIR. The container's
    # inherited HOME (/home/<user>) and /tmp are not writable/roomy on the compute node, so
    # point both at the Nextflow task dir (on /scratch): writable, discarded after the task.
    export HOME="\$PWD"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$TMPDIR"

    ids=(${names.join(' ')})
    fastas=(${fastas})

    # ---- two-column seqfile: <SAMPLE.HAP>  <fasta> ; locate the reference fasta ----
    : > seqfile.txt
    REF_FA=""
    for i in "\${!ids[@]}"; do
        printf '%s\\t%s\\n' "\${ids[\$i]}" "\${fastas[\$i]}" >> seqfile.txt
        [ "\${ids[\$i]}" = "${ref_name}" ] && REF_FA="\${fastas[\$i]}"
    done
    if [ -z "\${REF_FA}" ]; then
        echo "[PANGENOME ${taxid}] ERROR: reference '${ref_name}' not among inputs" >&2; exit 1
    fi

    # ---- refContigs = reference's harmonized chromosome scaffolds (chrN_p, not composites) ----
    REFCONTIGS=\$(grep '^>' "\${REF_FA}" | sed 's/^>//; s/[[:space:]].*//' \\
        | awk '/^chr[0-9]+_[0-9]+\$/' | tr '\\n' ' ')
    if [ -z "\${REFCONTIGS}" ]; then
        echo "[PANGENOME ${taxid}] ERROR: no chr-named scaffolds in reference \${REF_FA}" >&2; exit 1
    fi
    echo "[PANGENOME ${taxid}] reference=${ref_name}; refContigs=\${REFCONTIGS}"

    # ---- run cactus (jobstore must not exist; workDir + jobstore on scratch task dir) ----
    rm -rf js cactus_work out
    mkdir -p cactus_work out

    cactus-pangenome \\
        ./js \\
        seqfile.txt \\
        --workDir ./cactus_work \\
        --outDir ./out \\
        --outName ${taxid} \\
        --reference ${ref_name} \\
        --refContigs \${REFCONTIGS} \\
        --vcf \\
        --haplo \\
        --gfa full clip \\
        --gbz full clip \\
        --viz full clip \\
        --odgi full clip \\
        --chrom-og full clip \\
        --maxCores ${task.cpus} \\
        ${gpu} \\
        ${extra}

    # ---- cull construction scratch (pre-join per-chromosome intermediates, superseded by
    # the joined graph; the bulk of the file count) and the duplicate seqfile. Everything
    # else cactus produced is kept and published (flattened) via the output block above.
    rm -rf out/chrom-subproblems out/chrom-alignments
    rm -f  out/seqfile.txt

    CV=\$(cactus --version 2>&1 | awk 'NR==1{print}')
    printf 'process\\ttool\\tversion\\n%s\\tcactus\\t%s\\n' "${task.process}" "\${CV}" > versions.tsv
    """

    stub:
    """
    mkdir -p out out/${taxid}.chroms out/${taxid}.viz
    for x in gbz og hapl snarls gfa.gz vcf.gz vcf.gz.tbi; do : > "out/${taxid}.\$x"; done
    : > out/${taxid}.chroms/chr1_1.og
    : > out/${taxid}.viz/chr1_1.viz.png
    : > seqfile.txt
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
