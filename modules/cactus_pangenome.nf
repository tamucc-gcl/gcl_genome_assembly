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

    publishDir "${params.outdir}/pangenome/${taxid}", mode: params.publish_dir_mode

    input:
    tuple val(taxid), val(ref_name), val(names), path(fastas)

    output:
    tuple val(taxid), path("out/*.gbz"),    emit: gbz
    tuple val(taxid), path("out/*.og"),     emit: og,   optional: true
    tuple val(taxid), path("out/*.gfa.gz"), emit: gfa,  optional: true
    tuple val(taxid), path("out/*.vcf.gz"), emit: vcf,  optional: true
    tuple val(taxid), path("seqfile.txt"),  emit: seqfile
    path("out/**"),                         emit: all
    path("versions.tsv"),                   emit: versions

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

    CV=\$(cactus --version 2>&1 | awk 'NR==1{print}')
    printf 'process\\ttool\\tversion\\n%s\\tcactus\\t%s\\n' "${task.process}" "\${CV}" > versions.tsv
    """

    stub:
    """
    mkdir -p out
    : > out/${taxid}.gbz
    : > seqfile.txt
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
