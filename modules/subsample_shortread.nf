process SUBSAMPLE_SHORTREAD {
    tag "${meta.sample}"
    label 'subsample_shortreads'
    publishDir "${params.outdir}/fastq/shortread/subsample", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(r1), path(r2), path(gsize)

    output:
    tuple val(meta), path("${meta.sample}.sub.R1.fastq.gz"),
                     path("${meta.sample}.sub.R2.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.sample}.subsample.log"),   emit: log
    path "versions.tsv",                                     emit: versions

    script:
    def target = params.shortread_target_depth ?: 0
    def seed   = params.shortread_subsample_seed ?: 100
    """
    G=\$(tr -dc '0-9' < ${gsize})

    if [ "${target}" -le 0 ] || [ -z "\$G" ]; then
        # disabled or genome size unknown (GenomeScope emitted NA) -> pass through
        ln -s \$(readlink -f ${r1}) ${meta.sample}.sub.R1.fastq.gz
        ln -s \$(readlink -f ${r2}) ${meta.sample}.sub.R2.fastq.gz
        echo "pass-through (target=${target}, genome_size=\${G:-NA})" > ${meta.sample}.subsample.log
    else
        rasusa reads -c ${target} -g \${G} -s ${seed} \\
            -o ${meta.sample}.sub.R1.fastq.gz \\
            -o ${meta.sample}.sub.R2.fastq.gz \\
            ${r1} ${r2} 2> ${meta.sample}.subsample.log
    fi

    printf 'rasusa\t%s\n' "\$(rasusa --version 2>&1 | awk '{print \$NF}')" > versions.tsv
    """

    stub:
    """
    touch ${meta.sample}.sub.R1.fastq.gz ${meta.sample}.sub.R2.fastq.gz
    touch ${meta.sample}.subsample.log versions.tsv
    """
}