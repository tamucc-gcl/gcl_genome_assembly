process COUNT_SHORTREAD {
    tag "${meta.sample}"
    label 'read_count' 

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path(r1), path(r2), path("${meta.sample}.bases.txt"), emit: counted
    path "versions.tsv", emit: versions

    script:
    """
    # sum_len is column 5 of `seqkit stats -T`; sum across both mates -> total input bases
    seqkit stats -T -j ${task.cpus} ${r1} ${r2} \\
        | awk 'NR>1 { s += \$5 } END { print s+0 }' > ${meta.sample}.bases.txt

    printf 'seqkit\t%s\n' "\$(seqkit version 2>&1 | awk '{print \$NF}')" > versions.tsv
    """

    stub:
    """
    echo 0 > ${meta.sample}.bases.txt
    touch versions.tsv
    """
}