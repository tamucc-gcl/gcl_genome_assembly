// Nextflow supplies the container; do not launch nested container runners.
process FCS_ADAPTOR {
    tag "${meta.id}"
    label 'fcs_adaptor'
    publishDir "${params.outdir}/${stage}/decontam/fcs", mode: params.publish_dir_mode,
        pattern: "${meta.id}.fcs_adaptor*"
    input:
    tuple val(meta), path(assembly_fa), val(mode), val(engine), val(stage)
    output:
    tuple val(meta), path("${meta.id}.fcs_adaptor_report.txt"), emit: action_report
    tuple val(meta), path("${meta.id}.fcsadaptor"), emit: out_dir
    tuple val(meta), path("${meta.id}.fcs_adaptor.log"), emit: log
    path 'versions.tsv', emit: versions
    script:
    if (!(mode in ['euk','prok'])) error 'FCS-adaptor mode must be euk or prok'
    """
    set -euo pipefail
    mkdir -p ${meta.id}.fcsadaptor
    /app/fcs/bin/av_screen_x -o "\$PWD/${meta.id}.fcsadaptor" --${mode} ${assembly_fa} \
        > ${meta.id}.fcs_adaptor.log 2>&1
    test -f ${meta.id}.fcsadaptor/fcs_adaptor_report.txt
    cp ${meta.id}.fcsadaptor/fcs_adaptor_report.txt ${meta.id}.fcs_adaptor_report.txt
    printf 'FCS-adaptor-image\\t%s\\n' '${params.fcs_adaptor_container}' > versions.tsv
    """
    stub:
    """
    mkdir -p ${meta.id}.fcsadaptor
    printf '#accession\\tlength\\taction\\trange\\tname\\n' > ${meta.id}.fcs_adaptor_report.txt
    touch ${meta.id}.fcs_adaptor.log versions.tsv
    """
}

process FCS_ADAPTOR_CLEAN {
    tag "${meta.id}"
    label 'fcs_clean'
    publishDir "${params.outdir}/${stage}/decontam/fcs", mode: params.publish_dir_mode,
        pattern: "${meta.id}.adaptor*"
    input:
    tuple val(meta), path(assembly_fa), path(report), val(stage)
    output:
    tuple val(meta), path("${meta.id}.adaptor_cleaned.fasta"), emit: cleaned_fasta
    tuple val(meta), path("${meta.id}.adaptor_removed.fasta"), emit: contaminants
    tuple val(meta), path("${meta.id}.adaptor_cleaning.log"), emit: log
    script:
    """
    set -euo pipefail
    /app/bin/gx clean-genome --action-report ${report} --input ${assembly_fa} \
        --output ${meta.id}.adaptor_cleaned.fasta --contam-fasta-out ${meta.id}.adaptor_removed.fasta \
        > ${meta.id}.adaptor_cleaning.log 2>&1
    """
    stub:
    """
    cp ${assembly_fa} ${meta.id}.adaptor_cleaned.fasta
    touch ${meta.id}.adaptor_removed.fasta ${meta.id}.adaptor_cleaning.log
    """
}
