// The NCBI image has no shell. Run this task on the host and execute its app directly.
process FCS_ADAPTOR {
    tag "${meta.id}"
    label 'fcs_adaptor'
    container false
    conda 'conda-forge::python=3.11'
    publishDir "${params.outdir}/assembly/${stage}/decontam/fcs", mode: params.publish_dir_mode,
        pattern: "${meta.id}.fcs_adaptor*"
    input:
    tuple val(meta), path(assembly_fa), val(mode), val(engine), val(stage)
    path launcher
    output:
    tuple val(meta), path("${meta.id}.fcs_adaptor_report.txt"), emit: action_report
    tuple val(meta), path("${meta.id}.fcsadaptor"), emit: out_dir
    tuple val(meta), path("${meta.id}.fcs_adaptor.log"), emit: log
    tuple val(meta), path("${meta.id}.fcs_adaptor_launch.json"), emit: provenance
    path 'versions.tsv', emit: versions
    script:
    if (!(mode in ['euk','prok'])) error 'FCS-adaptor mode must be euk or prok'
    def runtime = engine == 'nextflow' ? params.decon_container_engine : engine
    def q = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    """
    set -euo pipefail
    mkdir -p ${meta.id}.fcsadaptor
    python3 ${q(launcher)} --fasta ${q(assembly_fa)} --output ${meta.id}.fcsadaptor \
        --image ${q(params.fcs_adaptor_container)} --engine ${q(runtime)} --mode ${mode} \
        --log ${meta.id}.fcs_adaptor.log --provenance ${meta.id}.fcs_adaptor_launch.json
    test -f ${meta.id}.fcsadaptor/fcs_adaptor_report.txt
    cp ${meta.id}.fcsadaptor/fcs_adaptor_report.txt ${meta.id}.fcs_adaptor_report.txt
    printf 'FCS-adaptor-image\\t%s\\n' '${params.fcs_adaptor_container}' > versions.tsv
    """
    stub:
    """
    mkdir -p ${meta.id}.fcsadaptor
    printf '#accession\\tlength\\taction\\trange\\tname\\n' > ${meta.id}.fcs_adaptor_report.txt
    touch ${meta.id}.fcs_adaptor.log versions.tsv
    printf '{}\\n' > ${meta.id}.fcs_adaptor_launch.json
    """
}

process FCS_ADAPTOR_CLEAN {
    tag "${meta.id}"
    label 'fcs_clean'
    publishDir "${params.outdir}/assembly/${stage}/decontam/fcs", mode: params.publish_dir_mode,
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
