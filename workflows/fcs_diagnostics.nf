include { MAPPING_QC as FCS_MAPPING } from '../modules/mapping_qc.nf'

process FCS_PLOT_EVIDENCE {
    tag "${meta.id}"
    conda 'conda-forge::python=3.11'
    cpus 1
    memory '2 GB'
    time '2h'
    publishDir "${params.outdir}/${stage}/decontam/fcs", mode:'copy',pattern:"${meta.id}.fcs_metrics.*"
    input:
    tuple val(meta), path(fasta), path(taxonomy), path(mapping)
    val stage
    path helpers
    output:
    tuple val(meta.id), path("${meta.id}.fcs_metrics.tsv"), emit: metrics
    tuple val(meta.id), path("${meta.id}.fcs_metrics.json"), emit: provenance
    script:
    def q={ x -> "'"+x.toString().replace("'","'\"'\"'")+"'" }
    """
    python3 fcs_plot_evidence.py --fasta ${q(fasta)} --taxonomy ${q(taxonomy)} \
        --coverage ${q(mapping.resolve('coverage.txt'))} --assembly ${q(meta.id)} \
        --stage ${q(stage)} --out ${q(meta.id+'.fcs_metrics.tsv')}
    """
}

workflow FCS_DIAGNOSTICS {
    take:
    assemblies
    taxonomy
    reads
    stage
    shared_mapping
    reuse_mapping
    main:
    def ch_mapping = shared_mapping
    def ch_versions = Channel.empty()
    if(!reuse_mapping) {
        FCS_MAPPING(assemblies.map { m,fa -> tuple(m.sample,m,fa) }
            .combine(reads.map { m,r -> tuple(m.sample,[r]) },by:0)
            .map { sample,m,fa,wrapped -> tuple(m,fa,wrapped[0]) },false)
        ch_mapping = FCS_MAPPING.out.results
        ch_versions = FCS_MAPPING.out.versions
    }
    ch_plot_inputs = assemblies.map { m,fa -> tuple(m.id.toString(),m,fa) }
        .join(taxonomy.map { m,p -> tuple(m.id.toString(),p) },failOnDuplicate:true,failOnMismatch:true)
        .join(ch_mapping.map { m,p -> tuple(m.id.toString(),p) },failOnDuplicate:true,failOnMismatch:true)
        .map { id,m,fa,tax,map -> tuple(m,fa,tax,map) }
    FCS_PLOT_EVIDENCE(ch_plot_inputs,stage,['fcs_plot_evidence.py','contamination_review.py'].collect {
        file("${moduleDir}/../py_scripts/${it}",checkIfExists:true) })
    emit:
    metrics = FCS_PLOT_EVIDENCE.out.metrics
    versions = ch_versions
    mapping = ch_mapping
}
