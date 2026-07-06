/*
========================================================================================
    SPADES MODULE
========================================================================================
    Short-read (paired-end) de novo assembly with SPAdes.
    Repo location: modules/spades.nf

    Produces a single collapsed assembly (SPAdes does not phase), so a short-read sample
    flows downstream as one assembly (meta.n_hap == 1, the 'primary' path). Organism
    ploidy is tracked separately (meta.ploidy) for genome-size estimation.

    All key knobs are params (Elvis defaults here; mirror in nextflow.config for
    discoverability — see the Phase-4 param list). Based on Jason's 2-spades.sh.

    Requires a 'spades' process label in nextflow.config with a high-memory queue,
    e.g. queue = 'ultramem', clusterOptions = '--qos=highmem', memory = 1000.GB —
    SPAdes' --memory cap is derived from task.memory below.
========================================================================================
*/

process SPADES {
    tag "${meta.sample}"
    label 'spades'

    publishDir "${params.outdir}/assembly/contig/spades", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.sample}.contigs.fasta"),        emit: contigs
    tuple val(meta), path("${meta.sample}.scaffolds.fasta"),      emit: scaffolds, optional: true
    tuple val(meta), path("${meta.sample}.assembly_graph.gfa"),   emit: gfa,       optional: true
    tuple val(meta), path("${meta.sample}.spades.log"),           emit: log

    script:
    def kmers    = params.spades_kmers      ?: '21,33,55,77'
    def cov_opt  = params.spades_cov_cutoff ? "--cov-cutoff ${params.spades_cov_cutoff}" : '--cov-cutoff auto'
    def mode     = params.containsKey('spades_mode') ? params.spades_mode : '--isolate'  // '--isolate' | '--careful' | '--sc' | '' ...
    def extra    = params.spades_extra      ?: ''
    def mem_gb   = task.memory ? task.memory.toGiga() : 1000
    """
    spades.py \\
        -o ${meta.sample}_spades \\
        -k ${kmers} \\
        ${cov_opt} \\
        ${mode} \\
        --pe1-1 ${r1} \\
        --pe1-2 ${r2} \\
        --threads ${task.cpus} \\
        --memory ${mem_gb} \\
        ${extra}

    # Standardize output names (contigs = the assembly output consumed downstream;
    # scaffolds/gfa published for reference). params.spades_output_level selects which
    # FASTA becomes <sample>.contigs.fasta (default 'contigs'; 'scaffolds' uses SPAdes'
    # own PE scaffolding — useful for short-read-only where no Hi-C scaffolding follows).
    cp ${meta.sample}_spades/${(params.spades_output_level ?: 'contigs') == 'scaffolds' ? 'scaffolds' : 'contigs'}.fasta ${meta.sample}.contigs.fasta
    if [ -f ${meta.sample}_spades/scaffolds.fasta ]; then cp ${meta.sample}_spades/scaffolds.fasta ${meta.sample}.scaffolds.fasta; fi
    if [ -f ${meta.sample}_spades/assembly_graph_with_scaffolds.gfa ]; then cp ${meta.sample}_spades/assembly_graph_with_scaffolds.gfa ${meta.sample}.assembly_graph.gfa; fi
    cp ${meta.sample}_spades/spades.log ${meta.sample}.spades.log
    """

    stub:
    """
    touch ${meta.sample}.contigs.fasta
    touch ${meta.sample}.scaffolds.fasta
    touch ${meta.sample}.assembly_graph.gfa
    touch ${meta.sample}.spades.log
    """
}
