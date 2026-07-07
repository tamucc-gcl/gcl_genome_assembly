/*
========================================================================================
    PILON MODULE
========================================================================================
    Repo location: modules/pilon.nf

    Optional short-read polishing of the conditioned assembly: map the PE reads back with
    bwa mem, then run Pilon to correct base errors / small indels. Gated by
    params.run_pilon (default off) — when off, main.nf passes REDUNDANS output straight
    through. Applied on the short-read branch after REDUNDANS, before FINALIZE_ASSEMBLY.

    Base-polish only here (structural fixes come from scaffolding evidence). Pilon's heap
    scales with genome size, so -Xmx is derived from task.memory. Requires a 'pilon' label.
========================================================================================
*/

process PILON {
    tag "${meta.sample}"
    label 'pilon'

    publishDir "${params.outdir}/assembly/contig/pilon", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(assembly_fasta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.sample}.pilon.fasta"), emit: assembly
    tuple val(meta), path("${meta.sample}.pilon.changes"), emit: changes, optional: true

    script:
    def mem_gb = task.memory ? task.memory.toGiga() : 32
    def fix    = params.pilon_fix ?: 'all'
    def rounds = params.pilon_rounds ?: 1
    def extra  = params.pilon_extra ?: ''
    """
    ref=${assembly_fasta}
    for i in \$(seq 1 ${rounds}); do
        bwa index \$ref
        bwa mem -t ${task.cpus} \$ref ${r1} ${r2} \\
            | samtools sort -@ ${task.cpus} -o aln.\$i.bam -
        samtools index aln.\$i.bam

        pilon -Xmx${mem_gb}g \\
            --genome \$ref \\
            --frags aln.\$i.bam \\
            --output ${meta.sample}.pilon.\$i \\
            --outdir . \\
            --fix ${fix} \\
            --changes \\
            --threads ${task.cpus} ${extra}

        ref=${meta.sample}.pilon.\$i.fasta
    done

    # Final round output -> canonical name
    cp \$ref ${meta.sample}.pilon.fasta
    if [ -f ${meta.sample}.pilon.${rounds}.changes ]; then
        cp ${meta.sample}.pilon.${rounds}.changes ${meta.sample}.pilon.changes
    fi
    """

    stub:
    """
    touch ${meta.sample}.pilon.fasta
    touch ${meta.sample}.pilon.changes
    """
}
