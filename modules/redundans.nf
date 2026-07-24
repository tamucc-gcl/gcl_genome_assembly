/*
========================================================================================
    REDUNDANS MODULE
========================================================================================
    Repo location: modules/redundans.nf

    Redundans (Pryszcz & Gabaldón) assists assembly of heterozygous/polymorphic genomes.
    It runs three stages, each independently skippable:
      1. Reduction     — collapse redundant heterozygous contigs (the dedup role)
      2. Scaffolding    — short-read (paired-end) scaffolding via SSPACE
      3. Gap-closing    — fill scaffold gaps with the PE reads

    Used here as the short-read conditioner (the SPAdes analogue of purge_dups on the
    HiFi path). Two usage modes, controlled by params:
      - Short-read ONLY: let redundans finish the assembly itself
        (reduction + scaffolding + gap-closing all on — the defaults).
      - Short-read + downstream scaffolding (Hi-C / long-read / linked-read):
        turn the internal scaffolding/gap-closing off (reduction only) so the finishing
        chain does the scaffolding, OR leave them on to "scaffold the scaffolds".
        Set params.run_redundans_scaffolding / params.run_redundans_gapclosing = false.

    Single collapsed assembly in, single conditioned assembly out (meta.n_hap == 1,
    the 'primary' path). All knobs are params (Elvis / containsKey defaults here; mirror
    in nextflow.config for discoverability). Requires a 'redundans' process label.
========================================================================================
*/

process REDUNDANS {
    tag "${meta.sample}"
    label 'redundans'

    publishDir "${params.outdir}/assembly/contig/redundans", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(assembly_fasta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.sample}.redundans.fasta"), emit: assembly
    path("${meta.sample}_redundans"),                        emit: workdir

    script:
    // --- stage toggles: the user-facing on/off for redundans' internal steps ---
    def do_reduction   = params.containsKey('redundans_run_reduction')   ? params.run_redundans_reduction   : true
    def do_scaffolding = params.containsKey('redundans_run_scaffolding') ? params.run_redundans_scaffolding : true
    def do_gapclosing  = params.containsKey('redundans_run_gapclosing')  ? params.run_redundans_gapclosing  : true

    // --- reduction params ---
    def identity  = params.redundans_identity   ?: 0.51
    def overlap   = params.redundans_overlap    ?: 0.80
    def minlen    = params.redundans_min_contig_bp ?: 200

    // --- scaffolding params ---
    def joins     = params.redundans_joins      ?: 5
    def linkratio = params.redundans_linkratio  ?: 0.7
    def limit     = params.redundans_limit      ?: 0.2
    def mapq      = params.redundans_mapq        ?: 10
    def iters     = params.redundans_iters       ?: 2

    // --- general / minimap2 ---
    def preset    = params.redundans_preset ?: 'asm5'
    def index     = params.redundans_index  ?: '4G'
    def mem_gb    = task.memory ? task.memory.toGiga() : 16

    // --- boolean flags + stage skips assembled into one string ---
    def flags = []
    if (!do_reduction)                                  flags << '--noreduction'
    if (!do_scaffolding)                                flags << '--noscaffolding'
    if (!do_gapclosing)                                 flags << '--nogapclosing'
    if (params.redundans_minimap2reduce    ?: false)    flags << '--minimap2reduce'
    if (params.redundans_usebwa            ?: false)    flags << '--usebwa'
    if (params.redundans_minimap2scaffold  ?: false)    flags << '--minimap2scaffold'
    if (params.redundans_populate_scaffolds?: false)    flags << '--populateScaffolds'
    if (params.redundans_norearrangements  ?: false)    flags << '--norearrangements'
    if (params.redundans_extra)                         flags << params.redundans_extra
    def opt_flags = flags.join(' ')
    """
    # /tmp collides across concurrent SSPACE/GapCloser jobs; \$PWD is on scratch here.
    export TMPDIR="\$PWD"

    redundans.py \\
        --verbose \\
        --fastq ${r1} ${r2} \\
        --fasta ${assembly_fasta} \\
        --outdir ${meta.sample}_redundans \\
        --threads ${task.cpus} \\
        --mem ${mem_gb} \\
        --tmp "\$PWD" \\
        --identity ${identity} \\
        --overlap ${overlap} \\
        --minLength ${minlen} \\
        --joins ${joins} \\
        --linkratio ${linkratio} \\
        --limit ${limit} \\
        --mapq ${mapq} \\
        --iters ${iters} \\
        --index ${index} \\
        --preset ${preset} \\
        ${opt_flags}

    # Grab the most-processed FASTA redundans produced (which one exists depends on
    # which stages ran): all stages -> scaffolds.reduced.fa; no final reduce ->
    # scaffolds.filled.fa; no gap-closing -> scaffolds.fa; reduction only -> contigs.reduced.fa
    outdir=${meta.sample}_redundans
    if   [ -s \$outdir/scaffolds.reduced.fa ]; then final=\$outdir/scaffolds.reduced.fa
    elif [ -s \$outdir/scaffolds.filled.fa  ]; then final=\$outdir/scaffolds.filled.fa
    elif [ -s \$outdir/scaffolds.fa         ]; then final=\$outdir/scaffolds.fa
    elif [ -s \$outdir/contigs.reduced.fa   ]; then final=\$outdir/contigs.reduced.fa
    else echo "ERROR: no redundans output FASTA found in \$outdir" >&2; exit 1
    fi
    cp \$final ${meta.sample}.redundans.fasta
    """

    stub:
    """
    mkdir -p ${meta.sample}_redundans
    touch ${meta.sample}_redundans/scaffolds.reduced.fa
    touch ${meta.sample}.redundans.fasta
    """
}
