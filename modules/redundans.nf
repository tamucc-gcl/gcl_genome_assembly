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
    path "versions.tsv", emit: versions

    script:
    // --- stage toggles: the user-facing on/off for redundans' internal steps ---
    // Boolean knobs can't use the Elvis default (`params.x ?: true` can never yield
    // false), so existence has to be tested separately from the value. boolParam
    // mirrors meta.nf's `pick`: the key is written ONCE, as a String, and read via
    // params[key], so the existence check and the value read cannot drift apart.
    // (They had: the check named 'redundans_run_*' while the read named
    // 'run_redundans_*', so containsKey was always false, all three stages were
    // pinned on, and --noreduction/--noscaffolding/--nogapclosing were unreachable.)
    // Also normalises the String 'false' a CLI `--x false` can deliver.
    def boolParam = { String key, boolean dflt ->
        if (!params.containsKey(key) || params[key] == null) return dflt
        def v = params[key]
        return (v instanceof Boolean) ? v
             : !(v.toString().trim().toLowerCase() in ['false', 'no', '0', 'off'])
    }

    def do_reduction   = boolParam('run_redundans_reduction',   true)
    def do_scaffolding = boolParam('run_redundans_scaffolding', true)
    def do_gapclosing  = boolParam('run_redundans_gapclosing',  true)

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

    # --- stage completion guards -----------------------------------------------------
    # redundans exits 0 even when an internal stage stops early (GapCloser in
    # particular processes what it can and still returns success), and the ladder
    # below selects by file EXISTENCE -- so a partial FASTA gets promoted and
    # published as though the stage had finished. These are conservation invariants,
    # not quality thresholds: each asserts that a stage which was REQUESTED actually
    # did its work. A violation exits non-zero so errorStrategy/maxRetries can act,
    # rather than feeding a truncated assembly to every downstream process.
    nseq() { grep -c '^>' "\$1" || true; }

    if [ "${do_reduction}" = "true" ] && [ ! -s "\$outdir/contigs.reduced.fa" ]; then
        echo "ERROR: reduction requested but \$outdir/contigs.reduced.fa is missing or empty" >&2
        exit 1
    fi

    if [ "${do_scaffolding}" = "true" ] && [ ! -s "\$outdir/scaffolds.fa" ]; then
        echo "ERROR: scaffolding requested but \$outdir/scaffolds.fa is missing or empty" >&2
        exit 1
    fi

    if [ "${do_gapclosing}" = "true" ]; then
        if [ ! -s "\$outdir/scaffolds.filled.fa" ]; then
            echo "ERROR: gap-closing requested but \$outdir/scaffolds.filled.fa is missing or empty" >&2
            exit 1
        fi
        if [ "${do_scaffolding}" = "true" ]; then
            n_in=\$(nseq "\$outdir/scaffolds.fa")
            n_out=\$(nseq "\$outdir/scaffolds.filled.fa")
            if [ "\$n_in" -ne "\$n_out" ]; then
                echo "ERROR: gap-closing did not process every scaffold; output is partial." >&2
                echo "         scaffolds.fa        : \$n_in sequences" >&2
                echo "         scaffolds.filled.fa : \$n_out sequences" >&2
                echo "       Gap-closing fills gaps; it cannot add or drop sequences, so these" >&2
                echo "       counts must be equal. GapCloser returns 0 when it stops early, so" >&2
                echo "       this count is the only evidence that the stage was incomplete." >&2
                exit 1
            fi
        fi
    fi
    if   [ -s \$outdir/scaffolds.reduced.fa ]; then final=\$outdir/scaffolds.reduced.fa
    elif [ -s \$outdir/scaffolds.filled.fa  ]; then final=\$outdir/scaffolds.filled.fa
    elif [ -s \$outdir/scaffolds.fa         ]; then final=\$outdir/scaffolds.fa
    elif [ -s \$outdir/contigs.reduced.fa   ]; then final=\$outdir/contigs.reduced.fa
    else echo "ERROR: no redundans output FASTA found in \$outdir" >&2; exit 1
    fi
    cp \$final ${meta.sample}.redundans.fasta

    printf 'Redundans\t%s\n' "\$(redundans.py --version 2>&1 | head -n1)" > versions.tsv
    """

    stub:
    """
    mkdir -p ${meta.sample}_redundans
    touch ${meta.sample}_redundans/scaffolds.reduced.fa
    touch ${meta.sample}.redundans.fasta
    touch versions.tsv
    """
}
