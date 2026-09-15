/*
========================================================================================
    PANGENOME INPUT COVERAGE MODULE
========================================================================================
    Repo location: modules/pangenome_input_coverage.nf

    How much of each input haplotype actually reached the pangenome graph, per chromosome,
    and where a chromosome is missing, why.

    WHY THIS EXISTS
    ---------------
    Both halves of the answer were already being computed and neither was surfaced. Cactus
    writes <outName>.stats.tgz containing sample-stats.tsv -- bp of each haplotype present in
    each reference chromosome's graph. Harmonization writes its report, flagging composite
    scaffolds with a cross-haplotype concordance vote. Nothing joined them, so this went
    unnoticed in every pangenome run built so far:

        chr9_1  Sde-CTlk_104#1  1,281,066 bp   -- 1.9% of the 68 Mb chromosome median

    Sde-CTlk_104_hap1 has no chr9 scaffold. Its chromosome 9 is fused into a 111.6 Mb
    scaffold named chr5_1+chr9_1, and cactus-graphmap-split assigns each contig to a SINGLE
    chromosome -- so the whole scaffold went to chr5's graph and 73.6 Mb of genuine
    chromosome 9 is absent from every graph we have built.

    Nothing else catches it. Cactus does not break chimeric contigs. Inspector works from
    read-to-contig alignment and a Hi-C scaffold join is an N-gap with no reads spanning it,
    so there is nothing to disagree with -- it found two structural errors on this haplotype,
    neither on the fused scaffold.

    THE JOIN IS THE DELIVERABLE
    ---------------------------
    A depressed chromosome count alone says something is wrong. A composite flag alone says a
    scaffold spans two reference chromosomes. Together they say: this haplotype's chrB is
    missing BECAUSE it is fused to chrA, and N other haplotypes keep them separate.

    NO PANGENOME DEPENDENCY BEYOND THE STATS BUNDLE
    -----------------------------------------------
    Deliberate: the break-evidence workflow is two rounds. Round one emits evidence and cuts
    nothing; you review and write the breakpoints file; round two applies it. So the evidence
    must be cheap and must not wait on a 17-hour graph build. This process needs only the
    stats bundle, the harmonization report and the reference .fai.

    THE CHROMOSOME-SCALE FLOOR IS DERIVED, NOT SET
    ----------------------------------------------
    It reuses harmonization's own dropoff rule, applied to the REFERENCE, giving ONE floor for
    the cohort. Applied per assembly it would pick a lower cut in a fragmented assembly --
    Sde-CPla_115_hap1 has 76 composites of 1-14 Mb, some of which would qualify as
    "chromosome-scale for this assembly" while being nothing of the kind.

    Input : tuple(taxid, stats_tgz, harmonization, ref_fai), script
    Output: joined coverage table / break candidates / audit / versions
========================================================================================
*/

process PANGENOME_INPUT_COVERAGE {
    tag "${taxid}"
    label 'pangenome_input_coverage'

    publishDir "${params.outdir}/pangenome/${taxid}/input_coverage",
               mode: params.publish_dir_mode

    input:
    tuple val(taxid), path(stats_tgz), path(harmonization), path(ref_fai)
    path(script)

    output:
    tuple val(taxid), path("${taxid}.input_coverage.tsv"),       emit: coverage
    tuple val(taxid), path("${taxid}.break_candidates.tsv"),     emit: candidates
    tuple val(taxid), path("${taxid}.input_coverage_audit.tsv"), emit: audit
    path("versions.tsv"),                                        emit: versions

    script:
    def lowfrac = params.pangenome_input_coverage_low_frac ?: 0.5
    // harmonization's own defaults, so the derived floor matches what named the composites
    // harmonize_min_scaffold_bp defaults to null, meaning "use finalize_min_scaffold_bp".
    // Chain it the same way harmonization does, or the derived floor is computed from a
    // different starting point than the one that named the composites.
    def minscaf = params.harmonize_min_scaffold_bp ?: params.finalize_min_scaffold_bp ?: 1000000
    def csmeth  = params.harmonize_chromosome_method ?: 'dropoff'
    def dratio  = params.harmonize_dropoff_ratio     ?: 2.0
    def dfrac   = params.harmonize_dropoff_min_frac  ?: 0.5
    def cfrac   = params.harmonize_min_chrom_frac    ?: 0.0
    """
    set -euo pipefail

    # cactus v3.1.4 tars the stats bundle rather than writing a directory, and the paths
    # inside carry the whole scratch prefix, so extract flat and find the file by name.
    mkdir -p stats
    tar xzf ${stats_tgz} -C stats 2>/dev/null || {
        echo "[INPUT_COVERAGE ${taxid}] ERROR: could not extract ${stats_tgz}" >&2; exit 1; }

    SS=\$(find stats -name 'sample-stats.tsv' | head -1)
    if [ -z "\${SS:-}" ]; then
        echo "[INPUT_COVERAGE ${taxid}] ERROR: no sample-stats.tsv in the bundle. Contents:" >&2
        find stats -type f | sed 's/^/  /' >&2
        echo "  Newer cactus writes a <outName>.stats/ DIRECTORY with richer tables" >&2
        echo "  (clipped-by-genome.tsv, refgaps.bed.gz). If this is v3.2.1+, the layout" >&2
        echo "  changed and this process needs updating to prefer those." >&2
        exit 1
    fi
    echo "[INPUT_COVERAGE ${taxid}] using \$SS" >&2

    python3 ${script} \\
        --sample-stats "\$SS" \\
        --harmonization ${harmonization} \\
        --ref-fai ${ref_fai} \\
        --label ${taxid} \\
        --outdir . \\
        --low-frac ${lowfrac} \\
        --min-scaffold-bp ${minscaf} \\
        --chromosome-set-method ${csmeth} \\
        --dropoff-ratio ${dratio} \\
        --dropoff-min-frac ${dfrac} \\
        --min-chrom-frac ${cfrac}

    A=${taxid}.input_coverage_audit.tsv

    # A chromosome far below its cohort median with NO composite to explain it is a
    # DIFFERENT problem -- a genuinely absent chromosome, or a chromosome-assignment failure
    # that left no composite name behind. Surface it, because the explained case is the one
    # we know how to fix and the unexplained case is the one nobody is looking for.
    ue=\$(awk -F'\\t' '\$1=="cells_low_unexplained"{print \$2}' "\$A")
    if [ "\${ue:-0}" -gt 0 ]; then
        echo "[INPUT_COVERAGE ${taxid}] WARNING: \${ue} haplotype x chromosome cell(s) are" >&2
        echo "  far below the cohort median with NO composite explanation:" >&2
        awk -F'\\t' '\$9=="LOW_UNEXPLAINED"{printf "    %s %s  %s bp (%.3f of median)\\n", \$1, \$2, \$4, \$6}' \\
            ${taxid}.input_coverage.tsv >&2
    fi

    nc=\$(awk -F'\\t' '\$1=="break_candidates"{print \$2}' "\$A")
    ct=\$(awk -F'\\t' '\$1=="composites_total"{print \$2}' "\$A")
    echo "[INPUT_COVERAGE ${taxid}] \${nc:-0} break candidate(s) of \${ct:-0} composite(s)" >&2
    if [ "\${nc:-0}" -gt 0 ]; then
        echo "[INPUT_COVERAGE ${taxid}] candidates (nothing has been broken):" >&2
        awk -F'\\t' 'NR>1 && \$1!~/^#/{printf "    %-26s %-24s %12s bp  %s\\n", \$1, \$2, \$3, \$6}' \\
            ${taxid}.break_candidates.tsv >&2
    fi

    {
      printf 'process\\ttool\\tversion\\n'
      printf '%s\\tpython\\t%s\\n' "${task.process}" "\$(python3 --version 2>&1 | awk '{print \$2}')"
    } > versions.tsv
    """

    stub:
    """
    printf 'chromosome\\thaplotype\\tindividual\\tbp_in_graph\\tchrom_median_bp\\tfrac_of_median\\tref_chrom_bp\\tchromosome_scale\\tstatus\\tcomposite\\tcomposite_bp\\tvote\\tflags\\n' \\
      > ${taxid}.input_coverage.tsv
    printf 'assembly\\tcomposite\\tspan_bp\\tn_members\\tmembers_chromosome_scale\\tvote\\tflags\\n' \\
      > ${taxid}.break_candidates.tsv
    printf 'metric\\tvalue\\ncells_low_unexplained\\t0\\nbreak_candidates\\t0\\ncomposites_total\\t0\\n' \\
      > ${taxid}.input_coverage_audit.tsv
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
