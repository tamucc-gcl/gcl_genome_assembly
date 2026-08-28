/*
========================================================================================
    PANGENOME UNTANGLE MODULE
========================================================================================
    Repo location: modules/pangenome_untangle.nf

    `odgi untangle` projects every haplotype path segment into reference coordinate space
    for ONE per-chromosome graph. Output columns:

      query.name query.start query.end ref.name ref.start ref.end score inv self.cov nth.best

    WHY THIS IS THE PRIMARY REARRANGEMENT INSTRUMENT, NOT A SUPPLEMENT
    ------------------------------------------------------------------
    Inversions in a Minigraph-Cactus graph are overwhelmingly NOT representable as bubbles.
    Romain et al. (bioRxiv 2025.03.14.643331) distinguish path-explicit inversions -- where
    ancestral and inverted alleles traverse SHARED nodes in opposite directions, visible to
    topology alone -- from alignment-rescued ones, where the alleles are disjoint paths.
    Measured on Spratelloides:

      path-explicit, from AT traversals   :     29 alleles   (~11 Mb genome-wide)
      alignment-rescued, from revcomp aln :    411 alleles   (2.0 Mb, all <100 kb)
      untangle, chr10 alone               :  ~21.5 Mb inverted, in 5-7 blocks

    One chromosome from untangle exceeds both bubble-based detectors combined by an order of
    magnitude. `self.cov` > 1 is also the ONLY duplication signal available anywhere: AT
    found just 27 node re-visits across 3,268,312 SV alleles, so tandem duplications are
    structurally unrepresentable in this graph's bubbles.

    FULL GRAPH IS PRIMARY
    ---------------------
    Clipping cuts paths into subpaths -- chr10 has 556 paths in the clip graph against 394
    in the full graph, and every clip path carries a [start-end] suffix. A rearrangement
    straddling a subpath boundary is lost, so path projection belongs on the unfragmented
    arm. Clip runs only for comparability with the variant catalog.

    ALL QUERIES, NOT JUST CHROMOSOME-SCALE ONES
    -------------------------------------------
    Measured on chr10: 394 queries took 3m35s against 6m19s for 43. Unrestricted is FASTER,
    because cost is dominated by chromosome-scale paths. And the 346 unplaced scaffolds
    project coherently -- bp/span ratios 0.81-1.00 across 50.7 Mb -- which makes untangle a
    scaffold PLACEMENT instrument as well, at no extra cost.

    PARAMETERS
    ----------
    -e/--cut-every is the resolution knob and cost is flat in it: 6m29s / 6m45s / 6m19s at
    1 Mb / 100 kb / 10 kb on chr10, memory constant. Finer is also more ACCURATE: the union
    of inverted intervals shrinks 13.95 -> 13.74 -> 12.85 Mb as -e tightens, because coarse
    segments absorb non-inverted interstitial sequence, while the run span (gaps <200 kb
    merged) holds at 14.45-14.49 Mb. Hence the 10 kb default. -e is a resolution CEILING,
    not a discovery knob -- max segment size tracks it exactly and total inverted bp is flat.

    -j filters noise (46 of 572 segments scored <0.01 unfiltered). -n > 1 is required for
    dispersed duplications: at the default of 1, nth.best is uniformly 1 by construction.

    SCOPE LIMIT
    -----------
    Cactus splits by refContig, so ref.name has ONE value per chromosome graph and
    INTER-CHROMOSOMAL TRANSLOCATIONS ARE INVISIBLE here. Composite scaffolds do appear
    (e.g. Sde-CPla_115#1#chr10_17+chr11_5), giving a free cross-check against
    harmonization's fusion flags. Translocations need one whole-graph untangle with all
    reference paths as targets -- one bigmem job, not fifteen -- costed separately.

    Input : tuple(taxid, flavor, og, stpidx)
    Output: tuple(taxid, flavor, tsv) / versions
========================================================================================
*/

process PANGENOME_UNTANGLE {
    tag "${taxid}:${flavor}:${og.simpleName}"
    label 'pangenome_untangle'

    // Deliberately NOT 'ignore'. That was set assuming a failure would mean one oversized
    // graph, but the clip-arm failures were a systematic odgi bug hitting all 15 tasks
    // identically -- and 'ignore' let the run continue with an empty clip channel, so
    // PANGENOME_REARRANGE silently never ran and the run "succeeded" with no clip
    // rearrangement data. A missing chromosome must be loud.
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 2

    publishDir "${params.outdir}/pangenome/${taxid}/untangle", mode: params.publish_dir_mode

    input:
    tuple val(taxid), val(flavor), path(og), path(stpidx)

    output:
    tuple val(taxid), val(flavor), path("${og.simpleName}.${flavor}.untangle.tsv.gz"),
        emit: tsv, optional: true
    path("versions.tsv"), emit: versions, optional: true

    script:
    // og is now the COMPACTED graph from PANGENOME_STEPINDEX (<chrom>.opt.og). Nextflow's
    // simpleName strips ALL extensions (documented: file.tar.gz -> file), so this yields
    // <chrom> for chr8_1.og, chr8_1.full.og and chr8_1.opt.og alike, and published names
    // stay <chrom>.<flavor>.untangle.tsv.gz.
    //
    // Deliberately NOT stripping ".opt" defensively here: the output: block above must use
    // the same expression, and it cannot see a script-block variable. If the two ever
    // diverged the task would fail with a missing output -- so a "defensive" transform here
    // would create exactly the breakage it was meant to guard against.
    def base   = og.simpleName
    def cutev  = params.pangenome_untangle_cut_every   ?: 10000
    def minjac = params.pangenome_untangle_min_jaccard ?: 0.1
    def nbest  = params.pangenome_untangle_n_best      ?: 3
    def mdist  = params.pangenome_untangle_merge_dist
    def merge  = mdist ? "-m ${mdist}" : ''
    """
    set -euo pipefail
    export HOME="\$PWD"

    # Target = the reference path. Cactus names the reference WITHOUT a PanSN haplotype
    # field, so it is the path with the fewest '#'-separated fields. Derived rather than
    # hardcoded so this works for any species and any chosen reference.
    #
    # NB do NOT pipe `odgi paths -L` into head: SIGPIPE under `set -o pipefail` truncates
    # the file silently. That produced a 42-line path list from a 556-path graph during
    # development and sent untangle at a subpath as its target.
    odgi paths -i ${og} -L > all_paths.txt
    awk -F'#' '{print NF"\\t"\$0}' all_paths.txt | sort -n | head -n1 | cut -f2 > target.txt

    if [ ! -s target.txt ]; then
        echo "[UNTANGLE ${taxid}:${flavor}:${base}] no paths in graph; nothing to do" >&2
        printf 'process\\ttool\\tversion\\n' > versions.tsv
        exit 0
    fi
    echo "[UNTANGLE ${taxid}:${flavor}:${base}] \$(wc -l < all_paths.txt) paths, target=\$(cat target.txt)" >&2

    # No -Q: every path is a query. Unrestricted measured FASTER than filtering to
    # chromosome-scale paths, and the unplaced-scaffold projections are wanted.
    odgi untangle \\
        -i ${og} \\
        -a ${stpidx} \\
        -R target.txt \\
        -e ${cutev} \\
        -j ${minjac} \\
        -n ${nbest} \\
        ${merge} \\
        -t ${task.cpus} -P \\
        | gzip -c > ${base}.${flavor}.untangle.tsv.gz

    {
      printf 'process\\ttool\\tversion\\n'
      printf '%s\\todgi\\t%s\\n' "${task.process}" "\$(odgi version 2>&1 | head -n1)"
    } > versions.tsv
    """

    stub:
    """
    printf '#query.name\\tquery.start\\tquery.end\\tref.name\\tref.start\\tref.end\\tscore\\tinv\\tself.cov\\tnth.best\\n' \\
      | gzip -c > ${og.simpleName}.${flavor}.untangle.tsv.gz
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
