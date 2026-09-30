/* One same-species graph from an explicit eligible cohort and reference chromosome list.
 * CLIP supplies biological products. GREF(CLIP) supplies standard variants.
 * FULL GFA plus Cactus clipping statistics are support artifacts, with no parallel analysis.
 * Graph-role validation runs separately so a missing optional export does not discard a build.
 */

process CACTUS_PANGENOME {
    tag "taxid_${taxid}"
    label 'cactus_pangenome'

    // publish flat: cactus nests everything under out/; strip that prefix at publish time
    // (leaves the workdir untouched, so the graph itself is never moved/rewritten).
    publishDir "${params.outdir}/pangenome/${taxid}", mode: params.publish_dir_mode,
        saveAs: { fn -> fn.startsWith('out/') ? fn.substring(4) : fn }

    input:
    tuple val(taxid), val(ref_name), val(names), path(fastas), val(ref_contigs), path(capabilities)

    output:
    // Exact CLIP handles; the separate audit enforces required exports so a missing
    // export does not discard an otherwise successful expensive Cactus task.
    tuple val(taxid), path("out/${taxid}.gbz"),        emit: gbz, optional: true
    tuple val(taxid), path("out/${taxid}.og"),         emit: og, optional: true
    tuple val(taxid), path("out/${taxid}.hapl"),       emit: hapl, optional: true
    tuple val(taxid), path("out/${taxid}.snarls"),     emit: snarls, optional: true
    tuple val(taxid), path("out/${taxid}.gfa.gz"),     emit: gfa, optional: true
    tuple val(taxid), path("out/${taxid}.vcf.gz"),     emit: vcf, optional: true
    tuple val(taxid), path("out/${taxid}.vcf.gz.tbi"), emit: vcf_tbi, optional: true
    tuple val(taxid), path("out/${taxid}.raw.vcf.gz"), emit: raw_vcf, optional: true
    tuple val(taxid), path("out/${taxid}.chroms/*"),   emit: chrom_og, optional: true
    tuple val(taxid), path("out/${taxid}.gaf.gz"),      emit: gaf, optional: true
    tuple val(taxid), path("out/${taxid}.viz/*"),      emit: viz, optional: true
    // Optional construction products; FULL has no parallel biological analysis arm.
    tuple val(taxid), path("out/${taxid}.full.gbz"),        emit: gbz_full,      optional: true
    tuple val(taxid), path("out/${taxid}.full.og"),         emit: og_full,       optional: true
    tuple val(taxid), path("out/${taxid}.full.gfa.gz"),     emit: gfa_full,      optional: true
    tuple val(taxid), path("out/${taxid}.full.snarls"),     emit: snarls_full,   optional: true
    tuple val(taxid), path("out/${taxid}.full.vcf.gz"),     emit: vcf_full,      optional: true
    tuple val(taxid), path("out/${taxid}.full.vcf.gz.tbi"), emit: vcf_full_tbi,  optional: true
    tuple val(taxid), path("out/${taxid}.full.raw.vcf.gz"), emit: raw_vcf_full,  optional: true
    // per-chromosome graphs, both flavours. chrom_og above globs out/<taxid>.chroms/* and so
    // already includes the .full.og files; these give each flavour its own named channel so
    // consumers stop having to filter on the filename.
    tuple val(taxid), path("out/${taxid}.chroms/*.full.og"), emit: chrom_og_full, optional: true
    // ---- graph-reference handles (--gref) ---------------------------------------------
    // ADDITIONAL outputs: the clip graph is untouched, so nothing downstream changes. All
    // optional, matching this block's convention that a naming assumption must not fail a
    // multi-hour task.
    //
    // The VCF emit is a GLOB because --grefL embeds its value in the filename:
    // `.gref.vcf.gz` without it, `.gref95.vcf.gz` with 0.95. A fixed name would silently
    // emit nothing the moment that param is set.
    tuple val(taxid), path("out/${taxid}.gref.gbz"),        emit: gref_gbz,   optional: true
    tuple val(taxid), path("out/${taxid}.gref.gfa.gz"),     emit: gref_gfa,   optional: true
    tuple val(taxid), path("out/${taxid}.gref.hapl"),       emit: gref_hapl,  optional: true
    tuple val(taxid), path("out/${taxid}.gref*.vcf.gz"),    emit: gref_vcf,   optional: true
    tuple val(taxid), path("out/${taxid}.gref*.vcf.gz.tbi"), emit: gref_vcf_tbi, optional: true
    // catch-all: keeps + publishes everything cactus produced EXCEPT the construction scratch
    // removed in-script (full graphs, HAL, GAF/PAF, SV graph, raw VCF, stats bundle, ...).
    // '**' so future cactus outputs are retained automatically (general-purpose).
    tuple val(taxid), path("out/**"),                  emit: all
    tuple val(taxid), path("seqfile.txt"),             emit: seqfile
    path("versions.tsv"),                              emit: versions

    script:
    def extra = params.pangenome_cactus_extra ?: ''
    if (!ref_contigs || ref_contigs.any { !(it ==~ /[A-Za-z0-9_.+-]+/) })
        error 'Reference contigs must be explicit, nonempty finalized sequence IDs'
    if (names.size() != fastas.size() || names.toSet().size() != names.size() || !names.contains(ref_name))
        error 'Cactus input identity/cardinality mismatch'
    if (extra =~ /--(reference|refContigs|gref|grefL|vcf|vcfwave|vcfbub|haplo|gfa|gbz|odgi|chrom-og|viz|clip|collapse|outName|outDir)(?:\s|=|$)/)
        error 'pangenome_cactus_extra cannot override graph-role contract options'
    def chromog = params.pangenome_odgi_chromosomes ? '--chrom-og clip' : ''
    def vizopt = params.pangenome_odgi_visualization ? '--viz clip' : ''
    def vcfbub = params.pangenome_vcfbub_max_ref ?: 100000
    // batch 6 C1. A NAMED param rather than a line in pangenome_cactus_extra: this is a
    // controlled experiment whose whole value is that one variable changed, and a free-text
    // passthrough leaves no record of which run carried it. Named, it appears in the params
    // dump, can be asserted on, and is echoed into this task's log below -- so a graph can
    // always be traced back to the scoring it was built with.
    // --gref adds .gref.* outputs and leaves the clip graph alone. Source type defaults to
    // 'clip', the tool author's default and the flavour the rest of this pipeline analyses.
    def gref  = ( params.pangenome_gref ?: false )
                  ? "--gref ${params.pangenome_gref_source ?: 'clip'}" : ''
    // --grefL passes -L to vg deconstruct: cluster traversals whose handle Jaccard
    // coefficient is >= F. Allele MERGING, not nesting. Off by default so both VCFs
    // decompose identically and --gref is the only difference between the two arms.
    def grefl = ( (params.pangenome_gref ?: false) && params.pangenome_grefl )
                  ? "--grefL ${params.pangenome_grefl}" : ''
    def grefmin = ( (params.pangenome_gref ?: false) && params.pangenome_gref_min_len )
                  ? "--minGrefLen ${params.pangenome_gref_min_len}" : ''
    def lasttrain = params.pangenome_cactus_lasttrain ? '--lastTrain' : ''
    // -gpu image runs KegAlign automatically; --gpu 1 pins it to the single requested GPU
    // and --lastzMemory is the recommended cluster safeguard for the alignment jobs.
    def gpu   = params.pangenome_use_gpu ? '--gpu 1 --lastzMemory 100G' : ''
    """
    set -euo pipefail

    # cactus/Toil write config under ~/.toil and heavy scratch under TMPDIR. The container's
    # inherited HOME (/home/<user>) and /tmp are not writable/roomy on the compute node, so
    # point both at the Nextflow task dir (on /scratch): writable, discarded after the task.
    export HOME="\$PWD"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$TMPDIR"

    ids=(${names.join(' ')})
    fastas=(${fastas})

    # ---- two-column seqfile: <SAMPLE.HAP>  <fasta> ; locate the reference fasta ----
    : > seqfile.txt
    REF_FA=""
    for i in "\${!ids[@]}"; do
        printf '%s\\t%s\\n' "\${ids[\$i]}" "\${fastas[\$i]}" >> seqfile.txt
        [ "\${ids[\$i]}" = "${ref_name}" ] && REF_FA="\${fastas[\$i]}"
    done
    if [ -z "\${REF_FA}" ]; then
        echo "[PANGENOME ${taxid}] ERROR: reference '${ref_name}' not among inputs" >&2; exit 1
    fi

    # Reference chromosomes come from the final-assembly eligibility audit.
    # This also supports a chromosome-scale haploid assembly without harmonized chr names.
    REFCONTIGS="${ref_contigs.join(' ')}"
    echo "[PANGENOME ${taxid}] reference=${ref_name}; refContigs=\${REFCONTIGS}"

    # ---- run cactus (jobstore must not exist; workDir + jobstore in the task dir) ----
    # This label sets scratch = false, so the task dir is on /work: odgi_squeeze needs >110 GB
    # on top of a job store reaching 337 GB and the 447 GB node disk is not enough. Measured
    # twice, both failures at odgi_squeeze with ENOSPC, the second with the node to itself.
    # Both directories are removed explicitly at the end of a successful run, since nothing
    # discards the task dir now.
    rm -rf js cactus_work out
    mkdir -p cactus_work out

    # Scoring provenance in the graph build's own log. The default derives from HOXD70, which
    # the cactus docs describe as suited to VERY DIVERGED genomes and warn can produce "long
    # runs of transitions that really should be gaps" in a pangenome -- and this cohort is ten
    # haplotypes of one species. Which scoring built a given graph is not recoverable from the
    # graph afterwards, so it is recorded here.
    if [ -n "${lasttrain}" ]; then
        echo "[PANGENOME ${taxid}] alignment scoring: --lastTrain (trained on these inputs)" >&2
    else
        echo "[PANGENOME ${taxid}] alignment scoring: cactus default (HOXD70-derived)" >&2
    fi
    if [ -n "${gref}" ]; then
        echo "[PANGENOME ${taxid}] graph reference: ${gref} ${grefmin} ${grefl}" >&2
        echo "  .gref.* outputs are ADDITIONAL; the clip graph is unchanged, so no" >&2
        echo "  haplotype analysis needs gref_* exclusion. Confirm after the run:" >&2
        echo "    Biological sample identities are recorded in assembly_eligibility.tsv" >&2
    else
        echo "[PANGENOME ${taxid}] graph reference: disabled" >&2
    fi

    cactus-pangenome \\
        ./js \\
        seqfile.txt \\
        --workDir ./cactus_work \\
        --outDir ./out \\
        --outName ${taxid} \\
        --reference ${ref_name} \\
        --refContigs \${REFCONTIGS} \\
        --vcf clip \\
        --vcfbub ${vcfbub} \\
        --haplo \\
        --gfa full clip \\
        --gbz clip \\
        ${vizopt} \\
        --odgi clip \\
        ${chromog} \\
        --maxCores ${task.cpus} \\
        ${gpu} \\
        ${lasttrain} \\
        ${gref} \\
        ${grefmin} \\
        ${grefl} \\
        ${extra}

    # ---- cull construction scratch (pre-join per-chromosome intermediates, superseded by
    # the joined graph; the bulk of the file count) and the duplicate seqfile. Everything
    # else cactus produced is kept and published (flattened) via the output block above.
    rm -rf out/chrom-subproblems out/chrom-alignments
    rm -f  out/seqfile.txt
    rm -rf cactus_work js

    CV=\$(cactus --version 2>&1 | awk 'NR==1{print}')
    printf 'process\\ttool\\tversion\\n%s\\tcactus\\t%s\\n' "${task.process}" "\${CV}" > versions.tsv
    """

    stub:
    """
    mkdir -p out out/${taxid}.chroms out/${taxid}.viz
    for x in gbz og hapl snarls gfa.gz vcf.gz vcf.gz.tbi raw.vcf.gz; do : > "out/${taxid}.\$x"; done
    # full-graph twins, so -stub-run exercises the same channel topology as a real run
    for x in full.gbz full.og full.gfa.gz full.snarls full.vcf.gz full.vcf.gz.tbi full.raw.vcf.gz; do
        : > "out/${taxid}.\$x"
    done
    : > out/${taxid}.chroms/chr1_1.full.og
    : > out/${taxid}.chroms/chr1_1.og
    : > out/${taxid}.viz/chr1_1.viz.png
    : > seqfile.txt
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
