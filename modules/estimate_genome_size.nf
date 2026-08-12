/*
========================================================================================
    ESTIMATE_GENOME_SIZE MODULE
========================================================================================
    K-mer based genome size / heterozygosity estimation (jellyfish -> GenomeScope2).
    Repo location: modules/estimate_genome_size.nf

    Runs concurrently with assembly in BOTH branches, keyed on the assembly reads:
      - long-read branch: HiFi FASTQ (from BAM_TO_FASTQ)
      - short-read branch: PE R1 + R2
    The input `reads` is therefore one or more files (coerced to a list below).

    Ploidy for the GenomeScope model is the ORGANISM ploidy (meta.ploidy), not the output
    haplotype count (meta.n_hap) — e.g. a diploid organism assembled collapsed from short
    reads is p=2 here but n_hap=1 for the assembly fork. Falls back to n_hap then 2 until
    meta.ploidy lands in the parser.

    Params (Elvis defaults; mirror in nextflow.config). Ploidy is per-sample (meta.ploidy),
    not a global param. Based on Jason's 1-estimate_genome_size.sh. Requires a 'genomescope'
    process label with a conda env bundling genomescope2 + jellyfish.

    Publishing (est_genome_size/):
      - FLAT: every final output lands directly in est_genome_size/, prefixed with the
        sample ID — <sample>.summary.txt, <sample>.model.txt, and the four
        <sample>.{,transformed_}{linear,log}_plot.png profiles. No per-sample subdir.
      - Per-round fit artifacts (round*/) and the verbose fit log (progress.txt) are dropped.
      - <sample>.genome_size.txt is emitted on `size` (feeds assembly + report) but is
        NOT published; the user-facing size lives in <sample>.summary.txt.
      - `results` carries <sample>.summary.txt (parsed downstream for het/repeat).
========================================================================================
*/

process ESTIMATE_GENOME_SIZE {
    tag "${meta.sample}"
    label 'genomescope'

    // Flat publish into est_genome_size/. Everything is a top-level, ID-prefixed file, so
    // saveAs sees a plain filename: drop only <sample>.genome_size.txt (emitted, not published).
    publishDir "${params.outdir}/est_genome_size", mode: params.publish_dir_mode,
        saveAs: { fn -> fn.endsWith('.genome_size.txt') ? null : fn }

    input:
    tuple val(meta), path(reads), val(ploidy)

    output:
    tuple val(meta), path("${meta.sample}.summary.txt"),     emit: results
    tuple val(meta), path("${meta.sample}.genome_size.txt"), emit: size
    tuple val(meta), path("${meta.sample}.model.txt"),       emit: model
    tuple val(meta), path("${meta.sample}.*_plot.png"),      emit: plots
    path "versions.tsv", emit: versions

    script:
    def k        = params.kmer_size    ?: 21
    def jf_hash  = params.jellyfish_hash_size ?: '5G'
    def read_list = (reads instanceof List ? reads : [reads])
    def read_streams = read_list.collect { "<(zcat -f ${it})" }.join(' ')
    """
    # Count canonical k-mers across all provided reads (HiFi fastq, or PE R1+R2)
    jellyfish count -C -m ${k} -s ${jf_hash} -t ${task.cpus} \\
        -o reads_k${k}.jf \\
        ${read_streams}

    jellyfish histo -t ${task.cpus} reads_k${k}.jf > reads_k${k}.histo

    genomescope2 \\
        -i reads_k${k}.histo \\
        -o ${meta.sample}_genomescope \\
        -k ${k} \\
        -p ${ploidy} \\
        --verbose

    # Keep only the final model: strip per-round fit artifacts.
    rm -rf ${meta.sample}_genomescope/round*

    # Flatten: pull the final files up to the task root, ID-prefixed, then drop the
    # GenomeScope dir (which also discards the verbose progress.txt fit log).
    mv ${meta.sample}_genomescope/summary.txt                 ${meta.sample}.summary.txt
    mv ${meta.sample}_genomescope/model.txt                   ${meta.sample}.model.txt
    mv ${meta.sample}_genomescope/linear_plot.png             ${meta.sample}.linear_plot.png
    mv ${meta.sample}_genomescope/log_plot.png                ${meta.sample}.log_plot.png
    mv ${meta.sample}_genomescope/transformed_linear_plot.png ${meta.sample}.transformed_linear_plot.png
    mv ${meta.sample}_genomescope/transformed_log_plot.png    ${meta.sample}.transformed_log_plot.png
    rm -rf ${meta.sample}_genomescope

    # Parse a single haploid genome-size value (bp) for downstream use (channel only,
    # feeds assembly + report; NOT published — user-facing size is in summary.txt).
    # NOTE: verify against your GenomeScope summary.txt layout — this grabs the MAX of the
    # "Genome Haploid Length" min/max range; switch to \$4 for the min or compute a mean.
    if [ -f ${meta.sample}.summary.txt ]; then
        grep 'Genome Haploid Length' ${meta.sample}.summary.txt \\
            | sed 's/,//g' \\
            | awk '{ print \$(NF-1) }' > ${meta.sample}.genome_size.txt
    else
        echo "NA" > ${meta.sample}.genome_size.txt
    fi
    printf 'Jellyfish\t%s\n' "\$(jellyfish --version 2>&1 | sed 's/jellyfish //')" > versions.tsv
    printf 'GenomeScope\t%s\n' "\$( genomescope2 -v 2>&1 | sed -n 's/^GenomeScope //p' | head -n1 || true )" >> versions.tsv
    """

    stub:
    """
    touch ${meta.sample}.summary.txt
    touch ${meta.sample}.model.txt
    touch ${meta.sample}.linear_plot.png
    touch ${meta.sample}.log_plot.png
    touch ${meta.sample}.transformed_linear_plot.png
    touch ${meta.sample}.transformed_log_plot.png
    echo "415000000" > ${meta.sample}.genome_size.txt
    touch versions.tsv
    """
}
