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
========================================================================================
*/

process ESTIMATE_GENOME_SIZE {
    tag "${meta.sample}"
    label 'genomescope'

    publishDir "${params.outdir}/qc/genome_size", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.sample}_genomescope"),     emit: results
    tuple val(meta), path("${meta.sample}.genome_size.txt"), emit: size

    script:
    def k        = params.kmer_size    ?: 21
    def ploidy   = meta.ploidy ?: (meta.n_hap ?: 2)   // organism ploidy, per-sample from the sheet
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

    # Parse a single genome-size value (bp) for downstream use.
    # NOTE: verify against your GenomeScope summary.txt layout — this grabs the MAX of the
    # "Genome Haploid Length" min/max range; switch to \$4 for the min or compute a mean.
    if [ -f ${meta.sample}_genomescope/summary.txt ]; then
        grep 'Genome Haploid Length' ${meta.sample}_genomescope/summary.txt \\
            | sed 's/,//g' \\
            | awk '{ print \$(NF-1) }' > ${meta.sample}.genome_size.txt
    else
        echo "NA" > ${meta.sample}.genome_size.txt
    fi
    """

    stub:
    """
    mkdir -p ${meta.sample}_genomescope
    touch ${meta.sample}_genomescope/summary.txt
    touch ${meta.sample}_genomescope/linear_plot.png
    touch ${meta.sample}_genomescope/transformed_linear_plot.png
    touch ${meta.sample}_genomescope/model.txt
    echo "415000000" > ${meta.sample}.genome_size.txt
    """
}
