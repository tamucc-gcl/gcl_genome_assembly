/*
========================================================================================
    HI-C SCAFFOLDING MODULE (YaHS)
========================================================================================
    Input : filtered Hi-C BAM (valid pairs) + BAI + assembly FASTA
    Output: scaffolded FASTA + final AGP + YaHS BIN + log
========================================================================================
*/

process SCAFFOLD_HIC {
    tag "${haplotype_id}"
    label 'scaffold_hic'

    publishDir "${params.outdir}/scaffolding/yahs/", mode: params.publish_dir_mode

    input:
    tuple val(haplotype_id), path(bam), path(bai), path(assembly_fasta)

    output:
    tuple val(haplotype_id), path("${haplotype_id}_scaffolds.fa"),          emit: scaffolds
    tuple val(haplotype_id), path("${haplotype_id}_scaffolds_final.agp"),   emit: agp
    tuple val(haplotype_id), path("${haplotype_id}_scaffolds_final.bin"),   emit: bin
    tuple val(haplotype_id), path("${haplotype_id}.yahs.log"),              emit: log

    script:
    def min_contig_len    = params.yahs_min_contig_length ?: 10000
    def min_mapq          = params.yahs_min_mapq ?: 1

    // YaHS expects comma-separated resolutions; default mirrors YaHS docs
    def resolutions       = params.yahs_resolutions ?: '10000,20000,50000,100000,200000,500000,1000000,2000000,5000000,10000000,20000000,50000000,100000000,200000000,500000000'

    // Optional knobs
    def rounds_per_res    = params.yahs_rounds_per_resolution ?: null   // corresponds to -R if set
    def enzyme            = params.yahs_enzyme ?: null                  // corresponds to -e if set
    def no_contig_ec      = params.yahs_no_contig_ec ? true : false
    def no_scaffold_ec    = params.yahs_no_scaffold_ec ? true : false

    // Use a prefix that won't double "scaffolds"
    def prefix = "${haplotype_id}_yahs"

    """
    set -euo pipefail

    # Keep temp files on local scratch if Nextflow provides it
    export TMPDIR="\${TMPDIR:-\$PWD}"

    # Index assembly if needed
    if [[ ! -s "${assembly_fasta}.fai" ]]; then
        samtools faidx "${assembly_fasta}"
    fi

    # Build YaHS args
    yahs_args=()
    yahs_args+=("-o" "${prefix}")
    yahs_args+=("-l" "${min_contig_len}")
    yahs_args+=("-r" "${resolutions}")

    # NOTE: -q is suppressed for BAMs not name-sorted (YaHS behavior). Keep it,
    # but be aware it may not take effect depending on your BAM sort order.
    yahs_args+=("-q" "${min_mapq}")

    if [[ -n "${rounds_per_res}" && "${rounds_per_res}" != "null" ]]; then
        yahs_args+=("-R" "${rounds_per_res}")
    fi

    if [[ -n "${enzyme}" && "${enzyme}" != "null" ]]; then
        yahs_args+=("-e" "${enzyme}")
    fi

    if ${no_contig_ec}; then
        yahs_args+=("--no-contig-ec")
    fi

    if ${no_scaffold_ec}; then
        yahs_args+=("--no-scaffold-ec")
    fi

    # Run YaHS
    yahs "${assembly_fasta}" "${bam}" "\${yahs_args[@]}" 2>&1 | tee "${haplotype_id}.yahs.log"

    # Standardize names (no double "scaffolds" now)
    mv "${prefix}_scaffolds_final.fa"  "${haplotype_id}_scaffolds.fa"
    mv "${prefix}_scaffolds_final.agp" "${haplotype_id}_scaffolds_final.agp"

    # Capture the YaHS-generated BIN (name can vary by version/prefix)
    # Prefer a BIN matching our prefix if it exists.
    bin_candidate=""
    if ls "${prefix}"*.bin >/dev/null 2>&1; then
        bin_candidate=\$(ls -1 "${prefix}"*.bin | head -n 1)
    elif ls *.bin >/dev/null 2>&1; then
        bin_candidate=\$(ls -1 *.bin | head -n 1)
    fi

    if [[ -z "\${bin_candidate}" ]]; then
        echo "[ERROR] No .bin produced by YaHS was found in work directory." >&2
        exit 1
    fi

    cp -f "\${bin_candidate}" "${haplotype_id}_scaffolds_final.bin"
    """
    
    stub:
    """
    touch ${haplotype_id}_scaffolds.fa
    touch ${haplotype_id}_scaffolds_final.agp
    touch ${haplotype_id}_scaffolds_final.bin
    touch ${haplotype_id}.yahs.log
    """
}
