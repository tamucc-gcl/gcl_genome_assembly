/*
========================================================================================
    HI-C SCAFFOLDING MODULE (YaHS) - UNIFIED FOR ALL ROUNDS
========================================================================================
    Input : filtered Hi-C BAM (valid pairs) + BAI + assembly FASTA + round identifier + round_params
    Output: scaffolded FASTA + final AGP + YaHS BIN + log
    
    The round parameter controls output directory and file naming:
    - round = "round1" or "" → scaffolding/yahs/
    - round = "round2" → scaffolding/yahs_round2/
    - round = "roundN" → scaffolding/yahs_roundN/
    
    The round_params map contains all YaHS parameters for this specific round
========================================================================================
*/

process SCAFFOLD_HIC {
    tag "${haplotype_id}_${round ?: 'round1'}"
    label 'scaffold_hic'

    publishDir "${params.outdir}/scaffolding/yahs${round && round != 'round1' ? '_' + round : ''}/", 
        mode: params.publish_dir_mode

    input:
    tuple val(haplotype_id), path(bam), path(bai), path(assembly_fasta), val(round), val(round_params)

    output:
    tuple val(haplotype_id), path("${haplotype_id}${round && round != 'round1' ? '_' + round : ''}_scaffolds.fa"), emit: scaffolds
    tuple val(haplotype_id), path("${haplotype_id}${round && round != 'round1' ? '_' + round : ''}_scaffolds_final.agp"), emit: agp
    tuple val(haplotype_id), path("${haplotype_id}${round && round != 'round1' ? '_' + round : ''}_scaffolds_final.bin"), emit: bin
    tuple val(haplotype_id), path("${haplotype_id}${round && round != 'round1' ? '_' + round : ''}.yahs.log"), emit: log

    script:
    // Extract parameters from round_params map
    def min_contig_len    = round_params.min_contig_length ?: 10000
    def min_mapq          = round_params.min_mapq ?: 1
    def resolutions       = round_params.resolutions ?: '10000,20000,50000,100000,200000,500000,1000000,2000000,5000000,10000000,20000000,50000000,100000000,200000000,500000000'
    def rounds_per_res    = round_params.rounds_per_resolution ?: null
    def enzyme            = round_params.enzyme ?: null
    def no_contig_ec      = round_params.no_contig_ec ? true : false
    def no_scaffold_ec    = round_params.no_scaffold_ec ? true : false

    // Use a prefix that includes round identifier if present
    def round_suffix = round && round != 'round1' ? "_${round}" : ""
    def prefix = "${haplotype_id}_yahs${round_suffix}"

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
    echo "[YAHS ${round ?: 'round1'}] Starting scaffolding for ${haplotype_id}"
    echo "[YAHS ${round ?: 'round1'}] Input assembly: ${assembly_fasta}"
    echo "[YAHS ${round ?: 'round1'}] Hi-C BAM: ${bam}"
    echo "[YAHS ${round ?: 'round1'}] Parameters:"
    echo "  min_contig_length: ${min_contig_len}"
    echo "  min_mapq: ${min_mapq}"
    echo "  resolutions: ${resolutions}"
    echo "  rounds_per_resolution: ${rounds_per_res ?: 'default'}"
    echo "  enzyme: ${enzyme ?: 'none'}"
    echo "  no_contig_ec: ${no_contig_ec}"
    echo "  no_scaffold_ec: ${no_scaffold_ec}"
    
    yahs "${assembly_fasta}" "${bam}" "\${yahs_args[@]}" 2>&1 | tee "${haplotype_id}${round_suffix}.yahs.log"

    # Standardize names
    mv "${prefix}_scaffolds_final.fa"  "${haplotype_id}${round_suffix}_scaffolds.fa"
    mv "${prefix}_scaffolds_final.agp" "${haplotype_id}${round_suffix}_scaffolds_final.agp"

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

    cp -f "\${bin_candidate}" "${haplotype_id}${round_suffix}_scaffolds_final.bin"
    
    echo "[YAHS ${round ?: 'round1'}] Scaffolding complete for ${haplotype_id}"
    """
    
    stub:
    def round_suffix = round && round != 'round1' ? "_${round}" : ""
    """
    touch ${haplotype_id}${round_suffix}_scaffolds.fa
    touch ${haplotype_id}${round_suffix}_scaffolds_final.agp
    touch ${haplotype_id}${round_suffix}_scaffolds_final.bin
    touch ${haplotype_id}${round_suffix}.yahs.log
    """
}