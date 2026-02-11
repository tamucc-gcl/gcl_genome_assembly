/*
========================================================================================
    PAIRWISE GENOME ALIGNMENT MODULE
========================================================================================
    Performs pairwise whole-genome alignment between two assemblies using minimap2
    
    Input:
    - haplotype_id1: Identifier for reference assembly
    - assembly1: Reference assembly FASTA
    - haplotype_id2: Identifier for query assembly
    - assembly2: Query assembly FASTA
    
    Output:
    - paf: Filtered PAF alignment file (gzipped)
    - qc: QC summary statistics TSV
    - log: Minimap2 log file
    
    Use cases:
    - Comparing haplotypes within a sample
    - Cross-sample genome comparisons
    - Input for dotplot visualization
========================================================================================
*/

process PAIRWISE_ALIGNMENT {
    tag "${haplotype_id1}_vs_${haplotype_id2}"
    label 'pairwise_alignment'
    
    publishDir "${params.outdir}/pairwise_alignments", mode: params.publish_dir_mode
    
    input:
    tuple val(haplotype_id1), path(assembly1), val(haplotype_id2), path(assembly2)
    
    output:
    tuple val(haplotype_id1), val(haplotype_id2), path("${haplotype_id1}_vs_${haplotype_id2}.paf.gz"), emit: paf
    tuple val(haplotype_id1), val(haplotype_id2), path("${haplotype_id1}_vs_${haplotype_id2}.qc.tsv"), emit: qc
    tuple val(haplotype_id1), val(haplotype_id2), path("${haplotype_id1}_vs_${haplotype_id2}.log"), emit: log
    
    script:
    def preset = params.pairwise_alignment_preset ?: 'asm5'
    def min_mapq = params.pairwise_alignment_min_mapq ?: 5
    def min_aln_bp = params.pairwise_alignment_min_aln_bp ?: 10000
    def out_prefix = "${haplotype_id1}_vs_${haplotype_id2}"
    """
    set -euo pipefail
    
    # QC header
    echo -e "ref_id\\tqry_id\\tpreset\\tmin_mapq\\tmin_aln_bp\\tn_align\\tsum_aln_bp\\tsum_match_bp\\tmean_identity\\tmean_dv\\tmean_mapq" > "${out_prefix}.qc.tsv"
    
    # Run minimap2 alignment, filter, compute QC stats, and compress
    minimap2 \\
        -t ${task.cpus} \\
        -x ${preset} \\
        --eqx \\
        --secondary=no \\
        ${assembly1} \\
        ${assembly2} \\
        2> "${out_prefix}.log" \\
    | awk -v mq="${min_mapq}" -v ml="${min_aln_bp}" 'BEGIN{FS=OFS="\\t"} \$12 >= mq && \$11 >= ml' \\
    | tee >( \\
        awk -v ref="${haplotype_id1}" \\
            -v qry="${haplotype_id2}" \\
            -v preset="${preset}" \\
            -v mq="${min_mapq}" \\
            -v ml="${min_aln_bp}" '
        BEGIN {
            FS=OFS="\\t"
            n=0
            sum_aln=0
            sum_match=0
            sum_mapq=0
            sum_dv_w=0
            sum_id_w=0
        }
        {
            n++
            aln=\$11
            sum_aln   += aln
            sum_match += \$10
            sum_mapq  += \$12

            dv=""
            for(i=13;i<=NF;i++){
                if(\$i ~ /^dv:f:/){
                    dv=\$i
                    sub(/^dv:f:/,"",dv)
                    break
                }
            }
            if(dv=="") dv=0

            sum_dv_w += dv * aln
            sum_id_w += (1 - dv) * aln
        }
        END {
            mean_mapq = (n>0 ? sum_mapq/n : 0)
            mean_dv   = (sum_aln>0 ? sum_dv_w/sum_aln : 0)
            mean_id   = (sum_aln>0 ? sum_id_w/sum_aln : 0)

            print ref, qry, preset, mq, ml, n, sum_aln, sum_match, mean_id, mean_dv, mean_mapq
        }' >> "${out_prefix}.qc.tsv" \\
    ) \\
    | pigz -p ${task.cpus} -c > "${out_prefix}.paf.gz"
    
    # Ensure output exists even if no alignments pass filters
    if [ ! -s "${out_prefix}.paf.gz" ]; then
        echo "" | pigz -c > "${out_prefix}.paf.gz"
    fi
    """
    
    stub:
    def out_prefix = "${haplotype_id1}_vs_${haplotype_id2}"
    """
    touch "${out_prefix}.paf.gz"
    echo -e "ref_id\\tqry_id\\tpreset\\tmin_mapq\\tmin_aln_bp\\tn_align\\tsum_aln_bp\\tsum_match_bp\\tmean_identity\\tmean_dv\\tmean_mapq" > "${out_prefix}.qc.tsv"
    touch "${out_prefix}.log"
    """
}