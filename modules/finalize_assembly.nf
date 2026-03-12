/*
========================================================================================
    FINALIZE ASSEMBLY MODULE
========================================================================================
    Takes the final assembly from the pipeline and produces a clean, standardized
    output suitable for downstream analysis and submission.

    Operations performed:
      1. Classify sequences as chromosomal scaffolds vs unplaced contigs
         based on a minimum scaffold size threshold
      2. Rename chromosomal scaffolds by descending size: scaffold_1, scaffold_2, ...
      3. Rename unplaced contigs by descending size: contig_1, contig_2, ...
      4. Write all sequences in sorted order (scaffolds first, then contigs)
      5. Produce a name mapping file (old name → new name) for traceability
      6. Index the final FASTA

    The input channel is ch_final_assembly, which is set in main.nf to either
    TELOCLIP_EXTEND.out.extended_assembly or GAP_FILLING.out.filled_assembly.

    Input:
    - tuple(haplotype_id, assembly_fasta)

    Output:
    - assembly:     tuple(haplotype_id, final_fasta)        — renamed, sorted FASTA
    - name_map:     tuple(haplotype_id, name_map_tsv)       — old_name → new_name mapping
    - fai:          tuple(haplotype_id, fasta_index)         — samtools faidx index
========================================================================================
*/

process FINALIZE_ASSEMBLY {
    tag "${haplotype_id}"
    label 'finalize_assembly'

    publishDir "${params.outdir}/assembly/final", mode: params.publish_dir_mode

    input:
    tuple val(haplotype_id), path(assembly_fasta, stageAs: 'input/*')

    output:
    tuple val(haplotype_id), path("${haplotype_id}.fasta"),          emit: assembly
    tuple val(haplotype_id), path("${haplotype_id}.name_map.tsv"),   emit: name_map
    tuple val(haplotype_id), path("${haplotype_id}.fasta.fai"),      emit: fai

    script:
    // Sequences >= this size are treated as chromosomal scaffolds.
    // Anything smaller is classified as an unplaced contig.
    // Default 1 Mb; override with params.finalize_min_scaffold_size
    def min_scaffold = params.finalize_min_scaffold_size ?: 1000000
    """
    set -euo pipefail

    INPUT_FA="${assembly_fasta.name}"

    # ------------------------------------------------------------------
    # 1. Index the input so we have sequence names + lengths
    # ------------------------------------------------------------------
    samtools faidx "\${INPUT_FA}"

    # ------------------------------------------------------------------
    # 2. Build the name mapping: classify, sort by size, rename
    #
    #    Output columns: old_name  new_name  length  class
    # ------------------------------------------------------------------
    awk -v min_scaf=${min_scaffold} '
    BEGIN { OFS = "\\t" }
    {
        name = \$1; len = \$2
        if (len >= min_scaf) {
            scaff[++ns] = name
            scaff_len[ns] = len
        } else {
            ctg[++nc] = name
            ctg_len[nc] = len
        }
    }
    END {
        # Sort scaffolds by descending length (simple insertion sort — fine for <100 seqs)
        for (i = 2; i <= ns; i++) {
            j = i
            while (j > 1 && scaff_len[j] > scaff_len[j-1]) {
                t = scaff[j];     scaff[j] = scaff[j-1];     scaff[j-1] = t
                t = scaff_len[j]; scaff_len[j] = scaff_len[j-1]; scaff_len[j-1] = t
                j--
            }
        }
        # Sort contigs by descending length
        for (i = 2; i <= nc; i++) {
            j = i
            while (j > 1 && ctg_len[j] > ctg_len[j-1]) {
                t = ctg[j];     ctg[j] = ctg[j-1];     ctg[j-1] = t
                t = ctg_len[j]; ctg_len[j] = ctg_len[j-1]; ctg_len[j-1] = t
                j--
            }
        }
        # Emit scaffolds then contigs
        for (i = 1; i <= ns; i++)
            print scaff[i], "scaffold_" i, scaff_len[i], "scaffold"
        for (i = 1; i <= nc; i++)
            print ctg[i], "contig_" i, ctg_len[i], "unplaced"
    }' "\${INPUT_FA}.fai" > name_map_full.tsv

    # ------------------------------------------------------------------
    # 3. Create the name mapping output (old_name → new_name)
    # ------------------------------------------------------------------
    echo -e "old_name\\tnew_name\\tlength\\tclass" > ${haplotype_id}.name_map.tsv
    cat name_map_full.tsv >> ${haplotype_id}.name_map.tsv

    # ------------------------------------------------------------------
    # 4. Extract sequences in sorted order and rename headers
    #
    #    samtools faidx extracts in the order we request, then we pipe
    #    through awk to swap in the new names.
    # ------------------------------------------------------------------
    # Ordered list of old names (scaffolds first, then contigs — by descending size)
    cut -f1 name_map_full.tsv > extract_order.txt

    # Build sed rename commands from mapping
    awk 'BEGIN{OFS="\\t"} {
        # Escape any regex-special chars in old name
        old = \$1; new = \$2
        gsub(/[\\/.\\*\\[\\]\\^\\$]/, "\\\\\\\\&", old)
        print "s/^>" old "\\\\b/>" new "/"
    }' name_map_full.tsv > rename.sed

    samtools faidx "\${INPUT_FA}" \$(cat extract_order.txt | tr '\\n' ' ') \\
        | sed -f rename.sed \\
        > ${haplotype_id}.fasta

    # ------------------------------------------------------------------
    # 5. Index the final output
    # ------------------------------------------------------------------
    samtools faidx ${haplotype_id}.fasta

    # ------------------------------------------------------------------
    # 6. Summary log
    # ------------------------------------------------------------------
    N_SCAFF=\$(awk '\$4 == "scaffold"' name_map_full.tsv | wc -l)
    N_CTG=\$(awk '\$4 == "unplaced"' name_map_full.tsv | wc -l)
    echo "[FINALIZE] ${haplotype_id}: \${N_SCAFF} scaffolds + \${N_CTG} unplaced contigs"
    echo "[FINALIZE] Scaffold size threshold: ${min_scaffold} bp"
    """

    stub:
    """
    touch ${haplotype_id}.fasta
    touch ${haplotype_id}.fasta.fai
    echo -e "old_name\\tnew_name\\tlength\\tclass" > ${haplotype_id}.name_map.tsv
    """
}