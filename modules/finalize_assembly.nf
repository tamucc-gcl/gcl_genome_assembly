/*
========================================================================================
    FINALIZE ASSEMBLY MODULE
========================================================================================
    Repo location: modules/finalize_assembly.nf

    Classifies scaffolds vs unplaced contigs by size, renames by descending size
    (scaffold_1..N, then contig_1..N), writes sorted FASTA + name map + index.
    Deterministic ordering prevents contact-map display artifacts from length reorder.

    Input : tuple(meta, assembly_fasta)   (ch_final_assembly: teloclip- or gap-fill output)
    Output: assembly / name_map / fai      (per-haplotype, meta.id)
========================================================================================
*/

process FINALIZE_ASSEMBLY {
    tag "${meta.id}"
    label 'finalize_assembly'

    publishDir "${params.outdir}/assembly/final", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(assembly_fasta, stageAs: 'input/*')

    output:
    tuple val(meta), path("${meta.id}.fasta"),          emit: assembly
    tuple val(meta), path("${meta.id}.name_map.tsv"),   emit: name_map
    tuple val(meta), path("${meta.id}.fasta.fai"),      emit: fai

    script:
    // Sequences >= this size are treated as chromosomal scaffolds.
    // Anything smaller is classified as an unplaced contig.
    // Default 1 Mb; override with params.finalize_min_scaffold_bp
    def min_scaffold = params.finalize_min_scaffold_bp ?: 1000000
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
    echo -e "old_name\\tnew_name\\tlength\\tclass" > ${meta.id}.name_map.tsv
    cat name_map_full.tsv >> ${meta.id}.name_map.tsv

    # ------------------------------------------------------------------
    # 4. Extract sequences in sorted order and rename headers
    #
    #    samtools faidx extracts in the order we request, then we pipe
    #    through awk to swap in the new names.
    # ------------------------------------------------------------------
    # Ordered list of old names (scaffolds first, then contigs — by descending size)
    cut -f1 name_map_full.tsv > extract_order.txt

    # Rename headers using awk exact-match lookup (no regex escaping needed)
    samtools faidx "\${INPUT_FA}" \$(cat extract_order.txt | tr '\\n' ' ') \\
        | awk '
            BEGIN { while ((getline < "name_map_full.tsv") > 0) map[\$1] = \$2 }
            /^>/ {
                old = substr(\$1, 2)  # strip leading ">"
                if (old in map) print ">" map[old]
                else print
                next
            }
            { print }
          ' > ${meta.id}.fasta

    # ------------------------------------------------------------------
    # 5. Index the final output
    # ------------------------------------------------------------------
    samtools faidx ${meta.id}.fasta

    # ------------------------------------------------------------------
    # 6. Summary log
    # ------------------------------------------------------------------
    N_SCAFF=\$(awk '\$4 == "scaffold"' name_map_full.tsv | wc -l)
    N_CTG=\$(awk '\$4 == "unplaced"' name_map_full.tsv | wc -l)
    echo "[FINALIZE] ${meta.id}: \${N_SCAFF} scaffolds + \${N_CTG} unplaced contigs"
    echo "[FINALIZE] Scaffold size threshold: ${min_scaffold} bp"
    """

    stub:
    """
    touch ${meta.id}.fasta
    touch ${meta.id}.fasta.fai
    echo -e "old_name\\tnew_name\\tlength\\tclass" > ${meta.id}.name_map.tsv
    """
}
