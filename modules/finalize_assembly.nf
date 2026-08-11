/*
========================================================================================
    FINALIZE ASSEMBLY MODULE
========================================================================================
    Repo location: modules/finalize_assembly.nf

    Renames scaffolds and writes a sorted FASTA + name map + index (per-haplotype).

    Two naming modes, selected by the third input:
      * name_map basename == 'NO_HARMONIZE'  -> classify vs size and rename by descending
        size (scaffold_1..N, then contig_1..N). Original behaviour, unchanged.
      * otherwise (a harmonized map from HARMONIZE_SCAFFOLDS) -> rename + reorient to the
        homology-derived names (chrN_p / chrA_i+chrB_j / unplaced_n). Rename-and-reorient
        ONLY: scaffold count, boundaries and gaps are invariant. Single-chromosome
        scaffolds flagged 'rev' are reverse-complemented; composites are left native.

    Input : tuple(meta, assembly_fasta, name_map)
    Output: assembly / name_map / fai   (the full harmonized map is published by
    HARMONIZE_SPECIES; FINALIZE emits the 4-column old->new name_map)
========================================================================================
*/

process FINALIZE_ASSEMBLY {
    tag "${meta.id}"
    label 'finalize_assembly'

    publishDir "${params.outdir}/assembly/final", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(assembly_fasta, stageAs: 'input/*'), path(name_map)

    output:
    tuple val(meta), path("${meta.id}.fasta"),          emit: assembly
    tuple val(meta), path("${meta.id}.name_map.tsv"),   emit: name_map
    tuple val(meta), path("${meta.id}.fasta.fai"),      emit: fai

    script:
    // Sequences >= this size are treated as chromosomal scaffolds in the size-rank
    // (NO_HARMONIZE) path. Default 1 Mb; override with params.finalize_min_scaffold_bp
    def min_scaffold = params.finalize_min_scaffold_bp ?: 1000000
    """
    set -euo pipefail

    INPUT_FA="${assembly_fasta.name}"
    MAP="${name_map.name}"

    # ------------------------------------------------------------------
    # 1. Index the input so we have sequence names + lengths
    # ------------------------------------------------------------------
    samtools faidx "\${INPUT_FA}"

    if [ "\${MAP}" = "NO_HARMONIZE" ]; then
        # ==============================================================
        # SIZE-RANK NAMING (original behaviour)
        # ==============================================================
        # Classify (scaffold >= min_scaf, else unplaced), sort by descending length within
        # class (scaffolds first), assign sequential names.
        awk -v min_scaf=${min_scaffold} 'BEGIN { OFS = "\\t" }
            { cls = (\$2 >= min_scaf ? "0scaffold" : "1unplaced"); print cls, \$2, \$1 }' "\${INPUT_FA}.fai" \\
          | sort -k1,1 -k2,2nr \\
          | awk 'BEGIN { OFS = "\\t" }
            {
                if (\$1 == "0scaffold") { new = "scaffold_" ++ns; lab = "scaffold" }
                else                    { new = "contig_"   ++nc; lab = "unplaced" }
                print \$3, new, \$2, lab
            }' > name_map_full.tsv

        echo -e "old_name\\tnew_name\\tlength\\tclass" > ${meta.id}.name_map.tsv
        cat name_map_full.tsv >> ${meta.id}.name_map.tsv

        cut -f1 name_map_full.tsv > extract_order.txt

        samtools faidx -r extract_order.txt "\${INPUT_FA}" \\
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

        N_SCAFF=\$(awk '\$4 == "scaffold"' name_map_full.tsv | wc -l)
        N_CTG=\$(awk '\$4 == "unplaced"' name_map_full.tsv | wc -l)
        echo "[FINALIZE] ${meta.id}: \${N_SCAFF} scaffolds + \${N_CTG} unplaced contigs (size-rank)"
        echo "[FINALIZE] Scaffold size threshold: ${min_scaffold} bp"

    else
        # ==============================================================
        # HARMONIZED NAMING (rename + reorient; NO re-scaffolding)
        #   map cols: old_name new_name orient order length class ref_span flags
        # ==============================================================
        NROWS=\$(tail -n +2 "\${MAP}" | wc -l)
        NSEQ=\$(wc -l < "\${INPUT_FA}.fai")
        if [ "\${NROWS}" -ne "\${NSEQ}" ]; then
            echo "[FINALIZE ${meta.id}] ERROR: harmonized map rows (\${NROWS}) != input sequences (\${NSEQ})" >&2
            exit 1
        fi

        # rows in output order (col 4)
        tail -n +2 "\${MAP}" | sort -t\$'\\t' -k4,4n > rows.tsv

        : > ${meta.id}.fasta
        while IFS=\$'\\t' read -r old new orient order length cls span flags; do
            if [ "\${orient}" = "rev" ]; then
                samtools faidx "\${INPUT_FA}" "\${old}" \\
                  | seqkit seq -r -p -t DNA -w 60 2>/dev/null \\
                  | awk -v n=">\${new}" 'NR==1{print n; next} {print}' >> ${meta.id}.fasta
            else
                samtools faidx "\${INPUT_FA}" "\${old}" \\
                  | awk -v n=">\${new}" 'NR==1{print n; next} {print}' >> ${meta.id}.fasta
            fi
        done < rows.tsv

        # 4-column compat name map (chromosome/composite -> scaffold; unplaced -> unplaced)
        echo -e "old_name\\tnew_name\\tlength\\tclass" > ${meta.id}.name_map.tsv
        tail -n +2 "\${MAP}" \\
          | awk 'BEGIN{FS=OFS="\\t"} { c=(\$6=="unplaced"?"unplaced":"scaffold"); print \$1,\$2,\$5,c }' \\
          >> ${meta.id}.name_map.tsv

        N_CHR=\$(tail -n +2 "\${MAP}"  | awk -F'\\t' '\$6=="chromosome"' | wc -l)
        N_COMP=\$(tail -n +2 "\${MAP}" | awk -F'\\t' '\$6=="composite"'  | wc -l)
        N_UNPL=\$(tail -n +2 "\${MAP}" | awk -F'\\t' '\$6=="unplaced"'   | wc -l)
        echo "[FINALIZE] ${meta.id}: \${N_CHR} chromosome + \${N_COMP} composite + \${N_UNPL} unplaced (harmonized)"
    fi

    # ------------------------------------------------------------------
    # Index the final output
    # ------------------------------------------------------------------
    samtools faidx ${meta.id}.fasta
    """

    stub:
    """
    touch ${meta.id}.fasta
    touch ${meta.id}.fasta.fai
    echo -e "old_name\\tnew_name\\tlength\\tclass" > ${meta.id}.name_map.tsv
    """
}
