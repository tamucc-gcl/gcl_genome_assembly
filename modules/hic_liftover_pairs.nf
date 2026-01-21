
/*
========================================================================================
    HI-C PAIRS LIFTOVER MODULE (AGP-BASED; NO REMAPPING)
========================================================================================
    Converts Hi-C pairs.gz from contig coordinates to scaffold coordinates using a scaffold AGP.

    Intended use:
    - Map Hi-C reads once to the contig assembly (BAM)
    - Generate pairs.gz (pairtools split) on contig coordinates
    - Liftover those pairs to scaffold coordinates using the YaHS (or other) AGP
    - Run pairtools/cooler QC on the scaffold-coordinate pairs without remapping

    Notes:
    - Uses agptools transform to lift BED coordinates with AGP
    - Re-sorts the output pairs (pairtools sort) after coordinate translation
    - Flips strand (+/-) for reads mapped to contigs placed in reverse orientation in the AGP
========================================================================================
*/

nextflow.enable.dsl = 2

process HIC_LIFTOVER_PAIRS {
    tag "${haplotype_id}_${qc_label}"
    label 'hic_liftover_pairs'

    publishDir "${params.outdir}/qc/hic_mapping/${qc_label}/liftover", mode: params.publish_dir_mode

    input:
    tuple val(haplotype_id), val(qc_label), path(pairs_gz), path(scaffold_agp), path(scaffold_fasta)

    output:
    tuple val(haplotype_id), val(qc_label), path("${haplotype_id}_${qc_label}.scaffold.pairs.gz"), emit: pairs

    script:
    """
    set -euo pipefail
    export LC_ALL=C

    TMPDIR="\${TMPDIR:-\$PWD}"

    # -------------------------------------------------------------------------
    # 1) Ensure scaffold chrom sizes (for sorting + downstream cooler)
    # -------------------------------------------------------------------------
    if [ ! -s "${scaffold_fasta}.fai" ]; then
      samtools faidx ${scaffold_fasta}
    fi
    cut -f1,2 ${scaffold_fasta}.fai > scaffold.chrom.sizes

    # -------------------------------------------------------------------------
    # 2) Identify contigs/components placed in reverse orientation in the AGP
    # -------------------------------------------------------------------------
    # AGP format (W lines): object object_beg object_end part_number component_id comp_beg comp_end orientation
    awk '\$5=="W" && \$9=="-" {print \$6}' ${scaffold_agp} | sort -u > rev_components.txt

    # -------------------------------------------------------------------------
    # 3) Split header vs body
    # -------------------------------------------------------------------------
    zcat ${pairs_gz} | awk 'BEGIN{OFS="\\t"} /^#/ {print > "header.txt"; next} {print > "body.tsv"}'

    # -------------------------------------------------------------------------
    # 4) Make 1bp BED records for each end (keyed by line number)
    # pairs columns assumed: readID chr1 pos1 chr2 pos2 strand1 strand2 ...
    # IMPORTANT: Use 1-based coordinates to match AGP (no -1 subtraction)
    # -------------------------------------------------------------------------
    awk 'BEGIN{OFS="\\t"} {print \$2, \$3, \$3, NR}' body.tsv > end1.bed
    awk 'BEGIN{OFS="\\t"} {print \$4, \$5, \$5, NR}' body.tsv > end2.bed

    # -------------------------------------------------------------------------
    # 5) Lift each end with AGP
    # -------------------------------------------------------------------------
    agptools transform end1.bed ${scaffold_agp} > end1.scaf.bed
    agptools transform end2.bed ${scaffold_agp} > end2.scaf.bed

    # Sort by NR key to ensure consistent join
    sort -k4,4n end1.scaf.bed > end1.scaf.sorted
    sort -k4,4n end2.scaf.bed > end2.scaf.sorted

    # -------------------------------------------------------------------------
    # 6) Rebuild pairs body with scaffold coordinates + strand flip on reversed contigs
    # -------------------------------------------------------------------------
    paste end1.scaf.sorted end2.scaf.sorted body.tsv \\
      | awk -v OFS="\\t" '
          BEGIN {
            while ((getline < "rev_components.txt") > 0) rev[\$1]=1
          }
          {
            # end1.scaf: \$1 \$2 \$3 \$4 ; end2.scaf: \$5 \$6 \$7 \$8 ; body: starts at \$9
            # agptools outputs: chr start end name (1-based coordinates)
            # We want the end position (\$3 and \$7) for 1-based pairs format
            new_chr1=\$1; new_pos1=\$3
            new_chr2=\$5; new_pos2=\$7

            # Body columns (starting at \$9):
            #  \$9  readID
            #  \$10 chr1
            #  \$11 pos1
            #  \$12 chr2
            #  \$13 pos2
            #  \$14 strand1
            #  \$15 strand2
            old_chr1=\$10
            old_chr2=\$12
            s1=\$14
            s2=\$15

            if (rev[old_chr1]) s1 = (s1=="+" ? "-" : "+")
            if (rev[old_chr2]) s2 = (s2=="+" ? "-" : "+")

            # Rewrite in-place
            \$10=new_chr1; \$11=new_pos1
            \$12=new_chr2; \$13=new_pos2
            \$14=s1;      \$15=s2

            # Print body fields only (preserve any additional columns after strand2)
            for (i=9; i<=NF; i++) {
              printf "%s%s", \$i, (i==NF ? ORS : OFS)
            }
          }
        ' > lifted.body.tsv

    # -------------------------------------------------------------------------
    # 7) Re-sort pairs in scaffold coordinate space and bgzip
    #    NOTE: pairtools sort does NOT accept --chroms-path
    #    It determines chromosome order from the input pairs header
    # -------------------------------------------------------------------------
    cat header.txt lifted.body.tsv \\
      | pairtools sort --nproc ${task.cpus} --tmpdir "\${TMPDIR}" \\
      | bgzip -c > ${haplotype_id}_${qc_label}.scaffold.pairs.gz
    """

    stub:
    """
    touch ${haplotype_id}_${qc_label}.scaffold.pairs.gz
    """
}