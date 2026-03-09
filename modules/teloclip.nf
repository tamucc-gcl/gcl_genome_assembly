/*
========================================================================================
    TELOCLIP — TELOMERE EXTENSION MODULE
========================================================================================
    Uses teloclip (Taranto) to recover missing telomeric sequences from
    soft-clipped HiFi read alignments and extend scaffold ends.

    Workflow:
      1. Map raw HiFi reads to gap-filled scaffolds (minimap2 map-hifi)
      2. Filter for soft-clipped alignments at scaffold ends with telomeric
         motifs (teloclip filter)
      3. Extend scaffolds using overhang sequences (teloclip extend)

    Single process handles the full pipeline to avoid staging large BAMs
    between processes. The intermediate BAM is created and consumed in the
    same work directory.

    Input:
      tuple(haplotype_id, scaffold_fasta, hifi_fastq)

    Output:
      tuple(haplotype_id, extended_fasta)  — scaffolds with telomeres appended
      tuple(haplotype_id, stats_report)    — per-contig extension statistics
========================================================================================
*/

process TELOCLIP_EXTEND {
    tag "${haplotype_id}"
    label 'teloclip'

    publishDir "${params.outdir}/assembly/scaffold/teloclip", mode: params.publish_dir_mode

    input:
    tuple val(haplotype_id), path(scaffold_fasta), path(hifi_fastq)

    output:
    tuple val(haplotype_id), path("${haplotype_id}.teloclip_extended.fasta"), emit: extended_assembly
    tuple val(haplotype_id), path("${haplotype_id}.teloclip_stats.tsv"),      emit: stats
    tuple val(haplotype_id), path("${haplotype_id}.teloclip_overhangs.bam"),  emit: overhangs_bam
    tuple val(haplotype_id), path("${haplotype_id}.teloclip_overhangs.bam.bai"), emit: overhangs_bai

    script:
    def motif       = params.telomere_motif ?: 'TTAGGG'
    def min_clip    = params.teloclip_min_clip ?: 1
    def max_break   = params.teloclip_max_break ?: 50
    def min_anchor  = params.teloclip_min_anchor ?: 100
    def min_mapq    = params.teloclip_min_mapq ?: 20
    def min_ovh     = params.teloclip_min_overhangs ?: 1
    def max_homopol = params.teloclip_max_homopolymer ?: 500
    """
    set -euo pipefail

    echo "[TELOCLIP] Starting telomere extension for ${haplotype_id}"
    echo "[TELOCLIP] Motif: ${motif}  min_clip: ${min_clip}  max_break: ${max_break}"
    echo "[TELOCLIP] min_anchor: ${min_anchor}  min_mapq: ${min_mapq}"
    echo ""

    # ---- Step 1: Index the scaffold assembly ----
    samtools faidx ${scaffold_fasta}

    # ---- Step 2: Map raw HiFi reads, filter for overhangs with telomere motifs ----
    #   -ax map-hifi : HiFi preset for minimap2
    #   -F 0x100     : exclude secondary alignments
    #   --motifs     : only keep overhangs containing the telomeric motif
    minimap2 \\
        -ax map-hifi \\
        -t ${task.cpus} \\
        ${scaffold_fasta} \\
        ${hifi_fastq} \\
    | samtools view -h -F 0x100 -q ${min_mapq} \\
    | teloclip filter \\
        --ref-idx ${scaffold_fasta}.fai \\
        --motifs ${motif} \\
        --min-clip ${min_clip} \\
        --max-break ${max_break} \\
        --min-anchor ${min_anchor} \\
    | samtools sort -@ ${task.cpus} -o ${haplotype_id}.teloclip_overhangs.bam

    samtools index ${haplotype_id}.teloclip_overhangs.bam

    # ---- Step 3: Count overhang reads (QC) ----
    OVERHANG_COUNT=\$(samtools view -c ${haplotype_id}.teloclip_overhangs.bam)
    echo "[TELOCLIP] Found \${OVERHANG_COUNT} telomere-containing overhang alignments"

    # ---- Step 4: Extend scaffolds ----
    if [ "\${OVERHANG_COUNT}" -gt 0 ]; then
        teloclip extend \\
            ${haplotype_id}.teloclip_overhangs.bam \\
            ${scaffold_fasta} \\
            --output-fasta ${haplotype_id}.teloclip_extended.fasta \\
            --stats-report ${haplotype_id}.teloclip_stats.tsv \\
            --min-overhangs ${min_ovh} \\
            --max-homopolymer ${max_homopol} \\
            --max-break ${max_break} \\
            --min-anchor ${min_anchor} \\
            --count-motifs ${motif},\$(echo ${motif} | rev | tr 'ACGTacgt' 'TGCAtgca') \\
            --exclude-outliers

        
# ---- Step 5: Parse human-readable stats report into TSV ----
    # teloclip extend --stats-report produces a formatted text report, not TSV.
    # Convert it to the TSV format expected by generate_summary_report.R
    mv ${haplotype_id}.teloclip_stats.tsv ${haplotype_id}.teloclip_stats_raw.txt

    python3 <<'PYEOF'
import re, csv, sys

hap_id = "${haplotype_id}"
infile  = f"{hap_id}.teloclip_stats_raw.txt"
outfile = f"{hap_id}.teloclip_stats.tsv"

rows = []
current_contig = None
original_len = None

with open(infile) as fh:
    for line in fh:
        line = line.rstrip()

        # Match contig header: ### scaffold_2
        m = re.match(r'^### (\S+)', line)
        if m:
            current_contig = m.group(1)
            original_len = None
            continue

        if current_contig is None:
            continue

        # Original length: 91,837,727
        m = re.match(r'^\s+Original length:\s+([\d,]+)', line)
        if m:
            original_len = int(m.group(1).replace(',', ''))
            continue

        # Left extension: +699bp from read ...
        # Right extension: +17884bp from read ...
        m = re.match(r'^\s+(Left|Right) extension:\s+\+(\d+)bp', line)
        if m and original_len is not None:
            end = m.group(1).lower()
            ext_len = int(m.group(2))
            # Prefix contig with haplotype_id so the R script can extract it
            contig_name = f"{hap_id}_{current_contig}"
            rows.append({
                "contig": contig_name,
                "contig_length": original_len,
                "end": end,
                "extension_length": ext_len,
                "overhang_count": "NA",
                "motif_counts": "NA",
            })

with open(outfile, "w", newline="") as out:
    fields = ["contig", "contig_length", "end", "extension_length",
              "overhang_count", "motif_counts"]
    w = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
    w.writeheader()
    w.writerows(rows)

print(f"[TELOCLIP] Parsed {len(rows)} extensions into {outfile}")
PYEOF
    else
        echo "[TELOCLIP] No overhang reads found — copying input assembly unchanged"
        cp ${scaffold_fasta} ${haplotype_id}.teloclip_extended.fasta
        # Create empty stats report with header
        printf "contig\\tcontig_length\\tend\\textension_length\\toverhang_count\\tmotif_counts\\n" \\
            > ${haplotype_id}.teloclip_stats.tsv
    fi

    echo "[TELOCLIP] Done: ${haplotype_id}"
    """

    stub:
    """
    touch ${haplotype_id}.teloclip_extended.fasta
    printf "contig\\tcontig_length\\tend\\textension_length\\toverhang_count\\tmotif_counts\\n" > ${haplotype_id}.teloclip_stats.tsv
    touch ${haplotype_id}.teloclip_overhangs.bam
    touch ${haplotype_id}.teloclip_overhangs.bam.bai
    """
}

/*
========================================================================================
    COLLECT_TELOCLIP_STATS
========================================================================================
    Aggregates per-haplotype teloclip extension statistics into a single
    summary file for the report.
========================================================================================
*/
process COLLECT_TELOCLIP_STATS {
    tag "collect_teloclip"
    label 'summarize_assembly'

    publishDir "${params.outdir}/assembly/scaffold/teloclip", mode: params.publish_dir_mode

    input:
    path("stats/*")

    output:
    path("all_teloclip_stats.tsv"), emit: stats

    script:
    """
    set -euo pipefail

    first_file=\$(ls stats/*.tsv 2>/dev/null | head -1)
    if [[ -n "\$first_file" ]]; then
        head -1 "\$first_file" > all_teloclip_stats.tsv
        for f in stats/*.tsv; do
            tail -n +2 "\$f" >> all_teloclip_stats.tsv
        done
    else
        echo "No teloclip stats files found" > all_teloclip_stats.tsv
    fi
    """

    stub:
    """
    touch all_teloclip_stats.tsv
    """
}