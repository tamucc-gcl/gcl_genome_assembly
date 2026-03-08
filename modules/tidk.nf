/*
========================================================================================
    TIDK — TELOMERE IDENTIFICATION TOOLKIT
========================================================================================
    Replaces the custom scan_telomeres.py with the published tidk toolkit
    (Brown et al. 2025, Bioinformatics).

    Three processes:
      TIDK_EXPLORE  — de novo telomere motif discovery (once per sample)
      TIDK_SEARCH   — windowed telomere repeat quantification (per haplotype)
      TIDK_PLOT     — SVG visualisation of telomere density (per haplotype)

    Plus a collector:
      COLLECT_TIDK_RESULTS — aggregate summaries + SVGs across haplotypes

    All processes use the same conda label so they share an environment.
========================================================================================
*/

/*
========================================================================================
    TIDK_EXPLORE
========================================================================================
    De-novo discovery of the telomeric repeat unit using a kmer-based approach.
    Runs on the ASSEMBLY (not raw reads) — one call per haplotype, but the
    results for both haplotypes of a sample should agree. The caller in main.nf
    can deduplicate to one-per-sample if desired.

    Input : tuple(haplotype_id, assembly_fasta)
    Output: tuple(haplotype_id, explore_tsv)  — ranked candidate repeats
========================================================================================
*/
process TIDK_EXPLORE {
    tag "${haplotype_id}"
    label 'tidk'

    publishDir "${params.outdir}/telomeres/explore", mode: params.publish_dir_mode

    input:
    tuple val(haplotype_id), path(assembly_fasta)

    output:
    tuple val(haplotype_id), path("${haplotype_id}_tidk_explore.tsv"), emit: explore

    script:
    def min_len = params.tidk_explore_minimum ?: 5
    def max_len = params.tidk_explore_maximum ?: 12
    """
    set -euo pipefail

    tidk explore \\
        --minimum ${min_len} \\
        --maximum ${max_len} \\
        ${assembly_fasta} \\
        > ${haplotype_id}_tidk_explore.tsv
    """

    stub:
    """
    printf "id\\trepeat\\tcount\\n" > ${haplotype_id}_tidk_explore.tsv
    """
}

/*
========================================================================================
    TIDK_SEARCH
========================================================================================
    Windowed telomere repeat quantification across each scaffold.
    Uses params.telomere_motif (default TTAGGG) — override per-species as needed.

    Input : tuple(haplotype_id, assembly_fasta)
    Output: tuple(haplotype_id, search_tsv)
            tuple(haplotype_id, search_bedgraph)  — for downstream tools
========================================================================================
*/
process TIDK_SEARCH {
    tag "${haplotype_id}"
    label 'tidk'

    publishDir "${params.outdir}/telomeres/search", mode: params.publish_dir_mode

    input:
    tuple val(haplotype_id), path(assembly_fasta)

    output:
    tuple val(haplotype_id), path("${haplotype_id}_telomeric_repeat_windows.tsv"), emit: search_tsv
    tuple val(haplotype_id), path("${haplotype_id}_telomeric_repeat_windows.bedgraph"), emit: search_bedgraph

    script:
    def motif  = params.telomere_motif ?: 'TTAGGG'
    def window = params.tidk_search_window ?: 10000
    """
    set -euo pipefail

    # TSV output
    tidk search \\
        --string ${motif} \\
        --window ${window} \\
        --output ${haplotype_id} \\
        --dir . \\
        --extension tsv \\
        ${assembly_fasta}

    # bedgraph output (same search, different format)
    tidk search \\
        --string ${motif} \\
        --window ${window} \\
        --output ${haplotype_id} \\
        --dir . \\
        --extension bedgraph \\
        ${assembly_fasta}
    """

    stub:
    """
    touch ${haplotype_id}_telomeric_repeat_windows.tsv
    touch ${haplotype_id}_telomeric_repeat_windows.bedgraph
    """
}

/*
========================================================================================
    TIDK_PLOT
========================================================================================
    Generate an SVG showing telomere density along each scaffold.

    Input : tuple(haplotype_id, search_tsv)  — from TIDK_SEARCH
    Output: tuple(haplotype_id, svg)
========================================================================================
*/
process TIDK_PLOT {
    tag "${haplotype_id}"
    label 'tidk'

    publishDir "${params.outdir}/telomeres/plots", mode: params.publish_dir_mode

    input:
    tuple val(haplotype_id), path(search_tsv)

    output:
    tuple val(haplotype_id), path("${haplotype_id}_tidk.svg"), emit: plot

    script:
    def height = params.tidk_plot_height ?: 200
    def width  = params.tidk_plot_width  ?: 1000
    """
    set -euo pipefail

    tidk plot \\
        --tsv ${search_tsv} \\
        --height ${height} \\
        --width ${width} \\
        --output ${haplotype_id}_tidk
    """

    stub:
    """
    touch ${haplotype_id}_tidk.svg
    """
}

/*
========================================================================================
    TIDK_SUMMARIZE
========================================================================================
    Parse the TIDK_SEARCH windowed TSV to produce a per-scaffold telomere
    presence/absence summary comparable to the old SCAN_TELOMERES output.

    This keeps the summary report section backward-compatible while adding
    richer tidk data.

    Logic:
      - For each scaffold, check first and last windows for motif counts
        above a threshold (params.telomere_min_repeats).
      - Report 5'/3'/both presence and overall percentage.

    Input : tuple(haplotype_id, search_tsv)
    Output: tuple(haplotype_id, summary_tsv)  — one row per haplotype
            tuple(haplotype_id, detail_tsv)   — one row per scaffold
========================================================================================
*/
process TIDK_SUMMARIZE {
    tag "${haplotype_id}"
    label 'summarize_assembly'

    input:
    tuple val(haplotype_id), path(search_tsv)

    output:
    tuple val(haplotype_id), path("${haplotype_id}.telomere_summary.tsv"),  emit: summary
    tuple val(haplotype_id), path("${haplotype_id}.telomeres.tsv"),         emit: telomeres

    script:
    def min_count = params.telomere_min_repeats ?: 10
    """
    set -euo pipefail

    python3 << 'PYEOF'
import csv, sys, collections

min_count = ${min_count}
hap_id = "${haplotype_id}"

# ---- Parse windowed TSV from tidk search ----
# Columns: id, window, forward_repeat_number, reverse_repeat_number, telomeric_repeat_number
scaffolds = collections.OrderedDict()  # scaffold -> {windows: [(start, fwd, rev, total)], length: int}

with open("${search_tsv}") as fh:
    reader = csv.DictReader(fh, delimiter="\\t")
    for row in reader:
        scaf = row["id"]
        win_start = int(row["window"])
        total = int(row["telomeric_repeat_number"]) if "telomeric_repeat_number" in row else (
            int(row.get("forward_repeat_number", 0)) + int(row.get("reverse_repeat_number", 0))
        )
        if scaf not in scaffolds:
            scaffolds[scaf] = {"windows": [], "max_end": 0}
        scaffolds[scaf]["windows"].append((win_start, total))
        end = win_start + ${params.tidk_search_window ?: 10000}
        if end > scaffolds[scaf]["max_end"]:
            scaffolds[scaf]["max_end"] = end

# ---- Determine telomere presence per scaffold ----
detail_rows = []
for scaf, data in scaffolds.items():
    wins = sorted(data["windows"], key=lambda x: x[0])
    length = data["max_end"]
    # First window = 5' end, last window = 3' end
    has_5 = 1 if (wins and wins[0][1] >= min_count) else 0
    has_3 = 1 if (wins and wins[-1][1] >= min_count) else 0
    has_both = 1 if (has_5 and has_3) else 0
    detail_rows.append({
        "haplotype_id": hap_id, "scaffold": scaf, "length": length,
        "telomere_5prime": has_5, "telomere_3prime": has_3, "telomere_both": has_both
    })

# ---- Write per-scaffold detail ----
with open(f"{hap_id}.telomeres.tsv", "w") as out:
    fields = ["haplotype_id", "scaffold", "length", "telomere_5prime", "telomere_3prime", "telomere_both"]
    w = csv.DictWriter(out, fieldnames=fields, delimiter="\\t")
    w.writeheader()
    w.writerows(detail_rows)

# ---- Write summary ----
total = len(detail_rows)
total_len = sum(r["length"] for r in detail_rows)
five = sum(r["telomere_5prime"] for r in detail_rows)
three = sum(r["telomere_3prime"] for r in detail_rows)
both = sum(r["telomere_both"] for r in detail_rows)
five_only = five - both
three_only = three - both
none = total - five - three + both
pct = f"{(five + three - both) / total * 100:.2f}" if total > 0 else "0.00"

with open(f"{hap_id}.telomere_summary.tsv", "w") as out:
    fields = ["haplotype_id", "scaffolds", "total_length", "telomere_5prime", "telomere_3prime",
              "telomere_both", "telomere_5prime_only", "telomere_3prime_only", "telomere_none", "pct_with_telomere"]
    w = csv.DictWriter(out, fieldnames=fields, delimiter="\\t")
    w.writeheader()
    w.writerow({
        "haplotype_id": hap_id, "scaffolds": total, "total_length": total_len,
        "telomere_5prime": five, "telomere_3prime": three, "telomere_both": both,
        "telomere_5prime_only": five_only, "telomere_3prime_only": three_only,
        "telomere_none": none, "pct_with_telomere": pct
    })
PYEOF
    """

    stub:
    """
    printf "haplotype_id\\tscaffolds\\ttotal_length\\ttelomere_5prime\\ttelomere_3prime\\ttelomere_both\\ttelomere_5prime_only\\ttelomere_3prime_only\\ttelomere_none\\tpct_with_telomere\\n" > ${haplotype_id}.telomere_summary.tsv
    printf "haplotype_id\\tscaffold\\tlength\\ttelomere_5prime\\ttelomere_3prime\\ttelomere_both\\n" > ${haplotype_id}.telomeres.tsv
    """
}

/*
========================================================================================
    COLLECT_TIDK_RESULTS
========================================================================================
    Combines per-haplotype tidk outputs into single summary files.
    Drop-in replacement for the old COLLECT_TELOMERE_RESULTS.
========================================================================================
*/
process COLLECT_TIDK_RESULTS {
    tag "collect_tidk"
    label 'summarize_assembly'

    publishDir "${params.outdir}/telomeres", mode: params.publish_dir_mode

    input:
    path("summaries/*")
    path("telomeres/*")

    output:
    path("all_telomere_summaries.tsv"), emit: summary
    path("all_telomeres.tsv"),          emit: telomeres

    script:
    """
    set -euo pipefail

    # --- Combine summary files ---
    first_summary=\$(ls summaries/*.tsv 2>/dev/null | head -1)
    if [[ -n "\$first_summary" ]]; then
        head -1 "\$first_summary" > all_telomere_summaries.tsv
        for f in summaries/*.tsv; do
            tail -n +2 "\$f" >> all_telomere_summaries.tsv
        done
    else
        echo "No telomere summary files found" > all_telomere_summaries.tsv
    fi

    # --- Combine per-scaffold detail files ---
    first_telo=\$(ls telomeres/*.tsv 2>/dev/null | head -1)
    if [[ -n "\$first_telo" ]]; then
        head -1 "\$first_telo" > all_telomeres.tsv
        for f in telomeres/*.tsv; do
            tail -n +2 "\$f" >> all_telomeres.tsv
        done
    else
        echo "No telomere detail files found" > all_telomeres.tsv
    fi
    """

    stub:
    """
    touch all_telomere_summaries.tsv
    touch all_telomeres.tsv
    """
}