/*
========================================================================================
    PANGENOME QC MODULE  (workstream C, Tier 1 — graph-intrinsic)
========================================================================================
    Repo location: modules/pangenome_qc.nf

    Quality diagnostics computed from the finished graph alone (no external inputs):
      - acyclicity            : vg stats on the .gbz  (a pangenome with inversions is
                                normally NOT acyclic -- informational, not a defect)
      - node degree / depth   : odgi degree / odgi depth summaries (tangle / collapse signal)
      - tangling / linearity  : odgi stats -l -s after odgi sort -O (mean links length,
                                sum of path-node distances)
      - re-alignment edit rate: parsed from the cactus .gaf (identity = Smatches/Saln-length;
                                low edit rate => graph faithfully represents the inputs)

    All metrics are best-effort: each command's raw output is captured to <taxid>.qc_raw.txt
    (published for inspection / to finalize any parsing), and the parsed values land in
    <taxid>.qc_metrics.tsv with NA fallbacks so the process never fails on a single metric.

    Input : tuple(taxid, gbz, og, gaf)   (clip .gbz + clip .og + cactus .gaf.gz)
    Output: qc_metrics.tsv / qc_raw.txt / versions
========================================================================================
*/

process PANGENOME_QC {
    tag "taxid_${taxid}"
    label 'cactus_tools'

    publishDir "${params.outdir}/pangenome/${taxid}", mode: params.publish_dir_mode

    input:
    tuple val(taxid), path(gbz), path(og), path(gaf)

    output:
    tuple val(taxid), path("${taxid}.qc_metrics.tsv"), emit: metrics
    tuple val(taxid), path("${taxid}.qc_raw.txt"),     emit: raw, optional: true
    path("versions.tsv"),                              emit: versions

    script:
    """
    set -euo pipefail
    export HOME="\$PWD"
    raw=${taxid}.qc_raw.txt
    : > "\$raw"

    # ---- acyclicity (vg) -- flag per plan; verify against `vg stats --help` -------------
    echo "### vg stats -A (is-acyclic)" >> "\$raw"
    acyc=\$(vg stats -A ${gbz} 2>>"\$raw" | tr -d '\\n' || true)
    [ -z "\$acyc" ] && acyc=NA
    echo "\$acyc" >> "\$raw"

    # ---- linearity / degree / depth (odgi; need an optimized graph) --------------------
    if odgi sort -i ${og} -o sorted.og -O 2>>"\$raw"; then S=sorted.og; else S=${og}; fi

    echo "### odgi stats -l -s (linearity)" >> "\$raw"
    odgi stats -i "\$S" -l -s >> "\$raw" 2>>"\$raw" || echo "odgi stats linearity: failed" >> "\$raw"

    echo "### odgi degree -S (node degree summary)" >> "\$raw"
    odgi degree -i "\$S" -S >> "\$raw" 2>>"\$raw" || echo "odgi degree: failed" >> "\$raw"

    echo "### odgi depth (node depth; head)" >> "\$raw"
    odgi depth -i "\$S" 2>>"\$raw" | head -20 >> "\$raw" || echo "odgi depth: failed" >> "\$raw"

    # ---- re-alignment edit rate (GAF: col10 matches, col11 aln block length) ------------
    echo "### edit rate from gaf" >> "\$raw"
    read editrate alnbp naln < <(zcat ${gaf} 2>>"\$raw" | awk -F'\\t' '
        \$11+0>0 { m+=\$10; a+=\$11; n++ }
        END { if(a>0) printf "%.6f %d %d\\n", 1-m/a, a, n; else print "NA 0 0" }')
    echo "edit_rate=\$editrate aligned_bp=\$alnbp n_alignments=\$naln" >> "\$raw"

    identity=NA
    [ "\$editrate" != "NA" ] && identity=\$(awk -v e="\$editrate" 'BEGIN{printf "%.6f", 1-e}')

    # extract a couple of parseable odgi metrics from the raw dump (lenient) --------------
    mll=\$(grep -iE 'mean.?links.?length' "\$raw" | grep -oE '[0-9]+(\\.[0-9]+)?' | tail -1 || true); [ -z "\$mll" ] && mll=NA
    spd=\$(grep -iE 'sum.?of.?path.?node.?(distances|dist)' "\$raw" | grep -oE '[0-9]+(\\.[0-9]+)?' | tail -1 || true); [ -z "\$spd" ] && spd=NA

    {
      printf 'metric\\tvalue\\n'
      printf 'is_acyclic\\t%s\\n'          "\$acyc"
      printf 'edit_rate\\t%s\\n'           "\$editrate"
      printf 'graph_identity\\t%s\\n'      "\$identity"
      printf 'realigned_bp\\t%s\\n'        "\$alnbp"
      printf 'n_alignments\\t%s\\n'        "\$naln"
      printf 'mean_links_length\\t%s\\n'   "\$mll"
      printf 'sum_path_node_dist\\t%s\\n'  "\$spd"
    } > ${taxid}.qc_metrics.tsv

    rm -f sorted.og
    {
      printf 'process\\ttool\\tversion\\n'
      printf '%s\\tvg\\t%s\\n'   "${task.process}" "\$(vg version 2>&1 | head -n1)"
      printf '%s\\todgi\\t%s\\n' "${task.process}" "\$(odgi version 2>&1 | head -n1)"
    } > versions.tsv
    """

    stub:
    """
    printf 'metric\\tvalue\\n' > ${taxid}.qc_metrics.tsv
    : > ${taxid}.qc_raw.txt
    printf 'process\\ttool\\tversion\\n' > versions.tsv
    """
}
