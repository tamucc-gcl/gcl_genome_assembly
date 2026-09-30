/* Missing required roles fail here, independently of the expensive cached build. */
process PANGENOME_GRAPH_CONTRACT {
    cpus 1
    memory '1 GB'
    time '10m'
    scratch false
    publishDir "${params.outdir}/pangenome/${taxid}", mode: params.publish_dir_mode
    input:
    tuple val(taxid), val(roles), path(artifacts)
    output:
    tuple val(taxid), path('pangenome_manifest.tsv'), emit: manifest
    tuple val(taxid), path('pangenome_build_report.md'), emit: report
    tuple val(taxid), path('contract.ok'), emit: ready
    script:
    def required = ['biological_gbz', 'biological_gfa', 'biological_odgi', 'haplotype_index',
                    'variant_gbz', 'variant_gfa', 'variants_standard', 'variants_standard_index',
                    'full_qc_gfa', 'clipping_stats']
    def missing = required.findAll { !roles[it] }
    def checks = required.findAll { roles[it] }.collect { "test -s '${roles[it]}'" }.join('\n')
    def graphOf = { r -> r.startsWith('variant') ? 'gref_clip' :
        (r.startsWith('full_') || r == 'clipping_stats' ? 'full_support' : 'clip') }
    def lines = roles.findAll { k, v -> v }.collect { k, v -> "${k}\t${graphOf(k)}\t${v}" }.join('\n')
    def rows = roles.findAll { k, v -> v }.collect { k, v -> "| ${k} | ${graphOf(k)} | [${v}](${v}) |" }.join('\n')
    """
    set -euo pipefail
    if [ '${missing.size()}' -ne 0 ]; then
        echo 'Missing required graph products: ${missing.join(', ')}. Inspect cached Cactus outputs; no fallback graph will be substituted.' >&2
        exit 1
    fi
    ${checks}
    cat > pangenome_manifest.tsv <<'MANIFEST'
role	graph	file
${lines}
MANIFEST
    cat > pangenome_build_report.md <<'GRAPH_REPORT'
<details>
<summary>Pangenome construction: taxid ${taxid}</summary>

Graph-role contract passed. CLIP is the biological graph; GREF(CLIP) supplies the
standard variant coordinate representation. FULL is retained for compact QC support.
The standard GREF VCF is Cactus's vcfbub output; no automatic vcfwave/fine catalog is run.
Synthetic GREF paths must not enter biological haplotype or individual denominators.
Reference-coordinate and synthetic-coordinate records must not be pooled blindly as
independent sites. Refer to [the identity ledger](pangenome_identity.tsv) for
biological sample grouping and assembly checksums.

| Role | Graph | File |
|---|---|---|
${rows}

This batch provides construction contracts and basic graph statistics.
Panacus sharing, chromosome-by-haplotype summaries, ordination, and the revised
structural-variant analysis are not yet part of this replacement analysis workflow.
Their absence is not a zero result.

</details>
GRAPH_REPORT
    touch contract.ok
    """
}
