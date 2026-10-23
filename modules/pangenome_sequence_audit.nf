process PANGENOME_SEQUENCE_AUDIT {
    tag "${taxid}"
    label 'cactus_tools'
    cpus 1
    memory '32 GB'
    time '4h'
    publishDir "${params.outdir}/pangenome/${taxid}/representation", mode: params.publish_dir_mode,
        pattern: 'sequence_audit/*', saveAs: { it.replaceFirst('^sequence_audit/', '') }
    input:
    tuple val(taxid), path(full), path(clip), path(policy)
    path(adapter)
    path(helper)
    output:
    tuple val(taxid), path('sequence_audit/sequence_graph_audit.md'), emit: report
    tuple val(taxid), path('sequence_audit/*'), emit: files
    tuple val(taxid), path('sequence_audit/sequence_graph_audit.json'), emit: ready
    script:
    """
    set -euo pipefail
    python3 ${adapter} --full ${full} --clip ${clip} --policy ${policy} --out sequence_audit
    """
    stub:
    """
    mkdir sequence_audit
    printf '# Stub only: graph sequence checks were not executed.\\n' > sequence_audit/sequence_graph_audit.md
    printf '{"status":"STUB"}\\n' > sequence_audit/sequence_graph_audit.json
    """
}
