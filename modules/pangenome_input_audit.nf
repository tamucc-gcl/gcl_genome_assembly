process PANGENOME_INPUT_AUDIT {
    tag "taxid_${taxid}"
    cpus 1
    memory '2 GB'
    time '2h'
    conda 'conda-forge::python=3.11'
    publishDir "${params.outdir}/pangenome/${taxid}", mode: params.publish_dir_mode
    input:
    tuple val(taxid), val(ref), val(names), path(fastas, stageAs: 'input??/*'), val(contigs), val(cohort)
    path(auditor)
    output:
    tuple val(taxid), path('input.ok'), emit: ready
    tuple val(taxid), path('pangenome_identity.tsv'), emit: identities
    tuple val(taxid), path('pangenome_input_audit.json'), emit: audit
    tuple val(taxid), path('pangenome_input_report.md'), emit: report
    script:
    def paths = fastas instanceof List ? fastas : [fastas]
    def members = cohort.members.withIndex().collect { member, i -> member + [taxid: taxid, staged_fasta: paths[i].toString()] }
    def payload = groovy.json.JsonOutput.toJson(cohort + [members: members])
    """
    set -euo pipefail
    cat > cohort.json <<'COHORT_JSON'
${payload}
COHORT_JSON
    python3 ${auditor} --input cohort.json
    """
}
