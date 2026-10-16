process SHARED_REFERENCE {
    tag "taxid_${taxid}"
    cpus 1
    memory '2 GB'
    time '2h'
    conda 'conda-forge::python=3.11'
    publishDir "${params.outdir}/references/${taxid}", mode: params.publish_dir_mode
    input:
    tuple val(taxid), val(selection), path(fasta), path(fai)
    path verifier
    path fasta_inspector
    output:
    tuple val(taxid), path('reference.json'), emit: manifest
    tuple val(taxid), path('reference.md'), emit: report
    script:
    def payload = groovy.json.JsonOutput.toJson(selection)
    """
    set -euo pipefail
    cat > selection.json <<'REFERENCE_SELECTION'
${payload}
REFERENCE_SELECTION
    python3 ${verifier} --selection selection.json --fasta '${fasta}' --fai '${fai}'
    """
}
