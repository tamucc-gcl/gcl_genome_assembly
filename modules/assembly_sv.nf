process ASSEMBLY_SV_PLAN {
    cpus 1
    memory '2 GB'
    time '1h'
    conda 'conda-forge::python=3.11'
    publishDir "${params.outdir}/post_assembly/sv", mode: params.publish_dir_mode
    input:
    path manifest
    path planner
    output:
    path 'sv_plan.json', emit: plan
    path 'sv_plan.tsv', emit: table
    path 'sv_plan.md', emit: report
    script:
    def selected = (params.assembly_sv_queries ?: '').toString().split(',').collect { it.trim() }.findAll { it }
    def selection_json = groovy.json.JsonOutput.toJson(selected)
    """
    cat > selection.json <<'SV_SELECTION'
${selection_json}
SV_SELECTION
    python3 ${planner} ${manifest} --selection selection.json
    """
}

process ASSEMBLY_SV_ALIGN {
    tag "${meta.taxid}:${meta.reference}:${meta.query}"
    label 'assembly_sv'
    maxForks 1
    input:
    tuple val(meta), path(reference, stageAs: 'input_ref/*'), path(query, stageAs: 'input_query/*')
    output:
    tuple val(meta), path('reference.fa'), path('query.fa'), path('out.filtered.delta'), path('out.filtered.coords'), path('alignment_provenance.json'), path('alignment_inputs.sha256'), emit: alignment
    tuple val(meta), path('alignment_provenance.json'), emit: provenance
    path 'alignment_versions.tsv', emit: versions
    publishDir "${params.outdir}/post_assembly/sv/${meta.taxid}/${meta.key}", mode: params.publish_dir_mode,
        pattern: '*provenance.json'
    script:
    def payload = groovy.json.JsonOutput.toJson(meta)
    """
    set -euo pipefail
    cat > alignment_provenance.json <<'SV_META'
${payload}
SV_META
    python3 - <<'PY'
import json, subprocess
names = json.load(open('alignment_provenance.json'))['chromosomes']
for source, target in [('${reference}', 'reference.fa'), ('${query}', 'query.fa')]:
    # Index staged copies only; retain native final-assembly sequence names.
    subprocess.run(['samtools', 'faidx', source], check=True)
    with open(target, 'w') as out:
        subprocess.run(['samtools', 'faidx', source, *names], stdout=out, check=True)
PY
    # MUMmer 4 reference-unique anchors are the default (no --mumreference option).
    nucmer --threads ${task.cpus} -c 100 -b 500 -l 50 reference.fa query.fa
    delta-filter -m -i 90 -l 100 out.delta > out.filtered.delta
    show-coords -THrd out.filtered.delta > out.filtered.coords
    test -s out.filtered.coords
    sha256sum ${reference} ${query} reference.fa query.fa > alignment_inputs.sha256
    printf 'nucmer\\t%s\\n' "\$(nucmer --version 2>&1)" > alignment_versions.tsv
    """
}

process ASSEMBLY_SV_CALL {
    tag "${meta.taxid}:${meta.reference}:${meta.query}"
    label 'assembly_sv'
    maxForks 1
    publishDir "${params.outdir}/post_assembly/sv/${meta.taxid}/${meta.key}", mode: params.publish_dir_mode,
        pattern: 'sv/*', saveAs: { it.replaceFirst('^sv/', '') }
    input:
    tuple val(meta), path(reference), path(query), path(delta), path(coords), path(provenance), path(hashes)
    path auditor
    output:
    tuple val(meta), path('sv/sv_qc.md'), emit: report
    tuple val(meta), path('sv/*'), emit: files
    path 'sv/versions.tsv', emit: versions
    script:
    """
    set -euo pipefail
    syri -c ${coords} -d ${delta} -r ${reference} -q ${query}
    test -s syri.out
    test -s syri.summary
    test -s syri.vcf
    mkdir sv
    python3 ${auditor} . sv --near-bp ${params.sv_qc_near_bp} \\
      --min-assigned-fraction ${params.sv_qc_min_assigned_fraction}
    cp syri.out syri.summary syri.vcf sv/
    cp ${provenance} ${hashes} sv/
    printf 'syri\\t%s\\n' "\$(syri --version 2>&1)" > sv/versions.tsv
    """
}
