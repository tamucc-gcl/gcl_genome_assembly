process ASSEMBLY_RUN_SUMMARY {
    tag 'assembly_status'
    cpus 1
    memory '1 GB'
    time '10m'
    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    val expected
    val completed

    output:
    path 'assembly_run_summary.md', emit: report

    script:
    def clean = { value -> value.toString().replaceAll(/[\r\n\t]/, ' ').replace('|', '/') }
    def bySample = completed.groupBy { it.sample }
    def rows = expected.collect { record ->
        def observed = bySample[record.sample] ?: []
        def state = record.status == 'accepted'
            ? (observed.size() == record.expected ? 'finalized' : 'incomplete') : record.status
        "| ${clean(record.sample)} | ${state} | ${observed.size()}/${record.expected} | ${clean(record.reason)} |"
    }.join('\n')
    def artifacts = completed.sort { it.id }.collect { item ->
        "- [${item.id}](assembly/final/${item.name})"
    }.join('\n')
    def state = expected.any { r -> r.status == 'invalid' || (r.status == 'accepted' && (bySample[r.sample] ?: []).size() != r.expected) }
        ? 'Incomplete: inspect sample statuses and execution logs.' : 'All accepted samples have their expected finalized assemblies.'
    """
    cat > assembly_run_summary.md <<'ASSEMBLY_REPORT'
# Assembly run summary

${state}

QC mode: ${params.qc_mode}. Pangenome requested: ${params.run_pangenome}.
Finalized describes workflow completion, not biological quality or pangenome eligibility.

<details>
<summary>Sample status</summary>

| Sample | Status | Finalized/expected | Input decision |
|---|---|---:|---|
${rows}

</details>

<details>
<summary>Final assembly files</summary>

${artifacts}

</details>
ASSEMBLY_REPORT
    """
}
