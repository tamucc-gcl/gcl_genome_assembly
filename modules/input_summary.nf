process INPUT_SUMMARY {
    tag 'input_validation'
    cpus 1
    memory '1 GB'
    time '10m'
    publishDir "${params.outdir}/pipeline", mode: params.publish_dir_mode

    input:
    val records

    output:
    path 'input_validation.md', emit: report
    path 'input_validation.tsv', emit: table

    script:
    def clean = { value -> value.toString().replaceAll(/[\r\n\t]/, ' ').replace('|', '/') }
    def rows = records.collect { r -> "${clean(r.sample)}\t${r.status}\t${r.expected}\t${clean(r.reason)}" }.join('\n')
    def md = records.collect { r -> "| ${clean(r.sample)} | ${r.status} | ${r.expected} | ${clean(r.reason)} |" }.join('\n')
    """
    cat > input_validation.tsv <<'INPUT_TABLE'
sample_id\tstatus\texpected_assemblies\treason
${rows}
INPUT_TABLE
    cat > input_validation.md <<'INPUT_REPORT'
# Input validation

Accepted means input validation passed; it does not mean assembly or biological QC passed.

<details>
<summary>Sample decisions</summary>

| Sample | Status | Expected assemblies | Reason |
|---|---|---:|---|
${md}

</details>
INPUT_REPORT
    """
}
