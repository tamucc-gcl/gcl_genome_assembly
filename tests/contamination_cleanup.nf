nextflow.enable.dsl = 2

params.cleanup_fixture = null
params.cleanup_scenario = 'baseline'

include { CONTAMINATION_POST_FCS } from '../workflows/contamination_post_fcs.nf'

process CHECK_CLEANUP_INPUT {
    tag "${meta.id}"
    conda 'conda-forge::python=3.11'
    cpus 1
    memory '1 GB'
    time '10m'
    publishDir "${params.outdir}/downstream", mode: 'copy',
        pattern: 'consumer.json', saveAs: { "${meta.id}.json" }

    input:
    tuple val(meta), path(fasta), val(expectation)
    path checker

    output:
    tuple val(meta.id), path('consumer.json')

    script:
    """
    cat > expected.json <<'EXPECTED'
${groovy.json.JsonOutput.toJson(expectation)}
EXPECTED
    python3 ${checker} verify-one --fasta '${fasta}' --expected expected.json --out consumer.json
    """
}

workflow {
    if (params.run_contamination_evidence)
        error 'This acceptance fixture must run with contamination evidence disabled'
    def fixture = new groovy.json.JsonSlurper().parseText(file(params.cleanup_fixture).text)
    def inputs = Channel.fromList(fixture.cases).map { row ->
        // Exercise the same interpolated metadata IDs as real haplotype assemblies.
        tuple([id: "${row.id}", sample: 'fixture', taxid: 373251], file(row.source))
    }
    CONTAMINATION_POST_FCS(inputs, Channel.empty(), Channel.empty())
    def expectations = fixture.cases.collectEntries { row ->
        [(row.id): [assembly: row.id, source: row.source,
                    source_records: row.original,
                    retained: row.expected[params.cleanup_scenario],
                    changed: row.expected[params.cleanup_scenario] != row.original]]
    }
    def downstream = CONTAMINATION_POST_FCS.out.assembly.map { meta, fasta ->
        tuple(meta, fasta, expectations[meta.id.toString()])
    }
    CHECK_CLEANUP_INPUT(downstream,
        file("${projectDir}/../scripts/comparisons/contamination_cleanup_acceptance.py").normalize())
    CHECK_CLEANUP_INPUT.out.toList().view { rows ->
        assert rows.size() == fixture.cases.size(): 'An assembly was lost or duplicated in cleanup routing'
        assert rows.collect { it[0].toString() }.toSet() == expectations.keySet()
        'All cleanup outputs reached the downstream consumer with the expected sequences'
    }
}
