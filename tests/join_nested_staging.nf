/* Generate real wrappers for the staging layout used by the production modules. */
nextflow.enable.dsl=2
params.fixture_source = null

process NESTED_INPUT {
    tag "${mode}"
    stageInMode { mode }
    executor 'local'
    cpus 1
    memory '1 GB'
    input:
    tuple val(mode), path(source,stageAs:'agp/*')
    output:
    tuple val(mode), path('result.txt')
    script:
    """
    # NESTED_STAGING_FIXTURE
    cat '${source}' > result.txt
    """
}

workflow {
    if (!params.fixture_source) error 'Supply --fixture_source pointing to a small test file'
    def source=file(params.fixture_source,checkIfExists:true)
    NESTED_INPUT(Channel.fromList(['symlink','copy','link']).map { mode -> tuple(mode,source) })
}
