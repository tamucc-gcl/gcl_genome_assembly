nextflow.enable.dsl=2
params.pool_test_node='crest-c018'
params.pool_test_root=null
process POOL_ACCEPTANCE {
    executor 'local'
    conda 'conda-forge::python=3.11 bioconda::blast=2.16.0'
    cpus 1
    memory '1 GB'
    scratch false
    cache false
    errorStrategy 'terminate'
    time '1h'
    publishDir params.outdir, mode:'copy'
    input:
    path helpers
    output:
    path 'acceptance.json'
    script:
    def q = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    """
    python3 contamination_pool_acceptance.py --root ${q(params.pool_test_root)} \
        --node ${q(params.pool_test_node)} --worker-helper contamination_workers.py
    cp ${q(params.pool_test_root + '/acceptance.json')} acceptance.json
    """
}
workflow {
    if(!params.pool_test_root) error 'Supply --pool_test_root with a new small fixture directory'
    def base=file("${projectDir}/../..")
    POOL_ACCEPTANCE(['tests/contamination_pool_acceptance.py','tests/contamination_node_fixture.py',
        'py_scripts/contamination_workers.py','py_scripts/contamination_slurm_pool.py'].collect { file(base.resolve(it),checkIfExists:true) })
}
