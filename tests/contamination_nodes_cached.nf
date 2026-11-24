nextflow.enable.dsl = 2

params.fixture_root = '/tmp/contamination-nodes-cached'
params.contamination_megablast_nodes = 4
params.contamination_megablast_slots = 1
params.contamination_megablast_threads = 1
params.contamination_megablast_scratch = '/tmp'
params.contamination_megablast_time = '1h'
params.contamination_megablast_partition = 'normal'
params.contamination_megablast_nodelist = ''
params.outdir = '/tmp/contamination-nodes-cached/results'
include { CONTAMINATION_MEGABLAST_NODES } from '__REPO__/workflows/contamination_megablast_nodes.nf'

workflow {
    def root = file(params.fixture_root)
    CONTAMINATION_MEGABLAST_NODES(Channel.of(tuple('A', [], root.resolve('A')), tuple('B', [], root.resolve('B'))),
        root.resolve('db/nt').toString(), file('__REPO__'))
    CONTAMINATION_MEGABLAST_NODES.out.results.toList().view { rows ->
        assert rows.size() == 2
        'Completed cached chunks collected across two assemblies without allocating worker nodes'
    }
}
