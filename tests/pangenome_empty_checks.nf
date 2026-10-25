nextflow.enable.dsl=2
params.run_pangenome = (System.getenv('CORE_INPUT_CASE') != 'graph_disabled')
params.pangenome_validate_only = false
include { PANGENOME } from '__REPO__/workflows/pangenome.nf'

workflow {
    def manifest = file('/tmp/chimera-nxf-runtime/empty-eligibility.json')
    manifest.text = '{"cohorts":[],"assemblies":[],"settings":{"min_individuals":2}}\n'
    // Exercise the actual top-level workflow, not a reconstruction of its routing.
    PANGENOME(Channel.value(manifest))
    PANGENOME.out.report.toList().view { reports ->
        assert reports.isEmpty()
        'PASS: disabled or unadmitted graph emits no reports and closes cleanly'
    }
}
