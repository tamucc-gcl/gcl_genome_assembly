nextflow.enable.dsl=2
params.manifest=null
params.repo=null
params.chimera_plots=true
include { MISASSEMBLY_JUNCTION } from '../modules/misassembly.nf'
workflow {
    def jobs=new groovy.json.JsonSlurper().parseText(file(params.manifest).text).jobs
    def inputs=Channel.fromList(jobs).map { j ->
        tuple([id:j.assembly,sample:j.sample],file(j.packet),file(j.paths.fasta),file(j.mapping_dir),
            file(j.paths.agp),file(j.paths.pairs),file(j.paths.source),file(j.paths.libraries),file(j.catalog),
            j.alignments.collect { [a:it.a,b:it.b] },j.alignments.collect { file(it.path) },
            j.competitive_fastas.keySet().toList(),j.competitive_fastas.values().collect { file(it) })
    }
    def helpers=Channel.value(['misassembly_core.py','chimera_tracks.py','chimera_controls.py','chimera_intervals.py','misassembly_junction.py','misassembly_enrich.py','misassembly_plot.py'].collect { file(params.repo+'/py_scripts/'+it) })
    MISASSEMBLY_JUNCTION(inputs,helpers)
    MISASSEMBLY_JUNCTION.out.evidence.view { m,p -> "JUNCTION_OUTPUT ${m.id} ${p}" }
}
