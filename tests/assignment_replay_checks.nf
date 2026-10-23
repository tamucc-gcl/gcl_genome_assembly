nextflow.enable.dsl=2
include { PANGENOME_ASSIGNMENT_PREPARE; PANGENOME_ASSIGNMENT_REPLAY; PANGENOME_ASSIGNMENT_REPORT } from '__REPO__/modules/pangenome_assignment_replay.nf'
workflow {
    def root=file('/tmp/chimera-nxf-runtime/assignment-fixtures')
    def adapters=['pangenome_assignment_replay.py','pangenome_source_audit.py'].collect { file('__REPO__/py_scripts/'+it) }
    PANGENOME_ASSIGNMENT_PREPARE(Channel.of(tuple('1',root.resolve('sv.gfa.gz'),root.resolve('input.paf'),
        root.resolve('input.gaf.gz'), ['chrA','chrB'],root.resolve('identities.tsv'),root.resolve('sources.tsv'))),adapters)
    PANGENOME_ASSIGNMENT_REPLAY(PANGENOME_ASSIGNMENT_PREPARE.out.prepared.flatMap { taxid,graph,paf,config,selection,identities,provenance,contigs ->
        ['baseline','whole_relaxed','regional_1Mb','regional_3Mb'].collect { name -> tuple(taxid,graph,paf,config,selection,identities,provenance,contigs,name) }
    },adapters)
    PANGENOME_ASSIGNMENT_REPORT(PANGENOME_ASSIGNMENT_REPLAY.out.results.groupTuple()
        .map { taxid,names,dirs -> tuple(taxid,dirs) }
        .join(PANGENOME_ASSIGNMENT_PREPARE.out.coverage),adapters)
    PANGENOME_ASSIGNMENT_REPORT.out.report.view { taxid,report ->
        assert report.text.contains('regional_3Mb')
        assert report.text.contains('baseline')
        assert report.text.contains('coverage_changes.tsv')
        'PASS: four cached-data assignment cases and comparison report'
    }
}
