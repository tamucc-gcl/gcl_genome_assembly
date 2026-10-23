nextflow.enable.dsl=2
include { PANGENOME_DELETION_DIAGNOSTIC } from '__REPO__/modules/pangenome_deletion_replay.nf'
workflow {
    def root=file('/tmp/chimera-nxf-runtime/deletion-fixtures')
    def adapters=['pangenome_deletion_replay.py','pangenome_assignment_replay.py','pangenome_source_audit.py']
        .collect { file('__REPO__/py_scripts/'+it) }
    PANGENOME_DELETION_DIAGNOSTIC(Channel.of(tuple('1',root.resolve('sv.gfa.gz'),root.resolve('input.paf'),
        root.resolve('input.paf.unfiltered.gz'),root.resolve('input.paf.filter.log'),
        ['chrA','chrB'],root.resolve('identities.tsv'),root.resolve('sources.tsv'))),adapters,
        file('__REPO__/vendor/cactus-gfa-tools-regional'),file('__REPO__/py_scripts/pangenome_regional_patch.py'))
    PANGENOME_DELETION_DIAGNOSTIC.out.report.view { taxid,report ->
        if (report.name == 'regional_gate_report.md') {
            assert report.text.contains('byte-identical')
            assert report.text.contains('chrA') && report.text.contains('chrB')
            return 'PASS: regional gate evidence report'
        }
        assert report.text.contains('chromosome_isolated') || report.text.contains('Corrected regional')
        if (!report.text.contains('Corrected regional')) assert report.text.contains('byte for byte')
        assert report.text.contains('chrA,chrB')
        'PASS: production deletion workflow, three policies, six assignments and summary'
    }
}
