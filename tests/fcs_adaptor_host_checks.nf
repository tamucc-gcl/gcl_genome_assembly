nextflow.enable.dsl=2
include { FCS_ADAPTOR as SCREEN } from '__REPO__/modules/fcs_adaptor.nf'
workflow {
    def root=file('/tmp/fcs-adaptor-host-check'); root.mkdirs()
    def fa=root.resolve('input.fasta');fa.text='>c1\nACGT\n'
    SCREEN(Channel.value(tuple([id:'fixture'],fa,'euk','singularity','contig')),
        file('__REPO__/py_scripts/run_fcs_adaptor.py'))
    SCREEN.out.provenance.view { m,p ->
        def manifest=new groovy.json.JsonSlurper().parse(p)
        assert manifest.runtime=='singularity fixture 1'
        assert manifest.command.contains('/app/fcs/bin/av_screen_x')
        assert !manifest.command.contains('/bin/bash')
        assert !manifest.command.contains('/bin/sh')
        assert manifest.command[-1]=='/sample-volume/input.fasta'
        'PASS: host task calls the application directly with native Singularity wrapping enabled globally'
    }
}
