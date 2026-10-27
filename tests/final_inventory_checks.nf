nextflow.enable.dsl=2
params.run_pangenome=false
params.qc_mode='none'
include { ASSEMBLY_ELIGIBILITY } from '__REPO__/modules/assembly_eligibility.nf'
include { SHARED_REFERENCE } from '__REPO__/modules/shared_reference.nf'
include { ASSEMBLY_RUN_SUMMARY } from '__REPO__/modules/assembly_run_summary.nf'

workflow FINALIZE_ASSEMBLY {
    take:
    records
    emit:
    assembly=records.map { meta,fa,fai -> tuple(meta,fa) }
    fai=records.map { meta,fa,fai -> tuple(meta,fai) }
}
workflow HARMONIZE_SCAFFOLDS {
    take:
    references
    emit:
    reference_id=references
}
workflow CHIMERA {
    main:
    absent=file('__REPO__/assets/NO_PAF')
    emit:
    coordinate_report=Channel.value(absent)
    report=Channel.value(absent)
}
workflow {
    def root=file('/tmp/chimera-nxf-runtime/inventory-fixtures');root.mkdirs()
    def pinned=System.getenv('CORE_INPUT_CASE')=='inventory_pinned_reference'
    def fixtureMeta=[sample:'sample',id:'sample_primary',taxid:'1',haplotype:'primary',
              assembler:(pinned ? 'hifiasm' : 'spades'),ploidy:(pinned ? 1 : 2),n_hap:1,hic:false]
    def fixtureFa=file("${root}/sample_primary.fasta");fixtureFa.text='>contig_1\nACGT\n'
    def fixtureFai=file("${root}/sample_primary.fasta.fai");fixtureFai.text='contig_1\t4\t10\t4\t5\n'
    FINALIZE_ASSEMBLY(Channel.of(tuple(fixtureMeta,fixtureFa,fixtureFai)))
    HARMONIZE_SCAFFOLDS(pinned ? Channel.of(tuple('1','sample_primary')) : Channel.empty())
    CHIMERA()
    ch_finalized_assembly=FINALIZE_ASSEMBLY.out.assembly
    ch_declared_assemblies=Channel.of(tuple('sample',fixtureMeta))
    ch_sample_identity=Channel.of(tuple('sample',[taxid:'1',name:'Species one']))
    ch_input_status=Channel.value([[sample:'sample',status:'accepted',expected:1,reason:'fixture']])
    __FINAL_INVENTORY__
    ASSEMBLY_ELIGIBILITY.out.manifest.toList().view { files ->
        assert files.size()==1
        def data=new groovy.json.JsonSlurper().parse(files[0])
        assert data.assemblies.size()==1
        assert data.linear_references*.id==['sample_primary']
        assert data.assemblies[0].representation==(pinned ? 'haploid' : 'collapsed_or_primary')
        'PASS: actual final inventory wiring with '+(pinned ? 'pinned' : 'empty')+' harmonization reference'
    }
    SHARED_REFERENCE.out.manifest.toList().view { rows ->
        assert rows.size()==1
        assert new groovy.json.JsonSlurper().parse(rows[0][1]).bases==4
        'PASS: final shared reference produced'
    }
    ASSEMBLY_RUN_SUMMARY.out.report.toList().view { files ->
        assert files.size()==1 && files[0].text.contains('| sample | finalized | 1/1 |')
        'PASS: final run summary produced'
    }
}
