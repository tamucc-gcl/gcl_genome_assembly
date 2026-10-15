nextflow.enable.dsl=2
include { FINAL_HIC_MAPS } from '__REPO__/workflows/final_hic_maps.nf'
workflow READ_PREPARATION {
    main:
    def root='/tmp/chimera-nxf-runtime'
    reads=Channel.of(
        tuple([sample:'single'],[file(root+'/one_R1.fq')],[file(root+'/one_R2.fq')]),
        tuple([sample:'multi'],[file(root+'/one_R1.fq'),file(root+'/two_R1.fq')],
                              [file(root+'/one_R2.fq'),file(root+'/two_R2.fq')]))
    emit:
    hic=reads
}
workflow {
    def root='/tmp/chimera-nxf-runtime'
    READ_PREPARATION()
    def post_on=false
    def qc_on=params.qc_mode != 'none'
    def capabilities=[hic:true]
    def assemblies=['single','multi','nohic'].collectMany { sample ->
        ['hap1','hap2'].collect { hap ->
            def sets=sample=='multi' ? [[library_id:'A',readset_id:'A1'],[library_id:'B',readset_id:'B1']] :
                                      [[library_id:'A',readset_id:'A1']]
            tuple([id:"${sample}_${hap}",sample:sample,haplotype:hap,hic:sample!='nohic',hic_readsets:sets],
                  file(root+(hap=='hap1' ? '/a.fa' : '/b.fa')))
        }
    }
    ch_finalized_assembly=Channel.fromList(assemblies)
    ch_compartments_script=file(root+'/names.tsv')
    ch_tad_book_script=file(root+'/names.tsv')
    ch_all_bam_metrics=Channel.empty()
    ch_all_pairs_metrics=Channel.empty()
    ch_versions=Channel.empty()
__FINAL_HIC_QC__
    ch_all_bam_metrics.toList().view { rows ->
        assert rows.size()==(qc_on ? 4 : 0)
        assert rows.every { it[1]=='final_raw_map' && it[0].sample!='nohic' }
        'PASS: final Hi-C mapping metrics follow QC independently of post-assembly toggle'
    }
    ch_all_pairs_metrics.toList().view { rows ->
        assert rows.size()==(qc_on ? 4 : 0)
        assert rows.every { it[1]=='final_filtered' }
        'PASS: final filtered-pair metrics include both haplotypes with one/two libraries'
    }
    ch_final_contact_maps.toList().view { rows ->
        assert rows.size()==(qc_on && params.run_final_contact_maps ? 4 : 0)
        assert rows.every { it[1]=='final' }
        'PASS: maps share final evidence and honor the plotting switch'
    }
    ch_versions.toList().view { rows ->
        assert rows.size()==(qc_on ? (params.run_final_contact_maps ? 12 : 8) : 0)
        'PASS: final Hi-C software versions reach the main workflow'
    }
}
