nextflow.enable.dsl=2
params.pangenome_validate_only = false
params.pangenome_gref = true
params.pangenome_gref_source = 'clip'
params.pangenome_grefl = null
params.pangenome_gref_min_len = null
params.pangenome_cactus_lasttrain = true
params.pangenome_use_gpu = false
params.pangenome_cactus_extra = ''
params.pangenome_vcfbub = 100000
params.pangenome_viz = false
params.pangenome_chrom_og = false
params.mitohifi_perc_identity = 50
params.mitohifi_cov_cutoff = 'auto'
params.mitohifi_bloom_filter = false
include { PANGENOME_CONSTRUCTION } from '__REPO__/workflows/pangenome_construction.nf'
include { ORGANELLE } from '__REPO__/workflows/organelle.nf'

workflow {
    def root = file('/tmp/chimera-nxf-runtime/core-fixtures')
    root.mkdirs()
    def cohorts = ['1','2'].collect { taxid ->
        def names = ["REF_${taxid}", "IND_${taxid}"]
        def fastas = ['a','b'].collect { sample ->
            def fa = file("${root}/${taxid}_${sample}.fa")
            fa.text = '>chr1\nACGT\n'
            fa
        }
        def members = ['a','b'].withIndex().collect { sample, i ->
            [id:"${taxid}_${sample}",sample:"${taxid}_${sample}",taxid:taxid,
             haplotype:'primary', graph_name:names[i],fasta:fastas[i].toString()]
        }
        def cohort = [taxid:taxid, reference_id:"${taxid}_a",reference_name:names[0],
            reference_contigs:['chr1'], min_individuals:2, members:members]
        tuple(taxid,names[0],names,fastas,['chr1'],cohort)
    }
    if (params.contract_case == 'cross_species') cohorts[0][5].members[1].taxid = '999'
    PANGENOME_CONSTRUCTION(Channel.fromList(cohorts))
    PANGENOME_CONSTRUCTION.out.report.collect(flat: false).view { reports ->
        assert reports.size() == 4  // input AND build report, independently for each species
        'PASS: both species retain input and build reports'
    }
    PANGENOME_CONSTRUCTION.out.variant_input.collect(flat: false).view { inputs ->
        assert inputs.size() == 2
        assert inputs*.getAt(0).sort() == ['1','2']
        'PASS: validated graph roles and variant/identity handles for two species'
    }
    def reads = file("${root}/reads.fq")
    reads.text = '@read\nACGT\n+\nIIII\n'
    def gb = file("${root}/reference.gb")
    gb.text = 'reference tool double\n'
    ORGANELLE(Channel.of(
        tuple([sample:'one',taxid:'1',hifi:true,mito_tool:'mitohifi'],reads,reads,reads),
        tuple([sample:'skip_hifi',taxid:'1',hifi:true,mito_tool:'none'],reads,reads,reads),
        tuple([sample:'skip_sr',taxid:'2',hifi:false,mito_tool:'none'],reads,reads,reads)),
        Channel.of(tuple('1',cohorts[0][3][0],gb)), Channel.of(tuple('1',2)), Channel.empty(),
        [hifi:true,shortread_organelle:false])
    ORGANELLE.out.baits.collect(flat:false).view { baits ->
        assert baits*.getAt(0).sort() == ['one','skip_hifi','skip_sr']
        assert baits.every { it[1].isEmpty() }
        'PASS: organelle no-recovery and disabled HiFi/short-read routes release every nuclear sample'
    }
}
