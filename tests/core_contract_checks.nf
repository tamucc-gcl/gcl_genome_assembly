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
include { PANGENOME_STATS } from '__REPO__/modules/pangenome_stats.nf'
include { PANGENOME_BUILD_QC } from '__REPO__/modules/pangenome_build_qc.nf'
include { PANGENOME_SOURCE_AUDIT } from '__REPO__/modules/pangenome_source_audit.nf'
include { PANGENOME_REGIONAL_CONTEXT } from '__REPO__/modules/pangenome_regional_plots.nf'

workflow {
    def root = file('/tmp/chimera-nxf-runtime/core-fixtures')
    root.mkdirs()
    def cohorts = ['1','2'].collect { taxid ->
        def names = ["REF_${taxid}", "IND_${taxid}"]
        def fastas = ['a','b'].collect { sample ->
            def fa = file("${root}/${taxid}_${sample}.fa")
            if (System.getenv('CORE_RESUME_CHECK') != 'true' || !fa.exists()) fa.text = '>chr1\nACGT\n'
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
    if (params.pangenome_assignment_replay) {
        PANGENOME_CONSTRUCTION.out.assignment_input.view { taxid,graph,paf,gaf,contigs ->
            assert graph.name == taxid+'.sv.gfa.gz'
            assert paf.name == taxid+'.paf'
            assert gaf.name == taxid+'.gaf.gz'
            assert contigs == ['chr1']
            'PASS: cached assignment inputs retain species-specific graph, PAF and reference contigs'
        }
    }
    PANGENOME_STATS(PANGENOME_CONSTRUCTION.out.clip_stats_input)
    PANGENOME_BUILD_QC(PANGENOME_CONSTRUCTION.out.support_qc
        .join(PANGENOME_STATS.out.vg_stats).join(PANGENOME_CONSTRUCTION.out.identities)
        .map { taxid,bundle,full,versions,stats,identities -> tuple(taxid,bundle,full,stats,identities,versions) },
        file('__REPO__/py_scripts/pangenome_build_qc.py'))
    PANGENOME_BUILD_QC.out.report.collect(flat:false).view { reports ->
        assert reports.size() == 2
        assert reports.every { it[1].text.contains('FULL support: 4 graph bp') }
        'PASS: per-species post-build clipping/provenance joins and outputs'
    }
    def sourceIndexes = cohorts.collect { row ->
        def indexes = row[5].members.collect { member ->
            def fai = file("${member.fasta}.fai")
            if (System.getenv('CORE_RESUME_CHECK') != 'true' || !fai.exists()) fai.text = 'chr1\t4\n'
            fai
        }
        def specs = row[5].members.withIndex().collect { member,i ->
            [id:member.id, fai:indexes[i].name, chromosomes:['chr1']]
        }
        tuple(row[0],indexes,specs)
    }
    PANGENOME_SOURCE_AUDIT(PANGENOME_BUILD_QC.out.files
        .map { taxid,files -> tuple(taxid,files.find { it.name == 'full_clip_sequence_paths.tsv' }) }
        .join(PANGENOME_CONSTRUCTION.out.identities).join(Channel.fromList(sourceIndexes))
        .join(PANGENOME_CONSTRUCTION.out.source_evidence)
        .map { taxid,paths,identities,indexes,specs,mappings,evidence ->
            tuple(taxid,identities,paths,indexes,specs,mappings,evidence)
        }, file('__REPO__/py_scripts/pangenome_source_audit.py'))
    PANGENOME_SOURCE_AUDIT.out.report.collect(flat:false).view { reports ->
        assert reports.size() == 2
        assert reports.every { it[1].text.contains('Source sequence representation') }
        assert reports.every { it[1].text.contains('Status: **PASS**') }
        'PASS: per-species source dictionaries, cached mappings and graph path joins'
    }
    def eligibility = file("${root}/eligibility.json")
    eligibility.text = groovy.json.JsonOutput.toJson([assemblies: cohorts.collectMany { row ->
        row[5].members.collect { m -> [id:m.id,taxid:row[0],chromosome_scale:true,reference_contigs:['chr1_1']] }
    }])
    def regional = cohorts.collect { row ->
        def sharing = file("${root}/${row[0]}.sharing.tsv")
        sharing.text = 'scope\thaplotype\tchromosome\n' + row[5].members.collect { m ->
            "haplotype\t${m.id}\t\nchromosome_haplotype\t${m.id}\tchr1\n"
        }.join('')
        tuple(row[0],sharing,eligibility)
    }
    PANGENOME_REGIONAL_CONTEXT(Channel.fromList(regional),file('__REPO__/py_scripts/pangenome_regional_context.py'))
    PANGENOME_REGIONAL_CONTEXT.out.context.collect(flat:false).view { records ->
        assert records.size() == 2
        assert records.every { it[2].text.contains('separate_group') }
        'PASS: production regional context process for two species'
    }
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
