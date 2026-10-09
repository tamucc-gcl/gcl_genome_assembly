/* Cohort misassembly review before finishing. One source-bound manual cut interface. */
include { MISASSEMBLY_CATALOG; MISASSEMBLY_ALIGN; MISASSEMBLY_DISCOVER; MISASSEMBLY_MAP_HIFI; MISASSEMBLY_ASSESS; MISASSEMBLY_REVIEW } from '../modules/misassembly.nf'
include { BREAK_CHIMERAS } from '../modules/break_chimeras.nf'
include { HARMONIZE_SPECIES as CHIMERA_REASSIGN_SPECIES } from '../modules/harmonize_species.nf'
include { harmonizerArgs } from './harmonize_scaffolds.nf'

workflow CHIMERA {
    take:
    ch_harmonized
    ch_shortread_finished
    ch_round1_agp
    ch_reference_id
    ch_peer_quality
    ch_contig_pairs_in
    ch_hic_evidence_inputs
    ch_native_graphs
    ch_telo_by_taxid
    ch_hifi_reads
    capabilities

    main:
    def absent = file("${projectDir}/assets/NO_PAIRS",checkIfExists:true)
    def noPaf = file("${projectDir}/assets/NO_PAF",checkIfExists:true)
    def detect_on = capabilities.harmonize && params.harmonize_scaffold_names && params.chimera_detect != false
    if (params.chimera_break?.toString() == 'auto') error 'Automatic cutting is deferred; supply a reviewed review.tsv.'
    def retainedBams = params.chimera_hifi_bam_manifest ? new groovy.json.JsonSlurper().parseText(file(params.chimera_hifi_bam_manifest,checkIfExists:true).text) : [:]
    def helpers = { names -> Channel.value(names.collect { file("${projectDir}/py_scripts/${it}",checkIfExists:true) }) }
    def core = ['misassembly_core.py','chimera_tracks.py','chimera_controls.py','chimera_intervals.py']
    ch_versions = Channel.empty()
    ch_broken = Channel.empty()
    ch_coordinate_report = Channel.value(noPaf)
    ch_registry_report = Channel.value(noPaf)
    ch_decision_report = Channel.value(noPaf)
    ch_evidence_report = Channel.value(noPaf)
    ch_review_input = Channel.empty()
    ch_action_review = Channel.value(absent)
    ch_cut_verifications = Channel.value([absent])
    ch_pre_finalize = ch_harmonized.mix(ch_shortread_finished.map { m,fa -> tuple(m,fa,file("${projectDir}/assets/NO_HARMONIZE",checkIfExists:true)) })

    if (detect_on) {
        ch_quality = ch_peer_quality.toList().map { files -> [quality:files.collectMany { f ->
            def lines=f.readLines().findAll { it && !it.startsWith('#') }; def header=lines[0].split('\t') as List
            lines.drop(1).collect { line -> def values=line.split('\t',-1)
                [(values[header.indexOf('id')]):[eligible:values[header.indexOf('role')]=='voter' && !values[header.indexOf('role_reasons')].toLowerCase().contains('forced'),reason:values[header.indexOf('role_reasons')]]] }
        }.collectEntries()] }
        ch_cohort = ch_harmonized.toList().map { records -> [records:records.sort { it[0].id }] }
            .combine(ch_quality)
            .map { cohort,quality -> [records:cohort.records.collect { m,fa,nm ->
                [meta:m,fasta:fa,name_map:nm,
                 catalog:[id:m.id.toString(),sample:m.sample.toString(),taxid:m.taxid.toString(),eligible:quality.quality[m.id.toString()]?.eligible == true,
                          role_reason:quality.quality[m.id.toString()]?.reason ?: 'quality unavailable']]
            }] }
        ch_pair_inputs = ch_cohort.flatMap { c ->
            c.records.groupBy { it.meta.taxid.toString() }.collectMany { taxid,records ->
                def pairs=[]
                for (int i=0;i<records.size();i++) for (int j=i+1;j<records.size();j++)
                    pairs << tuple([a:records[i].meta.id.toString(),b:records[j].meta.id.toString(),taxid:taxid,key:"${taxid}_pair_${i}_${j}"],records[i].fasta,records[j].fasta)
                pairs
            }
        }
        MISASSEMBLY_ALIGN(ch_pair_inputs)
        ch_alignments = MISASSEMBLY_ALIGN.out.alignment.toList().map { records -> [records:records.sort { it[0].key }] }
        MISASSEMBLY_CATALOG(ch_cohort.combine(ch_alignments).map { cohort,alignments ->
            tuple(cohort.records.collect { it.catalog },cohort.records.collect { it.name_map },alignments.records.collect { it[0] },alignments.records ? alignments.records.collect { it[1] } : [noPaf]) },
            helpers(core+['misassembly_discover.py']))
        ch_discovery_input = ch_harmonized.combine(ch_alignments).combine(MISASSEMBLY_CATALOG.out.catalog.first())
            .map { m,fa,nm,alignments,catalog ->
                def own=alignments.records.findAll { it[0].a==m.id.toString() || it[0].b==m.id.toString() }
                tuple(m,fa,catalog,own.collect { it[0] },own ? own.collect { it[1] } : [noPaf])
            }
        MISASSEMBLY_DISCOVER(ch_discovery_input,helpers(core+['misassembly_discover.py']))
        ch_reads = ch_hifi_reads.toList().map { records -> [reads:records.collectEntries { m,fq -> [(m.sample.toString()):fq] }] }
        ch_mapping_input = MISASSEMBLY_DISCOVER.out.discovery.join(ch_harmonized).combine(ch_reads)
            .map { m,discovery,fa,nm,reads ->
                def retained=params.chimera_hifi_context ? retainedBams[m.id.toString()] : null
                tuple(m,fa,discovery.resolve('discovery.json'),!retained && params.chimera_hifi_context ? reads.reads[m.sample.toString()] ?: absent : absent,
                    retained ? [file(retained.bam,checkIfExists:true),file(retained.index,checkIfExists:true)] : [absent],retained ? file(retained.provenance,checkIfExists:true) : absent)
            }
        MISASSEMBLY_MAP_HIFI(ch_mapping_input,helpers(core+['misassembly_map.py']))
        ch_assays = ch_hic_evidence_inputs.map { m,fa,libraries -> tuple(m.id.toString(),fa,libraries) }
            .join(ch_round1_agp.map { m,agp -> tuple(m.id.toString(),agp) })
            .join(ch_contig_pairs_in.map { m,stage,pairs -> tuple(m.id.toString(),pairs) }).toList()
            .map { records -> [assays:records.collectEntries { id,source,libraries,agp,pairs -> [(id):[source:source,libraries:libraries,agp:agp,pairs:pairs]] }] }
        ch_agps = ch_round1_agp.toList().map { records -> [agps:records.collectEntries { m,agp -> [(m.id.toString()):agp] }] }
        ch_graphs = ch_native_graphs.toList().map { records -> [graphs:records.collectEntries { m,graphs -> [(m.sample.toString()):graphs instanceof List ? graphs : [graphs]] }] }
        ch_motifs = ch_telo_by_taxid.toList().map { records -> [motifs:records.collectEntries { taxid,motif -> [(taxid.toString()):motif] }] }
        ch_assessment_input = MISASSEMBLY_DISCOVER.out.discovery.join(MISASSEMBLY_MAP_HIFI.out.mapping).join(ch_harmonized)
            .combine(ch_assays).combine(ch_agps).combine(ch_graphs).combine(ch_motifs)
            .map { m,discovery,mapping,fa,nm,assays,agps,graphs,motifs ->
                def assay=assays.assays[m.id.toString()] ?: [:]
                def suffix=m.id.toString().startsWith(m.sample.toString()+'_') ? m.id.toString().substring(m.sample.toString().length()+1) : 'primary'
                def graph=graphs.graphs[m.sample.toString()]?.find { it.name==m.sample.toString()+'.'+suffix+'.p_ctg.gfa' }
                tuple(m,fa,discovery.resolve('discovery.json'),mapping,agps.agps[m.id.toString()] ?: absent,
                      assay.pairs ?: absent,assay.source ?: absent,assay.libraries ?: absent,graph ?: absent,motifs.motifs[m.taxid.toString()] ?: 'CCCTAA')
            }
        MISASSEMBLY_ASSESS(ch_assessment_input,helpers(core+['misassembly_assess.py','chimera_graph_evidence.py']))
        ch_review_input = MISASSEMBLY_ASSESS.out.evidence.toList().map { packets -> [packets:packets] }.combine(ch_cohort).map { bundle,cohort ->
            def packets=bundle.packets
            def expected=cohort.records.collect { it.meta.id.toString() }.sort();def observed=packets.collect { it[0].id.toString() }.sort()
            if(expected!=observed) error "Misassembly evidence incomplete: expected ${expected}; received ${observed}"
            packets.sort { it[0].id }.collect { it[1] }
        }
    }
    if( capabilities.harmonize && params.harmonize_scaffold_names && params.chimera_break && params.chimera_break.toString() != 'false' ) {
        ch_break_script = Channel.value(['break_chimeras.py', 'chimera_actions.py']
            .collect { file("${projectDir}/py_scripts/${it}", checkIfExists: true) })
        // Review file is a staged input: changing it invalidates cutting/downstream tasks only.
        ch_known_assemblies = ch_harmonized.map { m, fa, nm -> m.id.toString() }.toList().map { ids -> [known:ids] }
        ch_review_file = Channel.fromPath(params.chimera_break.toString(), checkIfExists: true).first()
            .combine(ch_known_assemblies)
            .map { review, catalog ->
                def lines = review.readLines().findAll { it && !it.startsWith('#') }
                def header = lines[0].split('\t') as List
                if (!header.contains('selected') || !header.contains('assembly')) error 'Review file requires selected and assembly columns'
                lines.drop(1).each { line ->
                    def values = line.split('\t',-1)
                    def selection = values[header.indexOf('selected')]
                    if (!(selection in ['YES','NO'])) error 'Every review row requires selected=YES or NO'
                    if (selection=='YES' && !(values[header.indexOf('assembly')] in catalog.known))
                        error 'Selected review row names an unknown assembly: '+values[header.indexOf('assembly')]
                }
                review }
        ch_break_in = ch_harmonized.combine(ch_review_file)
        ch_break_in.branch { meta, fa, nm, cand ->
            def lines = cand.readLines().findAll { it && !it.startsWith('#') }
            def header = lines[0].split('\t') as List
            def selected = lines.drop(1).any { line ->
                def row = line.split('\t',-1)
                row[header.indexOf('selected')]=='YES' && row[header.indexOf('assembly')]==meta.id.toString() }
            cut: nm.name != 'NO_HARMONIZE' && selected
            passthrough: true
        }.set { ch_break_routes }
        BREAK_CHIMERAS(ch_break_routes.cut, ch_break_script)
        ch_action_review = ch_review_file
        ch_cut_verifications = BREAK_CHIMERAS.out.verification.toList().map { records -> records ? records.collect { it[1] } : [absent] }
        ch_versions = ch_versions.mix(BREAK_CHIMERAS.out.versions)
        ch_broken = BREAK_CHIMERAS.out.assemblies.map { m,fa,nm -> tuple(m.id.toString(),m,fa,nm) }
            .join(BREAK_CHIMERAS.out.verification.map { m,verification -> tuple(m.id.toString(),verification) })
            .filter { id,m,fa,nm,verification -> new groovy.json.JsonSlurper().parseText(verification.text).selected_actions.size() > 0 }
            .map { id,m,fa,nm,verification -> tuple(m,fa,nm) }
        ch_pre_finalize = BREAK_CHIMERAS.out.assemblies
            .mix(ch_break_routes.passthrough.map { meta, fa, nm, cand -> tuple(meta, fa, nm) })
            .mix(ch_shortread_finished.map { meta, fa -> tuple(meta, fa, file("${projectDir}/assets/NO_HARMONIZE", checkIfExists: true)) })
        // Reassign from actual corrected sequences, never from parent side labels.
        ch_reassign_lr = ch_pre_finalize.filter { m, fa, nm -> m.assembler == 'hifiasm' }
        ch_changed_taxids = BREAK_CHIMERAS.out.verification
            .filter { m, verification -> new groovy.json.JsonSlurper().parseText(verification.text).selected_actions.size() > 0 }
            .map { m, verification -> m.taxid.toString() }.toList().map { ids -> [taxids:ids] }
        ch_reassign_species = ch_reassign_lr.map { m, fa, nm -> tuple(m.taxid.toString(), m.id.toString(), fa) }
            .groupTuple()
            .filter { taxid, ids, fastas -> ids.size() >= 2 }
            .map { taxid, ids, fastas ->
                def order = (0..<ids.size()).toList().sort { ids[it] }
                tuple(taxid, order.collect { ids[it] }, order.collect { fastas[it] }) }
            .combine(ch_changed_taxids)
            .filter { taxid, ids, fastas, changed -> taxid in changed.taxids }
            .map { taxid, ids, fastas, changed -> tuple(taxid, ids, fastas) }
            .join(ch_reference_id.map { taxid, id -> tuple(taxid.toString(), id) })
        CHIMERA_REASSIGN_SPECIES(ch_reassign_species,
            file("${projectDir}/py_scripts/harmonize_names.py", checkIfExists: true), harmonizerArgs())
        ch_versions = ch_versions.mix(CHIMERA_REASSIGN_SPECIES.out.versions)
        ch_reassigned_maps = CHIMERA_REASSIGN_SPECIES.out.name_maps
            .flatMap { taxid, maps ->
                (maps instanceof List ? maps : [maps]).collect { nm ->
                    tuple(nm.name.replaceFirst(/\.harmonized_name_map\.tsv$/, ''), nm) } }
        ch_pre_finalize = ch_reassign_lr.map { m, fa, nm -> tuple(m.id.toString(), m, fa, nm) }
            .join(ch_reassigned_maps, remainder: true)
            .filter { row -> row[1] != null }
            .map { id, m, fa, oldMap, newMap -> tuple(m, fa, newMap ?: oldMap) }
            .mix(ch_shortread_finished.map { m, fa -> tuple(m, fa, file("${projectDir}/assets/NO_HARMONIZE", checkIfExists: true)) })
    }


    if (detect_on) {
        MISASSEMBLY_REVIEW(ch_review_input,helpers(core+['misassembly_report.py']),ch_action_review,ch_cut_verifications)
        ch_coordinate_report = MISASSEMBLY_REVIEW.out.coordinates
        ch_registry_report = MISASSEMBLY_REVIEW.out.registry
        ch_decision_report = MISASSEMBLY_REVIEW.out.decisions
        ch_evidence_report = MISASSEMBLY_REVIEW.out.report
        ch_versions = ch_versions.mix(MISASSEMBLY_REVIEW.out.versions)
    }

    emit:
    coordinate_report = ch_coordinate_report
    pre_finalize = ch_pre_finalize
    broken = ch_broken
    registry = ch_registry_report
    decisions = ch_decision_report
    report = ch_evidence_report
    versions = ch_versions
}
