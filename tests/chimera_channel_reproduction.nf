nextflow.enable.dsl=2
workflow {
    def seed_meta=[id:'sample_hap1',sample:'sample',taxid:'373251']
    ch_peer_quality=Channel.of(file('/tmp/chimera-nxf-runtime/quality.tsv'))
    ch_harmonized=Channel.of(tuple(seed_meta,file('/tmp/chimera-nxf-runtime/a.fa'),file('/tmp/chimera-nxf-runtime/names.tsv')))
    ch_hifi_reads=Channel.of(tuple([sample:'sample'],file('/tmp/chimera-nxf-runtime/reads.fq')))
    ch_native_graphs=Channel.of(tuple([sample:'sample'],[file('/tmp/chimera-nxf-runtime/sample.hap1.p_ctg.gfa')]))
    ch_hic_evidence_inputs=Channel.of(tuple(seed_meta,file('/tmp/chimera-nxf-runtime/a.fa'),file('/tmp/chimera-nxf-runtime/lib.tsv')))
    ch_round1_agp=Channel.of(tuple(seed_meta,file('/tmp/chimera-nxf-runtime/a.agp')))
    ch_contig_pairs_in=Channel.of(tuple(seed_meta,'scaffold',file('/tmp/chimera-nxf-runtime/a.pairs')))
    ch_chimeric_joins=Channel.of(tuple('373251','sample_hap1',file('/tmp/chimera-nxf-runtime/calls.tsv')))
    ch_telo_by_taxid=Channel.of(tuple('373251','CCCTAA'))
    def retainedBams=[:]
    params.chimera_hifi_context=true
            ch_context_quality = ch_peer_quality.toList().map { files ->
                log.info "[CHIMERA INPUT] Comparison quality ready: ${files.size()} tables"
                files.collectMany { quality ->
                    def lines = quality.readLines().findAll { it && !it.startsWith('#') }
                    def header = lines[0].split('\t') as List
                    lines.drop(1).collect { line ->
                        def fields = line.split('\t', -1)
                        [(fields[header.indexOf('id')]): fields[header.indexOf('role')] == 'voter' &&
                            !fields[header.indexOf('role_reasons')].toLowerCase().contains('forced')] }
                }.collectEntries() }
            ch_context_cohort = ch_harmonized.toList().map { records ->
                log.info "[CHIMERA INPUT] Assessment cohort ready: ${records.collect { it[0].id }}"
                [records:records] }
                .combine(ch_context_quality)
                .map { cohort, quality -> [records: cohort.records.collect { m, fa, nm ->
                    def labels = [:]
                    if (nm.name != 'NO_HARMONIZE') {
                        def lines = nm.readLines().findAll { it && !it.startsWith('#') }
                        def header = lines[0].split('\t') as List
                        lines.drop(1).each { line ->
                            def fields = line.split('\t', -1)
                            def name = fields[header.indexOf('new_name')]
                            if (header.contains('chromosome_member') && fields[header.indexOf('chromosome_member')]=='yes' && name.startsWith('chr') && !name.contains('+'))
                                labels[fields[header.indexOf('old_name')]] = name.split('_')[0]
                        }
                    }
                    tuple(m + [auto_evidence: quality[m.id.toString()] == true, chromosome_labels:labels], fa, nm) }] }
            ch_context_reads = ch_hifi_reads.toList()
                .map { records ->
                    log.info "[CHIMERA INPUT] HiFi reads ready: ${records.collect { it[0].sample }}"
                    records.collectEntries { m, fq -> [(m.sample.toString()): fq] } }
            ch_context_graphs = ch_native_graphs.toList().map { records ->
                log.info "[CHIMERA INPUT] Native graphs ready: ${records.collect { it[0].sample }}"
                records.collectEntries { m, paths -> [(m.sample.toString()): paths instanceof List ? paths : [paths]] } }
            ch_context_assays = ch_hic_evidence_inputs
                .map { m, fa, libraries -> tuple(m.id.toString(), fa, libraries) }
                .join(ch_round1_agp.map { m, agp -> tuple(m.id.toString(), agp) })
                .join(ch_contig_pairs_in.map { m, stage, pairs -> tuple(m.id.toString(), pairs) })
                .toList()
                .map { records ->
                    log.info "[CHIMERA INPUT] Hi-C assay collection closed: ${records.collect { it[0] }}"
                    records.collectEntries { id, source, libraries, agp, pairs ->
                        [(id): [source:source, libraries:libraries, agp:agp, pairs:pairs]] } }
            ch_context_in = ch_chimeric_joins
                .map { taxid, id, calls ->
                    log.info "[CHIMERA INPUT] Candidate calls received: ${id}"
                    tuple(id.toString(), calls) }
                .join(ch_harmonized.map { m, fa, nm -> tuple(m.id.toString(), m, fa) }, remainder: true)
                .filter { record -> record[1] != null }
                .map { id, calls, m, fa ->
                    if (m == null || fa == null) error "Chimera candidate calls have no matching assessed assembly: ${id}"
                    tuple(id, calls, m, fa) }
                .combine(ch_context_cohort)
                .combine(ch_context_reads)
                .combine(ch_context_assays)
                .combine(ch_context_graphs)
                .map { id, calls, m, fa, cohort, reads, assays, graphs ->
                    def peers = cohort.records.findAll { pm, pf, pn -> pm.taxid.toString() == m.taxid.toString() && pm.id != m.id }
                    def assay = assays[id] ?: [:]
                    def absent = file('/tmp/chimera-nxf-runtime/NO_PAIRS')
                    def retained = retainedBams[id]
                    def suffix = m.id.toString().startsWith(m.sample.toString()+'_') ? m.id.toString().substring(m.sample.toString().length()+1) : 'primary'
                    def nativeGraph = graphs[m.sample.toString()]?.find { it.name == m.sample.toString()+'.'+suffix+'.p_ctg.gfa' }
                    tuple(m.taxid.toString(), m, fa, calls,
                          peers.collect { it[0] }, peers ? peers.collect { it[1] } : [file('/tmp/chimera-nxf-runtime/NO_PAIRS')],
                          !retained && params.chimera_hifi_context && reads[m.sample.toString()] ? reads[m.sample.toString()] : absent,
                          assay.agp ?: absent, assay.pairs ?: absent, assay.source ?: absent, assay.libraries ?: absent,
                          retained ? [file(retained.bam,checkIfExists:true),file(retained.index,checkIfExists:true)] : [absent],
                          retained ? file(retained.provenance,checkIfExists:true) : absent, nativeGraph ?: absent)
                }
                .combine(ch_telo_by_taxid.map { taxid, motif ->
                    log.info "[CHIMERA INPUT] Motif received for taxid ${taxid}"
                    tuple(taxid.toString(), motif) }, by: 0)
                .map { taxid, m, fa, calls, pm, pf, reads, agp, pairs, source, libraries, bams, provenance, graph, motif ->
                    log.info "[CHIMERA INPUT] Sequence-context task ready: ${m.id}"
                    tuple(m, fa, calls, pm, pf, reads, motif ?: 'CCCTAA', agp, pairs, source, libraries, bams, provenance, graph) }
                .ifEmpty { error 'Chimera sequence context received no assessment inputs: inspect [CHIMERA INPUT] readiness messages and cohort/taxid joins' }
    ch_context_in.view { row -> 'READY '+row[0].id }
}
