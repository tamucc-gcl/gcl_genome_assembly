/* Validate both tables before scheduling samples. Nextflow 23.10 DSL2.
 * Relative paths resolve against the containing table. No duplicated sample rows.
 */
include { buildMeta } from './meta.nf'

def inputTable(tablePath, required) {
    def sheet = file(tablePath, checkIfExists: true)
    // Nextflow's own CSV parser, evaluated before registering workflow calls.
    // Synchronous list API verified against Nextflow 23.10.1.
    def records = new nextflow.splitter.CsvSplitter()
        .target(sheet)
        .options(header: false, sep: ',', quote: '"', strip: true)
        .list()
        .findAll { cells -> cells.any { it?.toString()?.trim() } }
    if (!records) error "Empty input table: ${sheet}"
    def header = records[0].collect { it?.toString()?.trim()?.replace('\uFEFF', '') }
    if (header.any { !it } || header.unique(false).size() != header.size())
        error "Empty or duplicate column names in ${sheet}"
    if (!header.containsAll(required)) error "${sheet} requires columns: ${required.join(', ')}"
    def rows = records.drop(1).withIndex().collect { cells, index ->
        if (cells.size() != header.size())
            error "${sheet}: data record ${index + 1} has ${cells.size()} fields; expected ${header.size()}"
        def row = [:]
        header.eachWithIndex { key, i -> row[key] = cells[i]?.toString()?.trim() ?: '' }
        row
    }
    return [sheet: sheet, rows: rows]
}

def parseSampleSheet(sampleSheetPath, hicReadsetsPath = null) {
    def sampleTable = inputTable(sampleSheetPath, ['sample_id'])
    def libraryTable = hicReadsetsPath
        ? inputTable(hicReadsetsPath, ['sample_id', 'library_id', 'readset_id', 'hic_r1', 'hic_r2'])
        : [sheet: null, rows: []]

    def table = sampleTable
    def libraries = libraryTable
    def safeId = { value, label ->
        if (!(value ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/))
            error "Invalid ${label} '${value}': use a letter/number followed by letters, numbers, underscore, dot or hyphen"
    }
    def resolve = { value, sheet ->
        if (!value) return null
        def absolute = value.startsWith('/') || value ==~ /^[A-Za-z]:[\\\/].*/ || value.contains('://')
        absolute ? file(value) : file("${sheet.parent}/${value}")
    }
    def seenSamples = [] as Set
    table.rows.each { row ->
        safeId(row.sample_id, 'sample_id')
        if (!seenSamples.add(row.sample_id))
            error "Duplicate sample_id '${row.sample_id}': use --hic_readsets for additional libraries"
    }
    def seenReadsets = [] as Set
    def seenOutputIds = [] as Set
    def seenPaths = [:]
    def checkPath = { path, owner ->
        if (path == null) return
        def key = path.toAbsolutePath().normalize().toString()
        if (seenPaths.containsKey(key)) error "Read file assigned twice: ${path} (${seenPaths[key]}; ${owner})"
        seenPaths[key] = owner
    }
    def hicBySample = [:].withDefault { [] }
    libraries.rows.each { row ->
        safeId(row.sample_id, 'sample_id')
        safeId(row.library_id, 'library_id')
        safeId(row.readset_id, 'readset_id')
        if (!seenSamples.contains(row.sample_id)) error "Hi-C table references unknown sample '${row.sample_id}'"
        if (!seenReadsets.add([row.sample_id, row.readset_id]))
            error "Repeated Hi-C readset_id '${row.readset_id}' in sample '${row.sample_id}'"
        if (!row.hic_r1 || !row.hic_r2) error "Hi-C read set '${row.readset_id}' needs both mates"
        def r1 = resolve(row.hic_r1, libraries.sheet)
        def r2 = resolve(row.hic_r2, libraries.sheet)
        checkPath(r1, "${row.sample_id}/${row.readset_id}/R1")
        checkPath(r2, "${row.sample_id}/${row.readset_id}/R2")
        hicBySample[row.sample_id] << [library_id: row.library_id, readset_id: row.readset_id, r1: r1, r2: r2]
    }
    def samples = []
    def status = []
    def readColumns = ['hifi_bam', 'hic_r1', 'hic_r2', 'sr_r1', 'sr_r2', 'tellseq_r1', 'tellseq_r2']
    table.rows.each { row ->
        def sid = row.sample_id
        def reads = [:]
        readColumns.each { key ->
            reads[key] = resolve(row[key], table.sheet)
            checkPath(reads[key], "${sid}/${key}")
        }
        def sets = hicBySample[sid].sort { a, b -> a.readset_id <=> b.readset_id }
        if (sets && (reads.hic_r1 || reads.hic_r2))
            error "Sample '${sid}' has both inline Hi-C and --hic_readsets; leave its inline cells empty"
        try {
            [['hic_r1','hic_r2'], ['sr_r1','sr_r2'], ['tellseq_r1','tellseq_r2']].each { pair ->
                if ((reads[pair[0]] != null) != (reads[pair[1]] != null))
                    throw new IllegalArgumentException("${pair.join('/')} requires both mates")
            }
            if (!sets && reads.hic_r1)
                sets = [[library_id: 'library1', readset_id: 'readset1', r1: reads.hic_r1, r2: reads.hic_r2]]
            sets.each { rs ->
                if (!seenOutputIds.add("${sid}__${rs.readset_id}".toString()))
                    error "Hi-C output name collision: ${sid}__${rs.readset_id}; choose distinct identifiers"
            }
            def paths = reads.values().findAll { it != null } + sets.collectMany { [it.r1, it.r2] }
            def missing = paths.findAll { !it.exists() }.unique()
            if (missing) throw new IllegalArgumentException("Missing reads: ${missing.join(', ')}")
            if (!reads.hifi_bam && !reads.sr_r1) {
                def reason = 'No supported assembly reads (HiFi or shotgun); Hi-C-only row retained as an assembly skip'
                log.info "Skipping '${sid}': ${reason}"
                status << [sample: sid, status: 'skipped', reason: reason, expected: 0]
            } else {
                def meta = buildMeta(row + [hifi: reads.hifi_bam != null, hic: !sets.isEmpty(),
                    shortread: reads.sr_r1 != null, tellseq: reads.tellseq_r1 != null])
                meta = meta + [hic_readsets: sets.collect { [library_id: it.library_id, readset_id: it.readset_id] }]
                reads.hic_sets = sets
                samples << tuple(meta, reads)
                status << [sample: sid, status: 'accepted', reason: "${meta.assembler}; ${sets.size()} Hi-C read sets", expected: meta.n_hap]
                log.info "Accepted '${sid}': assembler=${meta.assembler}, outputs=${meta.n_hap}, Hi-C libraries=${sets*.library_id.unique().size()}, read sets=${sets.size()}"
            }
        } catch (IllegalArgumentException e) {
            log.warn "Skipping invalid sample '${sid}': ${e.message}"
            status << [sample: sid, status: 'invalid', reason: e.message, expected: 0]
        }
    }
    if (!samples) error 'No assembly-capable valid samples remain; see validation messages'
    if (params.run_pangenome && status.any { it.status == 'invalid' })
        error 'Invalid assembly input(s): correct or explicitly remove them before enabling the pangenome'
    def metas = samples.collect { it[0] }
    def capabilities = [
        harmonize: metas.findAll { it.assembler == 'hifiasm' }.groupBy { it.taxid?.toString() }
            .any { taxid, group -> group.sum { (it.n_hap ?: 1) as int } >= 2 },
        hifi: metas.any { it.hifi },
        hic: metas.any { it.hic },
        shortread: metas.any { it.shortread },
        hifiasm: metas.any { it.assembler == 'hifiasm' },
        spades: metas.any { it.assembler == 'spades' },
        purge: metas.any { it.assembler == 'hifiasm' && it.dedup == 'purge_dups' },
        scaffold: metas.any { it.assembler == 'hifiasm' && it.hic },
        shortread_organelle: metas.any { !it.hifi }
    ]
    return [samples: samples, status: status, capabilities: capabilities]
}
