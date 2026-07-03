/*
========================================================================================
    META SCHEMA + BUILDERS — single source of truth for the `meta` map
========================================================================================
    Repo location: functions/meta.nf

    `meta` is the per-sample (later per-haplotype) carrier threaded through every process
    as the first element of tuple(meta, ...). Built ONCE here so the shape lives in one place.

    SAMPLE-LEVEL meta (emitted by parseSampleSheet):
      id           unique id (== sample at sample level)
      sample       sample identifier
      haplotype    null until assembly; then 'hap1' | 'hap2' | 'primary'
      n_hap        1 (haploid/collapsed) | 2 (diploid)  -> drives groupKey(sample, n_hap)
      hifi/hic/tellseq/shortread  booleans: which read types are present
      long_reads   derived: hifi (|| future ONT) -> gates teloclip / long-read gap-fill
      assembler    'hifiasm' | 'spades'
      dedup        'purge_dups' | 'redundans' | 'none'   (selectable, input-independent)
      mito_tool    'mitohifi' | 'mitofinder' | 'none'
      (scaffolding round count + scaffolder ordering are intentionally NOT in meta — see CACHING)

    PER-HAPLOTYPE meta (forkHaplotypeMeta — wired during Phase 1 module threading):
      id           "${sample}_hap1" | "${sample}_hap2" | "${sample}_primary"
      haplotype    'hap1' | 'hap2' | 'primary'
      (all other fields inherited unchanged)

    STRATEGY FIELD PRECEDENCE (assembler, ploidy, dedup, mito_tool):
      explicit per-row column value  >  global params.<field>  >  derived-from-data default

    BEHAVIOR PRESERVATION (Phase 1):
      - dedup default for hifiasm honors the existing params.run_purge_dups switch
        (so a HiFi+Hi-C run with run_purge_dups=false still resolves to dedup='none').
      - ploidy defaults to 'diploid'; mito_tool defaults to 'mitohifi' when HiFi is present.
      => the legacy 4-column sheet resolves to a meta that reproduces current behavior.

    CACHING:
      meta becomes part of each task's hash, so only stable per-sample identity/strategy fields
      live here. Volatile knobs (scaffolding round count, scaffolder ordering) are deliberately
      EXCLUDED — tuning e.g. params.hic_scaffold_rounds must not re-hash (and re-run) upstream
      assembly/QC tasks. Read params.hic_scaffold_rounds at point of use, and derive scaffolder
      ordering from meta.hic / meta.tellseq in routing (Phase 5).
========================================================================================
*/

// ---- Build the canonical sample-level meta map (validates; throws on bad config) ----
def buildMeta(Map a) {
    def sample  = a.sample_id?.toString()?.trim()
    def hasHifi = a.hifi      as boolean
    def hasHic  = a.hic       as boolean
    def hasTell = a.tellseq   as boolean
    def hasSr   = a.shortread as boolean

    if (!sample) throw new IllegalArgumentException("missing sample_id")
    if (!(hasHifi || hasSr))
        throw new IllegalArgumentException(
            "sample '${sample}': no contig-capable reads (need HiFi or short-read shotgun; " +
            "Hi-C and TellSeq are scaffolding-only)")

    // precedence helper: row value > global param > derived default
    def pick = { rowVal, String paramKey, Closure derived ->
        def rv = rowVal?.toString()?.trim()
        if (rv) return rv
        if (params.containsKey(paramKey) && params[paramKey] != null) return params[paramKey].toString()
        return derived()
    }

    def assembler = pick(a.assembler, 'assembler') { hasHifi ? 'hifiasm' : 'spades' }
    def ploidyRaw = pick(a.ploidy,    'ploidy')    { 'diploid' }
    // accept numeric (1|2) or string (haploid|diploid), case-insensitive
    def ploidy    = ['1':'haploid','haploid':'haploid','2':'diploid','diploid':'diploid'][ploidyRaw?.toString()?.trim()?.toLowerCase()]
    def dedup     = pick(a.dedup,     'dedup') {
        if (assembler == 'hifiasm')
            (params.containsKey('run_purge_dups') && params.run_purge_dups) ? 'purge_dups' : 'none'
        else
            'redundans'
    }
    def mito = pick(a.mito_tool, 'mito_tool') { hasHifi ? 'mitohifi' : 'mitofinder' }

    // enum validation
    if (!(assembler in ['hifiasm','spades']))
        throw new IllegalArgumentException("sample '${sample}': invalid assembler '${assembler}' (allowed: hifiasm, spades)")
    if (!(ploidy in ['haploid','diploid']))
        throw new IllegalArgumentException("sample '${sample}': invalid ploidy '${ploidyRaw}' (allowed: 1/haploid, 2/diploid)")
    if (!(dedup in ['purge_dups','redundans','none']))
        throw new IllegalArgumentException("sample '${sample}': invalid dedup '${dedup}' (allowed: purge_dups, redundans, none)")
    if (!(mito in ['mitohifi','mitofinder','none']))
        throw new IllegalArgumentException("sample '${sample}': invalid mito_tool '${mito}' (allowed: mitohifi, mitofinder, none)")

    // consistency: assembler must have its contig-source reads
    if (assembler == 'hifiasm' && !hasHifi)
        throw new IllegalArgumentException("sample '${sample}': assembler=hifiasm but no HiFi reads provided")
    if (assembler == 'spades' && !hasSr)
        throw new IllegalArgumentException("sample '${sample}': assembler=spades but no short-read shotgun provided")

    def n_hap = (ploidy == 'haploid') ? 1 : 2

    return [
        id:          sample,
        sample:      sample,
        haplotype:   null,
        n_hap:       n_hap,
        hifi:        hasHifi,
        hic:         hasHic,
        tellseq:     hasTell,
        shortread:   hasSr,
        long_reads:  hasHifi,          // (|| ONT later)
        assembler:   assembler,
        dedup:       dedup,
        mito_tool:   mito
    ]
}

// ---- Fork sample-level meta into per-haplotype metas (wired during Phase 1 threading) ----
//      Use after the assembler emits per-haplotype FASTAs:
//        ch_asm.flatMap { meta, fastas -> [forkHaplotypeMeta(meta), fastas].transpose()... }
//      (exact wiring handled when threading HIFIASM / the assembler selector)
def forkHaplotypeMeta(Map meta) {
    if (meta.n_hap == 1)
        return [ meta + [ id: "${meta.sample}_primary", haplotype: 'primary' ] ]
    return [
        meta + [ id: "${meta.sample}_hap1", haplotype: 'hap1' ],
        meta + [ id: "${meta.sample}_hap2", haplotype: 'hap2' ]
    ]
}
