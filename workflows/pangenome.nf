/*
========================================================================================
    PANGENOME SUBWORKFLOW  (minigraph-cactus, per species)
========================================================================================
    Repo location: workflows/pangenome.nf

    Runs after FINALIZE, per species (meta.taxid). Gates the finalized assemblies to the
    graph-ready set, PanSN-names them, and builds one minigraph-cactus graph per species
    using the harmonization reference as the graph reference.

    Gate (all tunable):
      - long-read only (short-read excluded; nuclear-only is already guaranteed upstream)
      - exclude assemblies whose chromosome-scale scaffold count (>= finalize_min_scaffold_bp)
        exceeds pangenome_max_chrom_scaffold_mult x the reference chromosome count
        (drops fragmented / primary-contig intermediates like *_hifi_primary)
      - the reference is always kept
      - require >= pangenome_min_haplotypes kept assemblies

    PanSN naming (cactus seqfile convention SAMPLE.HAPLOTYPE):
      - sample  = meta.id with a trailing _hap<N> stripped ('.' -> '_' to protect the
        haplotype separator)
      - the reference is written WITHOUT a haplotype suffix (cactus reference convention
        -> becomes haplotype 0 in the graph)
      - every other haplotype is <sample>.<N>

    Take : ch_finalized     = tuple(meta, fasta, fai)   (FINALIZE assembly joined with fai)
           ch_reference_ids = tuple(taxid, ref_id)      (from HARMONIZE_SCAFFOLDS)
    Emit : gbz / vcf / stats / versions
========================================================================================
*/

include { CACTUS_PANGENOME } from '../modules/cactus_pangenome.nf'
include { PANGENOME_STATS  } from '../modules/pangenome_stats.nf'

workflow PANGENOME {

    take:
    ch_finalized       // tuple(meta, fasta, fai)
    ch_reference_ids   // tuple(taxid, ref_id)

    main:
    ch_versions = Channel.empty()

    if( !params.run_pangenome ) {
        ch_gbz   = Channel.empty()
        ch_vcf   = Channel.empty()
        ch_stats = Channel.empty()
    }
    else {
        def min_scaf = (params.finalize_min_scaffold_bp ?: 1000000) as long
        def mult     = (params.pangenome_max_chrom_scaffold_mult ?: 3)
        def min_hap  = (params.pangenome_min_haplotypes ?: 2)

        // per-assembly gate metrics from the .fai (long-read only)
        ch_metrics = ch_finalized
            .filter { meta, fa, fai -> !meta.shortread }
            .map { meta, fa, fai ->
                def rows = fai.text.readLines().findAll { it?.trim() }.collect { it.split('\t') }
                def chrom_scale = rows.count { (it[1] as long) >= min_scaf }
                def chr_named   = rows.count { it[0] ==~ /^chr[0-9]+_[0-9]+$/ }
                tuple(meta.taxid.toString(), [meta: meta, fa: fa, cs: chrom_scale, cn: chr_named])
            }
            .groupTuple()
            .join(ch_reference_ids)

        // gate + PanSN naming per species (flatMap: return [] to drop a species)
        ch_cactus_in = ch_metrics
            .flatMap { taxid, members, ref_id ->
                def refm = members.find { it.meta.id == ref_id }
                if( refm == null ) return []
                def ref_chrom = (refm.cn ?: 1) as int
                def kept = members.findAll { it.meta.id == ref_id || it.cs <= mult * ref_chrom }
                if( kept.size() < min_hap ) return []
                def sampleOf = { String id -> id.replaceFirst(/_hap[0-9]+$/, '').replaceAll(/\./, '_') }
                def hapOf    = { String id -> def m = (id =~ /_hap([0-9]+)$/); m ? m[0][1] : '0' }
                def flat     = { String id -> id.replaceAll(/\./, '_') }   // full id, dots -> _
                // cactus forbids a sample name from sharing a prefix with the reference, so
                // the reference cannot be grouped with its own sibling haplotype. Assemblies
                // from the reference's individual therefore get distinct FLAT names (full id);
                // every other individual's haplotypes group as <individual>.<hap> (PanSN).
                def refInd = sampleOf(ref_id)
                def named = kept.collect { mm ->
                    def id  = mm.meta.id
                    def nm  = (sampleOf(id) == refInd) ? flat(id) : "${sampleOf(id)}.${hapOf(id)}"
                    [name: nm, fa: mm.fa]
                }
                // key used for publishDir / --outName: sanitized species name from
                // meta.species, falling back to the taxid when species is unset
                def sp    = refm.meta.species ?: members[0].meta.species
                def label = (sp?.toString()?.trim())
                    ? sp.toString().trim().replaceAll(/[^A-Za-z0-9._-]+/, '_').replaceAll(/^_+|_+$/, '')
                    : taxid
                return [ tuple(label, flat(ref_id), named*.name, named*.fa) ]
            }

        CACTUS_PANGENOME( ch_cactus_in )
        ch_versions = ch_versions.mix( CACTUS_PANGENOME.out.versions )

        PANGENOME_STATS( CACTUS_PANGENOME.out.gbz.join( CACTUS_PANGENOME.out.og ) )

        ch_gbz   = CACTUS_PANGENOME.out.gbz
        ch_vcf   = CACTUS_PANGENOME.out.vcf
        ch_stats = PANGENOME_STATS.out.vg_stats
    }

    emit:
    gbz      = ch_gbz
    vcf      = ch_vcf
    stats    = ch_stats
    versions = ch_versions
}
