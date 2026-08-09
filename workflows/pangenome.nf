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
include { PANGENOME_STATS     } from '../modules/pangenome_stats.nf'
include { PANGENOME_REF_FASTA } from '../modules/pangenome_ref_fasta.nf'
include { PANGENOME_VARIANTS  } from '../modules/pangenome_variants.nf'
include { PANGENOME_MANIFEST  } from '../modules/pangenome_manifest.nf'
include { PANGENOME_GROWTH    } from '../modules/pangenome_growth.nf'
include { PANGENOME_PLOTS     } from '../modules/pangenome_plots.nf'
include { PANGENOME_2D_VIZ    } from '../modules/pangenome_2d_viz.nf'
include { PANGENOME_QC        } from '../modules/pangenome_qc.nf'
include { PANGENOME_ODGI_STATS_MQC } from '../modules/pangenome_odgi_stats_mqc.nf'
include { MULTIQC_PANGENOME   } from '../modules/multiqc_pangenome.nf'

workflow PANGENOME {

    take:
    ch_finalized       // tuple(meta, fasta, fai)
    ch_reference_ids   // tuple(taxid, ref_id)
    ch_species         // tuple(taxid, species_name)  resolved from the taxid (RESOLVE_TAXONOMY)

    main:
    ch_versions  = Channel.empty()
    ch_gbz       = Channel.empty()
    ch_gfa       = Channel.empty()
    ch_og        = Channel.empty()
    ch_snarls    = Channel.empty()
    ch_hapl      = Channel.empty()
    ch_vcf       = Channel.empty()
    ch_chrom_og  = Channel.empty()
    ch_viz       = Channel.empty()
    ch_ref_fasta = Channel.empty()
    ch_variants  = Channel.empty()
    ch_sv_sizes  = Channel.empty()
    ch_manifest  = Channel.empty()
    ch_stats     = Channel.empty()
    ch_growth    = Channel.empty()
    ch_growth_hist = Channel.empty()
    ch_core_acc    = Channel.empty()
    ch_figures     = Channel.empty()
    ch_growth_fit  = Channel.empty()
    ch_viz2d       = Channel.empty()
    ch_qc          = Channel.empty()
    ch_mqc         = Channel.empty()

    if( params.run_pangenome ) {
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
            .join(ch_species)

        // gate + PanSN naming per species (flatMap: return [] to drop a species)
        ch_cactus_in = ch_metrics
            .flatMap { taxid, members, ref_id, species ->
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
                def refName = flat(ref_id)
                def named = kept.collect { mm ->
                    def id  = mm.meta.id
                    def nm  = (sampleOf(id) == refInd) ? flat(id) : "${sampleOf(id)}.${hapOf(id)}"
                    [name: nm, fa: mm.fa]
                }.sort { a, b ->                          // deterministic order: reference first,
                    def ra = (a.name == refName) ? 0 : 1  // then by name. groupTuple() emits members
                    def rb = (b.name == refName) ? 0 : 1  // in task-completion order, so without this
                    ra != rb ? ra <=> rb : a.name <=> b.name  // the val names/fastas hash unstably ->
                }                                         // CACTUS_PANGENOME cache-misses every run.
                // key for publishDir / --outName: prefer the taxid-resolved species name
                // (RESOLVE_TAXONOMY), then a sample-sheet meta.species, then the taxid
                def sp    = species ?: refm.meta.species ?: members[0].meta.species
                def label = (sp?.toString()?.trim())
                    ? sp.toString().trim().replaceAll(/[^A-Za-z0-9._-]+/, '_').replaceAll(/^_+|_+$/, '')
                    : taxid
                return [ tuple(label, refName, named*.name, named*.fa) ]
            }

        CACTUS_PANGENOME( ch_cactus_in )
        ch_versions = ch_versions.mix( CACTUS_PANGENOME.out.versions )

        // graph statistics (clip graph)
        PANGENOME_STATS( CACTUS_PANGENOME.out.gbz.join( CACTUS_PANGENOME.out.og ) )

        // reference-path FASTA (needs the PanSN reference name from the cactus input)
        ch_ref_name = ch_cactus_in.map { taxid, ref_name, names, fastas -> tuple(taxid, ref_name) }
        PANGENOME_REF_FASTA( ch_ref_name.join( CACTUS_PANGENOME.out.gbz ) )

        // variant catalog (SNP / indel / SV counts + SV size spectrum) from the raw VCF
        PANGENOME_VARIANTS( CACTUS_PANGENOME.out.raw_vcf )
        ch_versions = ch_versions.mix( PANGENOME_VARIANTS.out.versions )

        // openness / growth (panacus on the finished clip GFA; workstream E)
        if( params.pangenome_growth != false ) {
            PANGENOME_GROWTH( CACTUS_PANGENOME.out.gfa )
            ch_versions    = ch_versions.mix( PANGENOME_GROWTH.out.versions )
            ch_growth      = PANGENOME_GROWTH.out.histgrowth
            ch_growth_hist = PANGENOME_GROWTH.out.hist
            ch_core_acc    = PANGENOME_GROWTH.out.core_accessory

            // report figures: growth/core + Heaps + band (from the coverage histogram),
            // SV size spectrum + variant-class bar (from the catalog) — workstream D
            def plots_script = file("${projectDir}/r_scripts/pangenome_plots.R", checkIfExists: true)
            PANGENOME_PLOTS(
                PANGENOME_GROWTH.out.hist
                    .join( PANGENOME_VARIANTS.out.sv_sizes )
                    .join( PANGENOME_VARIANTS.out.summary ),
                plots_script
            )
            ch_figures    = PANGENOME_PLOTS.out.figures
            ch_growth_fit = PANGENOME_PLOTS.out.growth_fit
        }

        // 2D per-chromosome layout (odgi; workstream D) — scatter so the ~25-min-per-chrom
        // path-guided SGD layouts run in parallel rather than one long serial task
        if( params.pangenome_2d_viz != false ) {
            ch_chrom_scatter = CACTUS_PANGENOME.out.chrom_og
                .flatMap { taxid, ogs ->
                    ogs.findAll { !it.name.endsWith('.full.og') }.collect { og -> tuple(taxid, og) }
                }
            PANGENOME_2D_VIZ( ch_chrom_scatter )
            ch_versions = ch_versions.mix( PANGENOME_2D_VIZ.out.versions )
            ch_viz2d    = PANGENOME_2D_VIZ.out.png.groupTuple()
        }

        // graph-intrinsic quality diagnostics (workstream C, Tier 1) — all default on
        if( params.pangenome_qc != false ) {
            PANGENOME_QC(
                CACTUS_PANGENOME.out.gbz
                    .join( CACTUS_PANGENOME.out.og )
                    .join( CACTUS_PANGENOME.out.gaf )
            )
            ch_versions = ch_versions.mix( PANGENOME_QC.out.versions )
            ch_qc       = PANGENOME_QC.out.metrics
        }

        // MultiQC report: odgi graph stats (whole + per-chrom) + bcftools variant stats
        if( params.pangenome_multiqc != false ) {
            PANGENOME_ODGI_STATS_MQC(
                CACTUS_PANGENOME.out.og.join( CACTUS_PANGENOME.out.chrom_og )
            )
            ch_versions = ch_versions.mix( PANGENOME_ODGI_STATS_MQC.out.versions )

            ch_mqc_in = PANGENOME_ODGI_STATS_MQC.out.yaml
                .join( PANGENOME_VARIANTS.out.bcftools_stats )
                .map { taxid, yamls, bcf -> tuple(taxid, [yamls, bcf].flatten()) }
            MULTIQC_PANGENOME( ch_mqc_in )
            ch_versions = ch_versions.mix( MULTIQC_PANGENOME.out.versions )
            ch_mqc      = MULTIQC_PANGENOME.out.report
        }

        // downstream manifest (role -> file) over the graph products
        PANGENOME_MANIFEST(
            CACTUS_PANGENOME.out.gbz
                .join( CACTUS_PANGENOME.out.gfa )
                .join( CACTUS_PANGENOME.out.og )
                .join( CACTUS_PANGENOME.out.snarls )
                .join( CACTUS_PANGENOME.out.hapl )
                .join( PANGENOME_VARIANTS.out.vcf )
                .join( PANGENOME_VARIANTS.out.tbi )
                .join( PANGENOME_REF_FASTA.out.ref_fasta )
                .join( PANGENOME_REF_FASTA.out.ref_fai )
        )

        ch_gbz       = CACTUS_PANGENOME.out.gbz
        ch_gfa       = CACTUS_PANGENOME.out.gfa
        ch_og        = CACTUS_PANGENOME.out.og
        ch_snarls    = CACTUS_PANGENOME.out.snarls
        ch_hapl      = CACTUS_PANGENOME.out.hapl
        ch_vcf       = CACTUS_PANGENOME.out.vcf
        ch_chrom_og  = CACTUS_PANGENOME.out.chrom_og
        ch_viz       = CACTUS_PANGENOME.out.viz
        ch_ref_fasta = PANGENOME_REF_FASTA.out.ref_fasta
        ch_variants  = PANGENOME_VARIANTS.out.summary
        ch_sv_sizes  = PANGENOME_VARIANTS.out.sv_sizes
        ch_manifest  = PANGENOME_MANIFEST.out.manifest
        ch_stats     = PANGENOME_STATS.out.vg_stats
    }

    emit:
    gbz       = ch_gbz          // graph + haplotype index (downstream mapping)
    gfa       = ch_gfa          // interchange (panacus / Bandage / odgi)
    og        = ch_og           // odgi graph (stats / report)
    snarls    = ch_snarls       // bubbles (downstream vg call)
    hapl      = ch_hapl         // haplotype-sampling index (downstream giraffe)
    vcf       = ch_vcf          // decomposed variants vs reference
    chrom_og  = ch_chrom_og     // per-chromosome odgi graphs (report)
    viz       = ch_viz          // 1D odgi viz PNGs (report)
    ref_fasta = ch_ref_fasta    // reference-path FASTA (downstream surjection)
    variants  = ch_variants     // variant summary table (report)
    sv_sizes  = ch_sv_sizes     // SV size spectrum (report)
    manifest  = ch_manifest     // role -> file downstream manifest
    stats     = ch_stats        // vg/odgi stats (report)
    growth       = ch_growth        // panacus growth/core curves (report)
    growth_hist  = ch_growth_hist   // panacus coverage histogram (report)
    core_accessory = ch_core_acc    // core/accessory/private partition (report)
    figures      = ch_figures       // rendered report PNGs (growth/SV/coverage/variants)
    growth_fit   = ch_growth_fit    // machine-readable Heaps gamma / open-closed / sizes
    viz2d        = ch_viz2d         // per-chromosome 2D layout PNGs
    qc           = ch_qc            // graph-intrinsic QC metrics (report)
    multiqc      = ch_mqc           // pangenome MultiQC report (odgi + bcftools)
    versions  = ch_versions
}
