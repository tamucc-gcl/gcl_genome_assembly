/*
========================================================================================
    PANGENOME SUBWORKFLOW  (minigraph-cactus, per species)
========================================================================================
    Repo location: workflows/pangenome.nf

    Runs after FINALIZE, per species (meta.taxid). Gates the finalized assemblies to the
    graph-ready set, PanSN-names them, and builds one minigraph-cactus graph per species
    using the harmonization reference as the graph reference.

    Gate:
      - long-read only (short-read excluded; nuclear-only is already guaranteed upstream)
      - OPTIONAL contiguity gate, OFF BY DEFAULT. When pangenome_max_chrom_scaffold_mult
        is set (> 0), exclude assemblies whose chromosome-scale scaffold count
        (>= finalize_min_scaffold_bp) exceeds that multiple of the reference chromosome
        count. OFF because minigraph-cactus does not require chromosome-scale INPUT:
        only the reference defines the chromosome components, and every other sample's
        contigs are assigned to a component by minigraph alignment. Contigs that map
        nowhere confidently are dropped PER CONTIG as ambiguous by
        cactus-graphmap-split, which is the correct granularity -- so a sub-chromosome-
        scale but otherwise sound assembly is a real individual and belongs in the
        graph. (NB a hifiasm primary assembly is never a redundant second copy of
        an individual already in the graph: forkHaplotypeMeta() returns EITHER
        [<sample>_primary] OR [<sample>_hap1, <sample>_hap2], never both.)
      - COMPARABILITY gate: drop a collapsed (n_hap == 1) assembly when the same species
        also has phased (n_hap == 2) assemblies. A collapsed assembly is a PHASE MOSAIC,
        not a haplotype -- it switches parental haplotype along the genome and drops one
        allele at heterozygous sites -- so it is not a comparable unit for the graph, for
        panacus --groupby-haplotype (Heaps' gamma, core/accessory), or for the ordination.
        When EVERY assembly of a species is collapsed the species is haploid (or uniformly
        collapsed), no phased peer exists, and all are kept. Ploidy is inferred at the
        GROUP level because main.nf strips meta.ploidy onto a sample-keyed side-channel.
        Override with pangenome_allow_collapsed = true.
      - the reference is always kept
      - require >= pangenome_min_haplotypes kept assemblies

    The kept/dropped set is written to the log as "[PANGENOME] taxid <n>: ...", so
    an exclusion appears in .nextflow.log rather than happening silently.

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
include { PANGENOME_CLASSIFY  } from '../modules/pangenome_classify.nf'
include { PANGENOME_INVERSION_RESCUE } from '../modules/pangenome_inversion_rescue.nf'
include { PANGENOME_STEPINDEX } from '../modules/pangenome_stepindex.nf'
include { PANGENOME_UNTANGLE  } from '../modules/pangenome_untangle.nf'
include { PANGENOME_REARRANGE } from '../modules/pangenome_rearrange.nf'
include { PANGENOME_MANIFEST  } from '../modules/pangenome_manifest.nf'
include { PANGENOME_GROWTH    } from '../modules/pangenome_growth.nf'
include { PANGENOME_PLOTS     } from '../modules/pangenome_plots.nf'
include { PANGENOME_HAP_COVERAGE } from '../modules/pangenome_hap_coverage.nf'
include { PANGENOME_2D_VIZ    } from '../modules/pangenome_2d_viz.nf'
include { PANGENOME_QC        } from '../modules/pangenome_qc.nf'
include { PANGENOME_ODGI_STATS_MQC } from '../modules/pangenome_odgi_stats_mqc.nf'
include { MULTIQC_PANGENOME   } from '../modules/multiqc_pangenome.nf'
include { PANGENOME_REPORT    } from '../modules/pangenome_report.nf'
include { PANGENOME_PCA_NJ    } from '../modules/pangenome_pca_nj.nf'
include { PANGENOME_POPSTRUCT } from '../modules/pangenome_popstruct.nf'
include { PANGENOME_PROGRESSIVE      } from '../modules/pangenome_progressive.nf'
include { PANGENOME_PROGRESSIVE_PLOT } from '../modules/pangenome_progressive_plot.nf'

workflow PANGENOME {

    take:
    ch_finalized       // tuple(meta, fasta, fai)
    ch_reference_ids   // tuple(taxid, ref_id)
    ch_species         // tuple(taxid, species_name)  resolved from the taxid (RESOLVE_TAXONOMY)
    ch_harm_report     // tuple(taxid, harmonization_report.tsv) -- PANGENOME_REARRANGE joins
                       // it to tell real inversions from chimeric joins. May be empty.

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
    ch_candidates = Channel.empty()
    ch_untangle   = Channel.empty()
    ch_manifest  = Channel.empty()
    ch_stats     = Channel.empty()
    ch_growth    = Channel.empty()
    ch_growth_hist = Channel.empty()
    ch_core_acc    = Channel.empty()
    ch_figures     = Channel.empty()
    ch_growth_fit  = Channel.empty()
    ch_hap_cov     = Channel.empty()
    ch_hap_priv    = Channel.empty()
    ch_viz2d       = Channel.empty()
    ch_qc          = Channel.empty()
    ch_mqc         = Channel.empty()
    ch_report      = Channel.empty()
    ch_pca_png     = Channel.empty()
    ch_nj_png      = Channel.empty()
    ch_progressive = Channel.empty()
    ch_prog_png    = Channel.empty()

    if( params.run_pangenome ) {
        def min_scaf = (params.finalize_min_scaffold_bp ?: 1000000) as long
        // contiguity-gate multiplier. null / '' / 'false' / <= 0  =>  GATE OFF (default).
        // Kept as a param so a run can re-enable it, but no longer a default-on filter.
        // A non-numeric, non-disabling value throws rather than falling back to a
        // default: a silent wrong-value run is the failure mode to avoid here.
        def multRaw  = params.pangenome_max_chrom_scaffold_mult
        def multStr  = multRaw?.toString()?.trim()
        def multOff  = (!multStr || multStr.toLowerCase() in ['null', 'false', 'off', 'none'])
        if( !multOff && !multStr.isNumber() )
            throw new IllegalArgumentException(
                "pangenome_max_chrom_scaffold_mult must be a number > 0, or one of " +
                   "null/false/off/0 to disable the contiguity gate (got '${multRaw}')")
        def mult     = (multOff || (multStr as BigDecimal) <= 0) ? null : (multStr as BigDecimal)
        def min_hap  = (params.pangenome_min_haplotypes ?: 2)
        // comparability-gate toggle. Parsed strictly: a non-boolean value is a typo, and a
        // silent fall-back to the default is the failure mode to avoid.
        def acRaw    = params.pangenome_allow_collapsed
        def acStr    = acRaw?.toString()?.trim()?.toLowerCase()
        if( acStr && !(acStr in ['true', 'false']) )
            throw new IllegalArgumentException(
                "pangenome_allow_collapsed must be true or false (got '${acRaw}')")
        def allowCol = (acStr == 'true')

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

                // ---- comparability gate -------------------------------------------------
                // A collapsed (n_hap == 1) assembly is a PHASE MOSAIC, not a haplotype: it
                // switches parental haplotype along the genome and drops one allele at every
                // heterozygous site. Mixed into a graph of phased haplotypes it is not a
                // comparable unit -- it lands at an artificial intermediate position in the
                // PCoA/NJ (spuriously similar to everything, interpretable as nothing), and
                // it breaks the one-unit-per-haploid-genome assumption behind panacus
                // --groupby-haplotype, biasing Heaps' gamma in an uncontrolled direction.
                //
                // meta.ploidy is NOT available here: main.nf strips it onto a sample-keyed
                // side-channel so ploidy tweaks stay off the task hash. Ploidy is therefore
                // inferred at the GROUP level -- if any assembly of this species is phased,
                // a collapsed assembly of the same species is a collapsed DIPLOID. If every
                // assembly in the group is collapsed the species is haploid (or uniformly
                // collapsed): no phased peer exists, nothing is incomparable, all are kept.
                def nhap       = { m -> ((m.meta.n_hap ?: 2) as int) }
                def anyPhased  = members.any { nhap(it) > 1 }
                def comparable = { m -> allowCol || !anyPhased || nhap(m) > 1 }
                if( anyPhased && nhap(refm) == 1 )
                    log.warn("[PANGENOME] taxid ${taxid}: graph reference ${ref_id} is a " +
                             "COLLAPSED (n_hap=1) assembly while the cohort also has phased " +
                             "haplotypes. It is kept as the reference (the harmonization frame " +
                             "was built against it), but the reference path is a phase mosaic, " +
                             "so the variant catalog and every reference-relative coordinate " +
                             "sit on a chimeric frame.")

                // ---- kept set: comparability AND (optional) contiguity -------------------
                // mult == null => contiguity gate off (the default): only the REFERENCE has to
                // be chromosome-scale, since non-reference contigs are assigned to reference
                // chromosome components by minigraph alignment.
                def kept = members.findAll { m ->
                    m.meta.id == ref_id ||
                        (comparable(m) && (mult == null || m.cs <= mult * ref_chrom))
                }

                // make the gate decision visible, with the REASON per dropped assembly: a
                // dropped haplotype used to vanish silently into a smaller graph.
                def keptIds = kept.collect { it.meta.id }.sort()
                def why     = { m ->
                    def r = []
                    if( !comparable(m) )                                r << "collapsed(n_hap=${nhap(m)})"
                    if( mult != null && m.cs > mult * ref_chrom )        r << "scaffolds(cs=${m.cs})"
                    r ? r.join('+') : 'unknown'
                }
                def dropMsg = members.findAll { !(it.meta.id in keptIds) }
                                     .sort { it.meta.id }
                                     .collect { "${it.meta.id}[${why(it)}]" }
                def gateMsg = (mult == null) ? 'contiguity gate OFF'
                                            : "contiguity gate cs <= ${mult} * ${ref_chrom}"
                def compMsg = allowCol ? 'comparability gate OFF (allow_collapsed)'
                                       : (anyPhased ? 'comparability gate ON'
                                                    : 'comparability gate n/a (no phased peers)')
                log.info("[PANGENOME] taxid ${taxid}: ref=${ref_id} (${ref_chrom} chrN_p names); " +
                         "kept ${keptIds.size()}/${members.size()} [${keptIds.join(' ')}]; " +
                         "${gateMsg}; ${compMsg}" +
                         (dropMsg ? "; DROPPED ${dropMsg.join(' ')}" : ''))
                if( kept.size() < min_hap ) {
                    log.warn("[PANGENOME] taxid ${taxid}: only ${kept.size()} assemblies survive " +
                             "the gate (pangenome_min_haplotypes=${min_hap}) -- NO GRAPH built " +
                             "for this species")
                    return []
                }
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

        // decomposition only: parent tier (LV==0) + fine tier (vcfbub -> vcfwave -> norm).
        // Classification moved to PANGENOME_CLASSIFY, which reads the graph's own allele
        // traversals instead of REF/ALT string lengths.
        PANGENOME_VARIANTS( CACTUS_PANGENOME.out.raw_vcf )
        ch_versions = ch_versions.mix( PANGENOME_VARIANTS.out.versions )

        // parents_vcf feeds BOTH classification and inversion rescue. Referencing it twice
        // starves one of them: on the first run PANGENOME_INVERSION_RESCUE got the item and
        // PANGENOME_CLASSIFY submitted only its `fine` task, so the PARENT tier -- the only
        // view where AT is interpretable, and therefore the only place topology is valid --
        // silently never ran. No error, no warning, just a missing analysis.
        //
        // multiMap forks a single read into named outputs so both consumers are fed. It sits
        // OUTSIDE both `if` blocks on purpose: defined inside the classify block it would be
        // undefined whenever pangenome_classify = false, breaking rescue with a missing
        // property error for a reason that has nothing to do with rescue.
        PANGENOME_VARIANTS.out.parents_vcf.view { "PARENTS_VCF: $it" }

        PANGENOME_VARIANTS.out.parents_vcf
            .multiMap { taxid, vcf ->
                classify: tuple(taxid, 'clip', 'parent', vcf)
                rescue:   tuple(taxid, 'clip', vcf)
            }
            .set { ch_parents }

        // ---- topological classification + allele frequencies ---------------------------
        // The PARENT tier is the only view where AT is interpretable: vcfwave and
        // `bcftools norm -m` rewrite REF/ALT while AT is inherited from the parent record,
        // so allele i stops corresponding to traversal i+1. The module refuses to compute
        // topology on a decomposed VCF rather than trusting this wiring to be right.
        if( params.pangenome_classify != false ) {
            def classify_script = file("${projectDir}/py_scripts/classify_variants.py",
                                       checkIfExists: true)

            ch_classify_in = ch_parents.classify
                .view { "PARENT_BRANCH: $it" }
                .mix( PANGENOME_VARIANTS.out.vcf
                        .map { taxid, vcf -> tuple(taxid, 'clip', 'fine', vcf) } )
                .view { "CLASSIFY_IN: $it" }

            // Pair each VCF with its OWN graph, by taxid, and hand the process two channels
            // of EQUAL cardinality.
            //
            // Passing CACTUS_PANGENOME.out.gfa directly was the bug: it is a queue channel
            // with one item, ch_classify_in has two (parent and fine), and Nextflow consumes
            // input channels in LOCKSTEP -- one item from each per task. Two against one
            // makes one task, and the parent tuple was discarded with no error and nothing
            // in the log. The parent tier is the only view where AT is interpretable, so
            // topological classification never ran at all while the fine tier published a
            // complete set of files and made the absence look like a downstream problem.
            //
            // .first() would also make both tasks run, but takes the first item regardless
            // of taxid -- correct for one species, silently wrong for two, and wrong in a
            // way that produces plausible numbers rather than an error.
            ch_classify_in
                .combine( CACTUS_PANGENOME.out.gfa, by: 0 )
                .multiMap { taxid, flavor, tier, vcf, gfa ->
                    vcfs: tuple(taxid, flavor, tier, vcf)
                    gfas: tuple(taxid, gfa)
                }
                .set { ch_cls }

            PANGENOME_CLASSIFY( ch_cls.vcfs, ch_cls.gfas, classify_script )
            ch_versions = ch_versions.mix( PANGENOME_CLASSIFY.out.versions )

            // canonical arm for the existing per-taxid consumers: clip + parent, mapped back
            // to tuple(taxid, file). The fine tier publishes alongside without entering the
            // report; making the report flavour-aware is the report-matrix work.
            ch_variants = PANGENOME_CLASSIFY.out.summary
                .filter { taxid, flavor, tier, f -> flavor == 'clip' && tier == 'parent' }
                .map    { taxid, flavor, tier, f -> tuple(taxid, f) }
            ch_sv_sizes = PANGENOME_CLASSIFY.out.spectrum
                .filter { taxid, flavor, tier, f -> flavor == 'clip' && tier == 'parent' }
                .map    { taxid, flavor, tier, f -> tuple(taxid, f) }
        }

        // ---- alignment-rescued inversions ----------------------------------------------
        // Recovers inversions whose alleles are DISJOINT node sets and therefore look like
        // allele replacement to topology. Separate process because it is minimap2 over
        // ~310k allele pairs while classification is minutes; its thresholds will be tuned.
        if( params.pangenome_inv_rescue != false ) {
            def rescue_script = file("${projectDir}/py_scripts/rescue_inversions.py",
                                     checkIfExists: true)
            // ch_parents.rescue, not a second read of parents_vcf -- see the multiMap above
            PANGENOME_INVERSION_RESCUE( ch_parents.rescue, rescue_script )
            ch_versions = ch_versions.mix( PANGENOME_INVERSION_RESCUE.out.versions )
        }

        // ---- rearrangement: odgi untangle ----------------------------------------------
        // The PRIMARY inversion / duplication instrument. chr10 alone yields ~21.5 Mb of
        // inverted sequence against ~11 Mb genome-wide from both bubble-based detectors
        // combined, and self.cov > 1 is the only duplication signal available at all.
        //
        // FULL graph first: clipping cuts paths into subpaths (556 vs 394 on chr10) and a
        // rearrangement straddling a boundary is lost to path projection. Per chromosome
        // because the whole-graph clip .og is 118 GB.
        if( params.pangenome_untangle != false ) {
            def rearr_script = file("${projectDir}/py_scripts/rearrange_from_untangle.py",
                                    checkIfExists: true)
            def unt_flavors = (params.pangenome_untangle_flavors ?: 'full,clip')
                                  .toString().split(',').collect { it.trim() }

            ch_unt_scatter = Channel.empty()
            if( unt_flavors.contains('full') ) {
                ch_unt_scatter = ch_unt_scatter.mix(
                    CACTUS_PANGENOME.out.chrom_og_full.flatMap { taxid, ogs ->
                        (ogs instanceof List ? ogs : [ogs]).collect { og -> tuple(taxid, 'full', og) }
                    } )
            }
            if( unt_flavors.contains('clip') ) {
                ch_unt_scatter = ch_unt_scatter.mix(
                    CACTUS_PANGENOME.out.chrom_og.flatMap { taxid, ogs ->
                        (ogs instanceof List ? ogs : [ogs])
                            .findAll { !it.name.endsWith('.full.og') }
                            .collect { og -> tuple(taxid, 'clip', og) }
                    } )
            }

            // stepindex is a SEPARATE process for two reasons: odgi v0.9.2 in the cactus
            // image crashes building the index inside untangle, and the index is invariant
            // to every untangle parameter, so it survives retuning -e / -j / -n.
            PANGENOME_STEPINDEX( ch_unt_scatter )
            PANGENOME_UNTANGLE( PANGENOME_STEPINDEX.out.stpidx )
            ch_versions = ch_versions.mix( PANGENOME_STEPINDEX.out.versions )
            ch_versions = ch_versions.mix( PANGENOME_UNTANGLE.out.versions )
            ch_untangle = PANGENOME_UNTANGLE.out.tsv

            // one NO_FILE placeholder per taxid so an absent harmonization report cannot
            // silently stop REARRANGE from ever running; the real report wins where present
            ch_harm_safe = ch_untangle
                .map { label, flavor, f -> label }
                .unique()
                .combine( ch_harm_report.map { taxid, rpt -> rpt }.ifEmpty( file('NO_FILE') ) )
                .map { label, rpt -> tuple(label, rpt) }
                .view { "HARM_SAFE: $it" }

            // groupTuple WITHOUT size: on purpose. The count per flavour is known (one per
            // chromosome), but PANGENOME_UNTANGLE carries errorStrategy 'ignore', so a
            // single failed chromosome would leave a fixed-size groupTuple waiting forever.
            // The barrier is the correct trade: collect whatever completed.
            //
            // combine by taxid, NOT a separate channel: the grouped channel emits once per
            // flavour while the report channel emits once per taxid, so passing them
            // separately would run this process once and silently drop a flavour.
            ch_rearr_in = ch_untangle.groupTuple( by: [0, 1] )
                .combine( ch_harm_safe, by: 0 )
                .view { "REARR_IN: $it" }

            PANGENOME_REARRANGE(
                ch_rearr_in
                    .map { taxid, flavor, tsvs, harm -> tuple(taxid, flavor, tsvs, harm) }
                    .view { "REARR_PASSED: $it" },
                rearr_script
            )
            ch_versions   = ch_versions.mix( PANGENOME_REARRANGE.out.versions )
            ch_candidates = PANGENOME_REARRANGE.out.candidates
        }

        // openness / growth (panacus on the finished clip GFA; workstream E)
        if( params.pangenome_growth != false ) {
            PANGENOME_GROWTH( CACTUS_PANGENOME.out.gfa )
            ch_versions    = ch_versions.mix( PANGENOME_GROWTH.out.versions )
            ch_growth      = PANGENOME_GROWTH.out.histgrowth
            ch_growth_hist = PANGENOME_GROWTH.out.hist
            ch_core_acc    = PANGENOME_GROWTH.out.core_accessory

            // WHO owns the private sequence. panacus's coverage histogram is a marginal
            // over haplotypes, so the private bar can never be attributed from it after the
            // fact; this reads the clip GFA directly and resolves the owning haplotype per
            // node. Same PanSN grouping as panacus, so the emitted hap_coverage_check.tsv
            // should reproduce the hist exactly. Gated with the rest of the growth analysis
            // -- no separate toggle.
            def hapcov_script = file("${projectDir}/py_scripts/gfa_hap_coverage.py", checkIfExists: true)
            PANGENOME_HAP_COVERAGE( CACTUS_PANGENOME.out.gfa, hapcov_script )
            ch_versions    = ch_versions.mix( PANGENOME_HAP_COVERAGE.out.versions )
            ch_hap_cov     = PANGENOME_HAP_COVERAGE.out.matrix
            ch_hap_priv    = PANGENOME_HAP_COVERAGE.out.hap_private

            // report figures: growth/core + Heaps + band (from the coverage histogram),
            // SV size spectrum + variant-class bar (from the catalog) — workstream D
            def plots_script = file("${projectDir}/r_scripts/pangenome_plots.R", checkIfExists: true)
            PANGENOME_PLOTS(
                PANGENOME_GROWTH.out.hist
                    .join( ch_sv_sizes )
                    .join( ch_variants )
                    .join( ch_hap_priv, remainder: true )
                    // remainder: true is right for hap_priv, which is genuinely optional --
                    // but it also emits UNMATCHED RIGHT-HAND entries when the left side is
                    // empty, as [taxid, null, hap_private] with a single null placeholder
                    // because Nextflow never saw the left channel's arity. Destructuring that
                    // into five parameters aborts the whole pipeline with a Groovy arity
                    // error that names a closure and says nothing about the cause. An absent
                    // variant catalog must mean "no plots", not "kill the run".
                    .filter { def l = it as List
                              l.size() >= 4 && l[1] != null && l[2] != null && l[3] != null }
                    .map { taxid, hist, sv, vs, hp ->
                        tuple(taxid, hist, sv, vs, hp ?: file('NO_HAP_PRIVATE')) },
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

        // graph-derived population structure: per-haplotype PCoA + NJ tree from odgi
        // similarity (grouped by PanSN haplotype -> 4 units incl. reference at the n=4 test,
        // so it plots for real). Workstream D.
        if( params.pangenome_popstruct != false ) {
            PANGENOME_PCA_NJ( CACTUS_PANGENOME.out.og )
            ch_versions = ch_versions.mix( PANGENOME_PCA_NJ.out.versions )

            def popstruct_script = file("${projectDir}/r_scripts/pangenome_popstruct.R", checkIfExists: true)
            PANGENOME_POPSTRUCT( PANGENOME_PCA_NJ.out.similarity, popstruct_script )
            ch_pca_png = PANGENOME_POPSTRUCT.out.pca_png
            ch_nj_png  = PANGENOME_POPSTRUCT.out.figures
        }

        // progressive (incremental-construction) growth via minigraph — opt-in (workstream H)
        if( params.pangenome_progressive ) {
            ch_prog_in = ch_cactus_in.map { taxid, ref, names, fastas -> tuple(taxid, names.join(' '), fastas) }
            PANGENOME_PROGRESSIVE( ch_prog_in )
            ch_versions = ch_versions.mix( PANGENOME_PROGRESSIVE.out.versions )

            def prog_script = file("${projectDir}/r_scripts/pangenome_progressive.R", checkIfExists: true)
            PANGENOME_PROGRESSIVE_PLOT( PANGENOME_PROGRESSIVE.out.growth, prog_script )
            ch_progressive = PANGENOME_PROGRESSIVE.out.growth
            ch_prog_png    = PANGENOME_PROGRESSIVE_PLOT.out.png
        }

        // full-graph product NAMES for the manifest. These are already declared and published
        // via CACTUS_PANGENOME's `all` catch-all; picking them by name here means that process
        // needs no new output declarations (its cache is expensive to lose). A product cactus
        // did not emit yields '' and the manifest omits its row rather than erroring.
        //
        // Matched on EXACT basename, not suffix: <label>.chroms/ contains per-chromosome files
        // named chrN_M.full.og, so endsWith('.full.og') would happily return a chromosome graph
        // instead of the whole-graph one.
        ch_full_names = CACTUS_PANGENOME.out.all
            .map { taxid, files ->
                def fl   = (files instanceof List) ? files : [files]
                def pick = { String sfx ->
                    def want = "${taxid}${sfx}"
                    def f = fl.find { it.name == want }
                    f ? f.name : ''
                }
                tuple(taxid, pick('.full.gbz'), pick('.full.gfa.gz'),
                             pick('.full.og'),  pick('.full.snarls'))
            }

        // downstream manifest (role -> file, per graph) over the graph products
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
                .join( ch_full_names )
        )

        // pangenome report fragment + stats JSON (workstream F) — always on when built.
        // Base on the always-present variant summary; feed NO_* sentinels for any sub-analysis
        // that was disabled (the R script skips sentinels).
        if( params.pangenome_report != false ) {
            def report_script = file("${projectDir}/r_scripts/pangenome_report.R", checkIfExists: true)
            ch_report_in = ch_variants
                .join( ch_qc,                          remainder: true )
                .join( ch_growth_fit,                  remainder: true )
                .join( PANGENOME_STATS.out.odgi_stats, remainder: true )
                .join( ch_pca_png,                     remainder: true )
                .join( ch_prog_png,                    remainder: true )
                .join( PANGENOME_MANIFEST.out.manifest, remainder: true )
                .join( ch_hap_priv,                     remainder: true )
                .map { taxid, vs, qc, gf, gs, pca, prog, mf, hp ->
                    tuple(taxid,
                          qc ?: file('NO_QC'),
                          gf ?: file('NO_GROWTH'),
                          vs ?: file('NO_VARIANTS'),
                          gs ?: file('NO_GRAPH'),
                          pca ?: file('NO_POPSTRUCT'),
                          prog ?: file('NO_PROGRESSIVE'),
                          mf ?: file('NO_MANIFEST'),
                          hp ?: file('NO_HAP_PRIVATE')) }
            PANGENOME_REPORT( ch_report_in, report_script )
            ch_versions = ch_versions.mix( PANGENOME_REPORT.out.versions )
            ch_report   = PANGENOME_REPORT.out.report
        }

        ch_gbz       = CACTUS_PANGENOME.out.gbz
        ch_gfa       = CACTUS_PANGENOME.out.gfa
        ch_og        = CACTUS_PANGENOME.out.og
        ch_snarls    = CACTUS_PANGENOME.out.snarls
        ch_hapl      = CACTUS_PANGENOME.out.hapl
        ch_vcf       = CACTUS_PANGENOME.out.vcf
        ch_chrom_og  = CACTUS_PANGENOME.out.chrom_og
        ch_viz       = CACTUS_PANGENOME.out.viz
        ch_ref_fasta = PANGENOME_REF_FASTA.out.ref_fasta
        // ch_variants / ch_sv_sizes are set by the PANGENOME_CLASSIFY block above; the
        // length-based awk that used to fill them here has been removed.
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
    sv_sizes  = ch_sv_sizes     // SV size spectrum, coarse-binned (report)
    candidates = ch_candidates  // rearrangement candidate loci, artifact-flagged
    untangle  = ch_untangle     // per-chromosome path projections
    manifest  = ch_manifest     // role -> file downstream manifest
    stats     = ch_stats        // vg/odgi stats (report)
    growth       = ch_growth        // panacus growth/core curves (report)
    growth_hist  = ch_growth_hist   // panacus coverage histogram (report)
    core_accessory = ch_core_acc    // core/accessory/private partition (report)
    figures      = ch_figures       // rendered report PNGs (growth/SV/coverage/variants)
    growth_fit   = ch_growth_fit    // machine-readable Heaps gamma / open-closed / sizes
    hap_coverage = ch_hap_cov       // joint haplotype x coverage-level bp matrix
    hap_private  = ch_hap_priv      // per-haplotype private-sequence breakdown
    viz2d        = ch_viz2d         // per-chromosome 2D layout PNGs
    qc           = ch_qc            // graph-intrinsic QC metrics (report)
    multiqc      = ch_mqc           // pangenome MultiQC report (odgi + bcftools)
    report       = ch_report        // pangenome report fragment (markdown) for the main report
    pca_png      = ch_pca_png        // haplotype PCoA (report indicator)
    nj_png       = ch_nj_png         // all popstruct figures (hap + individual PCoA/NJ)
    progressive  = ch_progressive    // progressive growth table (opt-in)
    progressive_png = ch_prog_png    // progressive growth curve (opt-in)
    versions  = ch_versions
}
