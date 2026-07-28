/*
========================================================================================
    TAXONOMY HELPERS — resolved lineage -> BUSCO odb10 dataset + kingdom flag
========================================================================================
    Pure functions. Match on the FULL lineage string (robust to NCBI's inconsistent rank
    assignment, e.g. monocots filed under class Magnoliopsida). These feed a per-sample
    side-channel (not meta), so edits re-run only BUSCO/mito, never assembly. Extend freely.
========================================================================================
*/

def organismName(Map r) {
    def s = r.species?.toString()?.trim()
    def g = r.genus?.toString()?.trim()
    if (s && s != 'NA') return s
    if (g && g != 'NA') return g
    return null
}

def kingdomFlag(Map r) {
    def K   = (r.kingdom ?: '').toString().toLowerCase()
    def lin = (r.lineage ?: '').toString().toLowerCase()
    if (K.contains('viridiplantae') || lin.contains('viridiplantae')) return 'plant'
    if (K.contains('metazoa')       || lin.contains('metazoa'))       return 'animal'
    if (K.contains('fungi')         || lin.contains('fungi'))         return 'fungi'
    return 'other'
}

// Curated lineage -> BUSCO v5 odb10. Most-specific match wins; eukaryota fallback.
def buscoLineageFor(Map r) {
    def lin = (r.lineage ?: '').toString().toLowerCase()

    // Plants
    if (lin.contains('viridiplantae')) {
        if (lin.contains('liliopsida'))      return 'liliopsida_odb10'    // monocots
        if (lin.contains('eudicotyledon'))   return 'eudicots_odb10'
        return 'embryophyta_odb10'
    }
    // Fungi
    if (lin.contains('fungi')) return 'fungi_odb10'
    // Animals
    if (lin.contains('metazoa')) {
        if (lin.contains('chordata')) {
            if (lin.contains('actinopteri') || lin.contains('actinopterygii')) return 'actinopterygii_odb10'
            if (lin.contains('aves'))         return 'aves_odb10'
            if (lin.contains('mammalia'))     return 'mammalia_odb10'
            if (lin.contains('amphibia'))     return 'tetrapoda_odb10'
            if (lin.contains('lepidosauria') || lin.contains('testudines') ||
                lin.contains('crocodylia')   || lin.contains('archelosauria') ||
                lin.contains('archosauria')  || lin.contains('sauropsida'))    return 'sauropsida_odb10'
            return 'vertebrata_odb10'
        }
        if (lin.contains('arthropoda')) {
            if (lin.contains('hymenoptera'))  return 'hymenoptera_odb10'
            if (lin.contains('diptera'))      return 'diptera_odb10'
            if (lin.contains('lepidoptera'))  return 'lepidoptera_odb10'
            if (lin.contains('hemiptera'))    return 'hemiptera_odb10'
            if (lin.contains('endopterygota') || lin.contains('holometabola')) return 'endopterygota_odb10'  // e.g. Coleoptera
            if (lin.contains('insecta'))      return 'insecta_odb10'
            if (lin.contains('arachnida'))    return 'arachnida_odb10'
            return 'arthropoda_odb10'
        }
        if (lin.contains('nematoda'))         return 'nematoda_odb10'
        if (lin.contains('mollusca'))         return 'mollusca_odb10'
        return 'metazoa_odb10'
    }
    return 'eukaryota_odb10'
}

// NCBI mitochondrial genetic code (MitoHiFi -o). Clade-specific; verify/extend for your taxa.
def geneticCodeFor(Map r) {
    def lin = (r.lineage ?: '').toString().toLowerCase()
    if (lin.contains('viridiplantae'))   return 1    // plant mito: standard code
    if (lin.contains('fungi'))           return 4    // mold/protozoan mito (yeast = 3)
    if (lin.contains('metazoa')) {
        if (lin.contains('vertebrata') || lin.contains('craniata'))            return 2   // vertebrate mito
        if (lin.contains('echinodermata') || lin.contains('platyhelminthes'))  return 9
        if (lin.contains('cnidaria') || lin.contains('porifera'))              return 4
        return 5    // invertebrate mito (arthropods, molluscs, nematodes, ...)
    }
    return 1        // standard-code fallback
}

// Canonical telomere repeat, C-rich strand to match the existing 'CCCTAA' default.
// Only well-characterised clades mapped; else the vertebrate-style fallback. Tools search
// both strands, but VERIFY this convention suits hifiasm --telo-m / tidk, and extend freely.
def telomereMotifFor(Map r) {
    def lin = (r.lineage ?: '').toString().toLowerCase()
    if (lin.contains('viridiplantae'))   return 'CCCTAAA'   // plants: (TTTAGGG)n
    if (lin.contains('metazoa')) {
        if (lin.contains('vertebrata') || lin.contains('craniata')) return 'CCCTAA'   // (TTAGGG)n
        if (lin.contains('insecta'))     return 'CCTAA'     // many insects: (TTAGG)n
        if (lin.contains('nematoda'))    return 'GCCTAA'    // e.g. C. elegans: (TTAGGC)n
        return 'CCCTAA'
    }
    return 'CCCTAA'
}

// ── GetOrganelle: short-read organelle targets + per-organelle defaults ──────────────
// organelleTypesFor takes either an already-resolved kingdom flag ('plant'/'animal'/'fungi'
// /'other', as stored in ch_taxonomy) or a raw lineage/kingdom map (falls back to kingdomFlag).
def organelleTypesFor(Map r) {
    def k = (r.kingdom ?: '').toString().toLowerCase()
    if( !(k in ['plant','animal','fungi','other']) ) k = kingdomFlag(r)
    switch( k ) {
        case 'plant':  return ['embplant_pt', 'embplant_mt']   // plastid + plant mito
        case 'fungi':  return ['fungus_mt']
        case 'animal': return ['animal_mt']
        default:       return ['animal_mt']                     // conservative: mito only
    }
}

// GetOrganelle -R (max extension rounds). Plastid converges fast; plant mito needs more.
def getorganelleRecursionFor(String organelle) {
    switch( organelle ) {
        case 'embplant_pt': return 15
        case 'embplant_mt': return 30
        default:            return 10
    }
}

// GetOrganelle -k ladder. Full ladder to 127 for the large embryophyte organelles.
def getorganelleKmersFor(String organelle) {
    switch( organelle ) {
        case 'embplant_pt':
        case 'embplant_mt': return '21,45,65,85,105,127'
        default:            return '21,45,65,85,105'
    }
}

// Target organelle coverage for --reduce-reads-for-coverage. The single most important lever:
// deep WGS over-covers organelles (30,000x observed) into error hairballs; capping to ~50-100x
// gives a clean plastome/mito-sized graph in minutes. NEVER pass 'inf'.
def getorganelleCoverageFor(String organelle) {
    switch( organelle ) {
        case 'embplant_pt': return 100
        case 'embplant_mt': return 50
        default:            return 100
    }
}