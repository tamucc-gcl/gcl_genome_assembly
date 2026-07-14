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