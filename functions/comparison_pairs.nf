/* Finalized assemblies only; biological identity comes from metadata, not filenames. */
def comparisonPlan(records, mode) {
    if (!(mode in ['within_sample','all'])) error "Invalid pairwise_alignment_mode '${mode}'"
    def ordered=records.sort { a,b -> a[0].id <=> b[0].id }
    def seen=[] as Set
    def speciesBySample=[:]
    ordered.each { meta,fa,fai ->
        if (!meta.id || !meta.sample || !meta.taxid || !fa || !fai)
            error 'Final comparison input lacks id/sample/taxid/FASTA/FAI'
        if (!(meta.id.toString() ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) error 'Unsafe final comparison identity'
        if (!seen.add(meta.id.toString())) error "Duplicate final comparison identity '${meta.id}'"
        def taxid=meta.taxid.toString()
        if (speciesBySample.containsKey(meta.sample) && speciesBySample[meta.sample]!=taxid)
            error "Sample '${meta.sample}' has conflicting species in final comparisons"
        speciesBySample[meta.sample]=taxid
    }
    def pairs=[]
    def inventory=[]
    def outputNames=[] as Set
    for (int i=0;i<ordered.size();i++) {
        for (int j=i+1;j<ordered.size();j++) {
            def (left,fa1,fai1)=ordered[i]
            def (right,fa2,fai2)=ordered[j]
            if (left.taxid.toString()!=right.taxid.toString()) continue
            if (mode=='within_sample' && left.sample!=right.sample) continue
            if (!outputNames.add("${left.id}_vs_${right.id}".toString())) error 'Comparison output filename collision'
            pairs << tuple(left.id,fa1,fai1,right.id,fa2,fai2)
            inventory << [left.taxid.toString(),left.sample,right.sample,left.id,right.id,
                          fa1.toString(),fa2.toString(),'final',mode]
        }
    }
    [pairs:pairs,inventory:inventory]
}
