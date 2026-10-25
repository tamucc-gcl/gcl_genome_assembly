/* Broadcast species context to every sample without silently dropping missing keys. */
def attachRequiredContext(records, context, label) {
    def lookup = context.toList().map { rows ->
        def byTaxid = [:]
        rows.each { row ->
            def taxid = row[0]?.toString()
            if (!taxid || byTaxid.containsKey(taxid) || row.size() < 2 || row.drop(1).any { it == null })
                error "Invalid or duplicate ${label} context for taxid '${taxid}'"
            byTaxid[taxid] = row.drop(1)
        }
        byTaxid
    }
    records.map { row -> tuple(row[0]?.toString(), row) }.combine(lookup)
        .map { taxid, row, byTaxid ->
            if (!taxid || !byTaxid.containsKey(taxid))
                error "Missing required ${label} context for taxid '${taxid}'"
            tuple(*(row.toList() + byTaxid[taxid]))
        }
}
