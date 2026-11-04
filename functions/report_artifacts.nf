// Published artifact locations are explicit; task scratch paths are never report URLs.
def graphReviewRows(source) {
    def lines=source.readLines()
    if (!lines) throw new IllegalArgumentException('Empty graph source audit table')
    def columns=lines[0].split('\t',-1) as List
    def flag=columns.indexOf('retention_review')
    if (flag<0) throw new IllegalArgumentException('Graph source audit lacks retention_review')
    ([lines[0]]+lines.drop(1).findAll { line ->
        def fields=line.split('\t',-1)
        fields.size()>flag && fields[flag]=='yes'
    }).join('\n')+'\n'
}
def graphReportPath(taxid, name) {
    def folders=[
        'clipping_report.md':'representation', 'source_report.md':'representation',
        'sequence_graph_audit.md':'representation', 'comparison_report.md':'comparison',
        'regional_report.md':'regional', 'regional_figures.md':'regional/figures',
        'sharing_report.md':'sharing', 'sharing_figures.md':'figures',
        'variant_report.md':'variants', 'pangenome_build_report.md':'',
        'pangenome_input_report.md':'', 'pangenome_report.md':''
    ]
    if (!folders.containsKey(name.toString())) throw new IllegalArgumentException("Register published report location for ${name}")
    ['pangenome',taxid.toString(),folders[name.toString()],name.toString()].findAll { it }.join('/')
}

def rebaseReportMarkdown(text, publishedPath) {
    def parent=java.nio.file.Paths.get(publishedPath).parent
    def target={ url ->
        if (url =~ /^(?:[a-zA-Z][a-zA-Z0-9+.-]*:|\/)/) return url
        if (url.startsWith('#')) return publishedPath+url
        parent.resolve(url).normalize().toString().replace('\\','/')
    }
    text.replaceAll(/(!?\[[^\]]*\]\()([^\s)]+)(\))/) { all,a,url,b -> a+target(url)+b }
        .replaceAll(/(<img\b[^>]*\bsrc=["'])([^"']+)(["'])/) { all,a,url,b -> a+target(url)+b }
        .replaceAll(/(?m)^(#{1,3}) /) { all,hashes -> '###'+hashes+' ' }
}
