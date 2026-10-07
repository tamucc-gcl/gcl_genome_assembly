import argparse,csv,html
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--table',action='append',required=True);a=p.parse_args()
rows=[];fields=None;ids=set()
for path in a.table:
    with open(path) as h:
        reader=csv.DictReader(h,delimiter='\t')
        if fields is None:fields=reader.fieldnames
        if reader.fieldnames!=fields:raise ValueError('Review table headers differ')
        for row in reader:
            identity=(row['assembly'],row['id'])
            if identity in ids:raise ValueError('Duplicate candidate identity')
            ids.add(identity);rows.append(row)
rows.sort(key=lambda r:(r['assembly'],r['scaffold'],r['id']))
with open('chimera_review.tsv','w',newline='',encoding='utf-8') as h:
    w=csv.DictWriter(h,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
assemblies=sorted({r['assembly'] for r in rows})
content='<!doctype html><meta charset="utf-8"><title>Chimera review index</title><h1>Chimera evidence and manual cuts</h1><p>Automated cutting is deferred. Download or copy chimera_review.tsv, keep unselected rows, and mark selected=YES only after reviewing a cut. Supply reviewer and reason. Run again with --chimera_break /absolute/path/to/edited.tsv -resume. Every cut is validated against the original pre-finishing FASTA.</p><a href="chimera_review.tsv">Editable review file</a><ul>'
for assembly in assemblies:content+='<li><a href="'+html.escape(assembly+'.review/report.html')+'">'+html.escape(assembly)+'</a></li>'
Path('index.html').write_text(content+'</ul>',encoding='utf-8')
