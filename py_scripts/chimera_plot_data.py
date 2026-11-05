"""Lossless contact-matrix presentation export; no assay or decision changes."""
import gzip
import json
from pathlib import Path


def export_contacts(packet, event):
    import numpy as np
    root=Path(packet)/'junction'
    summary=root/'hic/summary.json'
    records=[]
    if summary.is_file():
        for row in json.loads(summary.read_text()).get('maps',[]):
            if row['event']!=event['id'] or row['library']=='POOLED' or row['resolution'] not in (10000,25000):continue
            with np.load(root/'hic'/row['file'],allow_pickle=False) as source:
                matrix=source['matrix'];lo=int(source['lo']);resolution=int(source['resolution'])
            if matrix.ndim!=2 or matrix.shape[0]!=matrix.shape[1] or matrix.shape[0]==0 or not np.isfinite(matrix).all() or (matrix<0).any():
                raise ValueError('Invalid retained contact matrix: '+row['file'])
            if lo!=row['lo'] or resolution!=row['resolution'] or not np.equal(matrix,np.floor(matrix)).all():
                raise ValueError('Contact metadata or counts disagree: '+row['file'])
            filename=f"{event['id']}.{row['library']}.{resolution}.plot.tsv.gz"
            i,j=np.nonzero(matrix)
            with gzip.open(root/filename,'wt',encoding='utf-8',newline='') as handle:
                handle.write('row\tcolumn\tcount\n')
                np.savetxt(handle,np.column_stack((i,j,matrix[i,j])),fmt='%d',delimiter='\t')
            records.append(dict(event=event['id'],library=row['library'],resolution=resolution,
                lo=lo,bins=len(matrix),file=filename,maximum=int(matrix.max(initial=0)),
                total=int(matrix.sum()),nonzero_cells=len(i),zero_cells=int(matrix.size-len(i))))
    (root/(event['id']+'.contact_plot_spec.json')).write_text(json.dumps(records,indent=2)+'\n',encoding='utf-8')
    return records
