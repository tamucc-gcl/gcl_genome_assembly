"""Validate/filter PAF and write synchronous, explicit alignment summaries."""
import argparse
import csv
import gzip
import math
from pathlib import Path
import sys

FIELDS=['ref_id','qry_id','preset','min_mapq','min_aln_bp','n_align','sum_aln_bp',
        'sum_match_bp','mean_identity','mean_dv','mean_mapq','dv_observed_alignments',
        'mapq_observed_alignments']


def parse_paf(line):
    f=line.rstrip('\r\n').split('\t')
    if len(f)<12:raise ValueError('PAF record has fewer than 12 columns')
    try:
        qlen,qs,qe=int(f[1]),int(f[2]),int(f[3])
        tlen,ts,te=int(f[6]),int(f[7]),int(f[8])
        matches,length,mapq=int(f[9]),int(f[10]),int(f[11])
    except ValueError as exc:raise ValueError('Invalid PAF integer') from exc
    if (not f[0] or not f[5] or f[4] not in ['+','-'] or not 0<=qs<qe<=qlen
            or not 0<=ts<te<=tlen or length<=0 or not 0<=matches<=length or not 0<=mapq<=255):
        raise ValueError('Invalid PAF coordinates/strand/counts')
    dv=None
    tags=[s for s in f[12:] if s.startswith('dv:')]
    if len(tags)>1:raise ValueError('Duplicate PAF divergence tag')
    if tags:
        if not tags[0].startswith('dv:f:'):raise ValueError('Invalid PAF divergence tag type')
        dv=float(tags[0][5:])
        if not math.isfinite(dv) or not 0<=dv<=1:raise ValueError('Invalid PAF divergence')
    return dict(qname=f[0],qlen=qlen,qs=qs,qe=qe,strand=f[4],tname=f[5],tlen=tlen,
                ts=ts,te=te,matches=matches,length=length,mapq=mapq,dv=dv)


def filter_paf(stream, output, qc, ref, query, preset, min_mapq, min_length):
    if not 0<=min_mapq<=254 or min_length<0:raise ValueError('Invalid alignment filter thresholds')
    n=length=matches=dv_n=dv_length=mapq_n=0
    dv_sum=mapq_sum=0.0
    with gzip.open(output,'wt',compresslevel=3) as out:
        for line in stream:
            if not line.strip():continue
            r=parse_paf(line)
            # MAPQ 255 means unavailable, not exceptionally reliable.
            if r['mapq']==255 and min_mapq>0:continue
            if (r['mapq']!=255 and r['mapq']<min_mapq) or r['length']<min_length:continue
            out.write(line.rstrip('\r\n')+'\n')
            n+=1;length+=r['length'];matches+=r['matches']
            if r['dv'] is not None:
                dv_n+=1;dv_length+=r['length'];dv_sum+=r['dv']*r['length']
            if r['mapq']!=255:mapq_n+=1;mapq_sum+=r['mapq']
    row=dict(zip(FIELDS,[ref,query,preset,min_mapq,min_length,n,length,matches,
        matches/length if length else 'NA',dv_sum/dv_length if dv_length else 'NA',
        mapq_sum/mapq_n if mapq_n else 'NA',dv_n,mapq_n]))
    with open(qc,'w',newline='') as out:
        w=csv.DictWriter(out,FIELDS,delimiter='\t');w.writeheader();w.writerow(row)


def collect(paths,output):
    rows=[]
    for path in sorted(paths):
        with open(path) as handle:
            reader=csv.DictReader(handle,delimiter='\t')
            if reader.fieldnames!=FIELDS:raise ValueError(f'Invalid QC columns: {path}')
            data=list(reader)
            if len(data)!=1:raise ValueError(f'Expected one QC record: {path}')
            rows.extend(data)
    if len({(r['ref_id'],r['qry_id']) for r in rows})!=len(rows):raise ValueError('Duplicate comparison QC identity')
    with open(output,'w',newline='') as out:
        w=csv.DictWriter(out,FIELDS,delimiter='\t');w.writeheader()
        w.writerows(sorted(rows,key=lambda r:(r['ref_id'],r['qry_id'])))


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    f=sub.add_parser('filter')
    for name in ['output','qc','ref','query','preset']:f.add_argument('--'+name,required=True)
    f.add_argument('--min-mapq',type=int,required=True);f.add_argument('--min-length',type=int,required=True)
    c=sub.add_parser('collect');c.add_argument('--directory',required=True);c.add_argument('--output',required=True)
    a=p.parse_args()
    try:
        if a.command=='filter':filter_paf(sys.stdin,a.output,a.qc,a.ref,a.query,a.preset,a.min_mapq,a.min_length)
        else:collect(Path(a.directory).glob('*.tsv'),a.output)
    except (ValueError,OSError) as exc:p.exit(65,str(exc)+'\n')


if __name__=='__main__':main()
