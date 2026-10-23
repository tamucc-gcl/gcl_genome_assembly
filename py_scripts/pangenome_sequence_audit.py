"""Post-build FULL/CLIP node accounting and original-source coordinate ownership."""
import argparse
from collections import defaultdict
import gzip
import json
from pathlib import Path
import re

from pangenome_source_audit import read_table, write_table, union_bp


def graph_aliases(sources):
    result={}
    for r in sources:
        sample=r['sample'];base,dot,hap=sample.rpartition('.')
        alias=base+'#'+hap if dot and hap.isdigit() else sample+'#0'
        if alias in result and result[alias]!=sample:raise ValueError('Ambiguous production graph identity')
        result[alias]=sample
    return result


def resolve_sequence(sample,sequence,lengths):
    query=f'id={sample}|{sequence}'
    if query in lengths:return query,0,lengths[query]
    sequence=re.sub(r'#\d+$','',sequence)
    query=f'id={sample}|{sequence}'
    if query in lengths:return query,0,lengths[query]
    # Pinned hal2vg/vg can retain a fragment suffix in a sequence identifier.
    match=re.fullmatch(r'(.*?)(?:\[(\d+)-(\d+)\]|_sub_(\d+)_(\d+))',sequence)
    if match:
        query=f'id={sample}|{match[1]}';lo=int(match[2] or match[4]);hi=int(match[3] or match[5])
        if query in lengths and 0<=lo<hi<=lengths[query]:return query,lo,hi
    raise ValueError('Graph sequence absent from source manifest: '+query)


def subtract(left,right):
    result=[];j=0
    for lo,hi in sorted(left):
        while j<len(right) and right[j][1]<=lo:j+=1
        pos=lo;k=j
        while k<len(right) and right[k][0]<hi:
            start,end=right[k]
            if start>pos:result.append((pos,min(start,hi)))
            pos=max(pos,end)
            if pos>=hi:break
            k+=1
        if pos<hi:result.append((pos,hi))
    return result


def scan_graph(path,sources,graph,out):
    lengths={r['query']:int(r['source_bp']) for r in sources};aliases=graph_aliases(sources)
    intervals=defaultdict(list);nodes=set();node_bp=0;synthetic=0;records=0
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt',encoding='utf-8') as handle:
        for line in handle:
            fields=line.rstrip().split('\t')
            if fields[0]=='S':
                if fields[1] in nodes:raise ValueError('Duplicate graph segment ID')
                nodes.add(fields[1]);seq=fields[2]
                length=len(seq) if seq!='*' else next((int(t[5:]) for t in fields[3:] if t.startswith('LN:i:')),None)
                if length is None or length<0:raise ValueError('Graph segment lacks valid length')
                node_bp+=length
            elif fields[0]=='W':
                records+=1;sample=aliases.get('#'.join(fields[1:3]))
                if sample is None:
                    if fields[1].startswith(('_MINIGRAPH_','_GREF_')):synthetic+=1;continue
                    raise ValueError('Unrecognized biological graph walk sample: '+fields[1])
                query,offset,end=resolve_sequence(sample,fields[3],lengths)
                if not fields[4].isdigit() or not fields[5].isdigit():raise ValueError('Graph walk lacks source coordinates')
                lo,hi=int(fields[4])+offset,int(fields[5])+offset
                if not offset<=lo<=hi<=end:raise ValueError('Graph walk exceeds original source/fragment bounds')
                if hi>lo:intervals[query].append((lo,hi))
            elif fields[0]=='P':
                # Biological source accountability requires W coordinates from the pinned export.
                alias='#'.join(fields[1].split('#')[:2])
                if alias in aliases:raise ValueError('Biological P path lacks audited source W coordinates')
    if not nodes or not records:raise ValueError('Empty graph or missing walk records')
    rows=[]
    for r in sources:
        spans=intervals[r['query']];total=sum(hi-lo for lo,hi in spans);covered=union_bp(spans)
        if total!=covered:raise ValueError('Overlapping source walk coordinates: '+r['query'])
        rows.append(dict(graph=graph,**r,source_span_bp=covered,source_not_in_walk_spans_bp=int(r['source_bp'])-covered,
            walk_records=len(spans),overlap_bp=0,status='PASS'))
    write_table(out/(graph.lower()+'_source_coordinates.tsv'),rows,list(rows[0]))
    fields=['query','start_0based','end_exclusive']
    with gzip.open(out/(graph.lower()+'_source_intervals.tsv.gz'),'wt',encoding='utf-8',newline='') as handle:
        import csv
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader()
        for query,spans in sorted(intervals.items()):
            for lo,hi in sorted(spans):writer.writerow(dict(query=query,start_0based=lo,end_exclusive=hi))
    return dict(graph=graph,distinct_nodes=len(nodes),distinct_node_bp=node_bp,biological_source_span_bp=sum(r['source_span_bp'] for r in rows),
        source_bp=sum(int(r['source_bp']) for r in sources),synthetic_walk_records=synthetic,coordinate_integrity='PASS'),rows


def audit(full,clip,policy,out):
    policy=Path(policy);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    meta=json.loads((policy/'sequence_accounting.json').read_text());sources=read_table(policy/'source_lengths.tsv')
    stages=[];by={}
    for label,path in [('FULL',full),('CLIP',clip)]:
        stage,rows=scan_graph(path,sources,label,out);stages.append(stage);by[label]={r['query']:r for r in rows}
    partitions=list(policy.glob('split_*.sources.tsv'))
    if len(partitions)!=1:raise ValueError('Missing/ambiguous production partition source audit')
    assigned={r['query']:r for r in read_table(partitions[0])}
    if set(assigned)!=set(by['FULL']):raise ValueError('Source identities differ between graph and partition audit')
    regions=list(policy.glob('split_*.regions.tsv.gz'))
    if len(regions)!=1:raise ValueError('Missing production partition region coordinates')
    masks={stage:defaultdict(list) for stage in ('partition','FULL','CLIP')}
    ambiguous=json.loads((policy/'policy.json').read_text())['ambiguous_name']
    for r in read_table(regions[0]):
        if r['chromosome']!=ambiguous:masks['partition'][r['query']].append((int(r['start_0based']),int(r['end_exclusive'])))
    for stage in ('FULL','CLIP'):
        for r in read_table(out/(stage.lower()+'_source_intervals.tsv.gz')):
            masks[stage][r['query']].append((int(r['start_0based']),int(r['end_exclusive'])))
    rows=[];differences=[]
    for r in sources:
        q=r['query'];p=assigned[q];f=by['FULL'][q];c=by['CLIP'][q]
        changes={}
        for label,a,b in [('assigned_missing_from_FULL','partition','FULL'),
                          ('FULL_outside_assigned','FULL','partition'),('CLIP_outside_FULL','CLIP','FULL')]:
            spans=subtract(masks[a][q],sorted(masks[b][q]));changes[label+'_bp']=union_bp(spans)
            for lo,hi in spans:differences.append(dict(query=q,category=label,start_0based=lo,end_exclusive=hi))
        rows.append(dict(**r,partition_assigned_bp=p['assigned_bp'],partition_ambiguous_bp=p['ambiguous_bp'],
            partition_unreported_bp=p['unreported_bp'],full_source_span_bp=f['source_span_bp'],clip_source_span_bp=c['source_span_bp'],
            full_walk_records=f['walk_records'],clip_walk_records=c['walk_records'],**changes))
    write_table(out/'source_stage_accounting.tsv',rows,list(rows[0]))
    plain=out/'source_stage_discrepancies.tsv';write_table(plain,differences,['query','category','start_0based','end_exclusive'])
    with plain.open('rb') as src,gzip.open(str(plain)+'.gz','wb') as dst:
        import shutil
        shutil.copyfileobj(src,dst)
    plain.unlink()
    write_table(out/'distinct_graph_nodes.tsv',stages,list(stages[0]))
    discrepancy_totals={name:sum(r[name] for r in rows) for name in changes}
    (out/'sequence_graph_audit.json').write_text(json.dumps(dict(status='PASS',policy=meta['policy'],sources=len(sources),graphs=stages,
        representation_review='required' if any(discrepancy_totals.values()) else 'no_interval_discrepancies',
        coordinate_discrepancies=discrepancy_totals),indent=2)+'\n')
    text=['# Graph sequence acceptance checks','',
        'Original-source walk coordinates are in bounds and non-overlapping in both graphs. Each graph segment is counted once, '
        'regardless of how many haplotype paths traverse it. Graph-node totals and source-coordinate spans measure different things.','',
        '| Graph | Distinct nodes | Distinct graph bases | Biological source coordinate spans |','|---|---:|---:|---:|']
    for r in stages:text.append(f'| {r["graph"]} | {r["distinct_nodes"]:,} | {r["distinct_node_bp"]:,} | {r["biological_source_span_bp"]:,} |')
    text+=['','[Every source across preprocessing, FULL and CLIP](source_stage_accounting.tsv)','',
        f'Assigned source spans missing from FULL: {discrepancy_totals["assigned_missing_from_FULL_bp"]:,} bp; '
        f'FULL spans outside assigned partitions: {discrepancy_totals["FULL_outside_assigned_bp"]:,} bp; '
        f'CLIP spans outside FULL: {discrepancy_totals["CLIP_outside_FULL_bp"]:,} bp. '
        'Nonzero differences require review; exact coordinates are retained in `source_stage_discrepancies.tsv.gz`.','',
        'Walk spans do not establish that every base within the span survives alignment or clipping. '
        'Read this audit alongside clipping/path-length reports and `sequence_policy/sequence_accounting.md`. '
        'Preprocessing exclusions are explicit; original assemblies retain the excluded sequence.','']
    (out/'sequence_graph_audit.md').write_text('\n'.join(text),encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('full','clip','policy','out'):p.add_argument('--'+k,required=True)
    audit(**vars(p.parse_args()))
