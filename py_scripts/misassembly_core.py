"""Cohort misassembly discovery and source-bound evidence primitives.

This module contains no sample names or known breakpoint coordinates.
"""
import csv
import gzip
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from chimera_tracks import tracks
from chimera_controls import molecule_chains, qualified_molecule

POLICY='cohort-review-v1'


def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as handle:
        for block in iter(lambda:handle.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def table(path):
    with open(path,encoding='utf-8') as handle:
        return list(csv.DictReader((l for l in handle if l.strip() and not l.startswith('#')),delimiter='\t'))


def write_table(path,rows,fields):
    with open(path,'w',newline='',encoding='utf-8') as handle:
        w=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)


def chromosome(name):
    return name.split('_')[0] if name.startswith('chr') and '+' not in name else None


def label_catalog(records):
    catalog={}
    for record in records:
        rows=table(record['name_map']) if Path(record['name_map']).name!='NO_HARMONIZE' else []
        labels={r['old_name']:chromosome(r['new_name']) for r in rows if chromosome(r['new_name'])}
        record=dict(record,chromosome_labels=labels)
        catalog[record['id']]=dict(record,scope=[r['old_name'] for r in rows if r.get('chromosome_member')=='yes'],
            scope_lengths={r['old_name']:int(r['length']) for r in rows if str(r.get('length','')).isdigit()},
            scope_unresolved=any(r.get('chromosome_member')=='unresolved' for r in rows))
    return catalog


def cohort_labels(catalog,alignments):
    """Resolve composite blocks using independent, plainly labelled cohort chromosomes.

    No iterative voting: composite-derived identities never label other composites.
    Both haplotypes of one individual count once, and conflicting identities stay blank.
    """
    from chimera_intervals import regions
    observations=defaultdict(lambda:defaultdict(lambda:defaultdict(list)))
    for pair in alignments:
        for own_id,peer_id,reverse in ((pair['a'],pair['b'],False),(pair['b'],pair['a'],True)):
            own,peer=catalog[own_id],catalog[peer_id]
            if own['taxid']!=peer['taxid']:raise ValueError('Cross-species cohort labels')
            if own['sample']==peer['sample'] or not peer['eligible']:continue
            unlabelled=set(own['scope'])-set(own['chromosome_labels'])
            if not unlabelled:continue
            grouped=defaultdict(list)
            for f in read_paf(pair['path'],reverse,unlabelled):
                grouped[f[0]].append(f)
            for scaffold,hits in grouped.items():
                for block in tracks(hits,peer['chromosome_labels'],mode='chromosome',min_mapq=20,minimum_target=0,compact=True):
                    if block['chrom']:observations[own_id][scaffold][peer['sample']].append((block['lo'],block['hi'],block['chrom']))
    for own_id,scaffolds in observations.items():
        own=catalog[own_id]
        available={v['sample'] for v in catalog.values() if v['taxid']==own['taxid'] and v['sample']!=own['sample'] and v['eligible']}
        minimum=min(2,len(available));own['composite_label_minimum_individuals']=minimum
        for scaffold,individuals in scaffolds.items():
            events=defaultdict(list)
            for sample,values in individuals.items():
                for lo,hi,chrom in regions(values,1):events[lo].append((sample,chrom,1));events[hi].append((sample,chrom,-1))
            active=defaultdict(Counter);coordinates=sorted(events);labels=[]
            for i,lo in enumerate(coordinates[:-1]):
                for sample,chrom,change in events[lo]:
                    active[sample][chrom]+=change
                    if active[sample][chrom]==0:del active[sample][chrom]
                votes=[next(iter(c)) for c in active.values() if len(c)==1]
                if len(votes)>=minimum and len(set(votes))==1:
                    block=dict(lo=lo,hi=coordinates[i+1],chrom=votes[0],individuals=len(votes),source_individuals=sorted(sample for sample,c in active.items() if len(c)==1))
                    if labels and labels[-1]['hi']==lo and all(labels[-1][k]==block[k] for k in ('chrom','individuals','source_individuals')):
                        labels[-1]['hi']=block['hi']
                    else:labels.append(block)
            own['chromosome_labels'][scaffold]=labels
        own['block_label_source']='CIGAR-projected cohort chromosomes; independent-individual agreement; conflicts and unmapped bases unlabelled'
    return catalog


def independent_labels(peer,focal_sample):
    """Do not count a focal individual's propagated label as confirmation of itself."""
    return {name:value if isinstance(value,str) else [b for b in value if len(set(b.get('source_individuals',[]))-{focal_sample})>=peer.get('composite_label_minimum_individuals',2)]
            for name,value in peer['chromosome_labels'].items()}


def invert_paf(f):
    """Swap query/target coordinates and CIGAR, including reverse strand."""
    out=[f[5],f[6],f[7],f[8],f[4],f[0],f[1],f[2],f[3],f[9],f[10],f[11]]
    cigar=next((v[5:] for v in f[12:] if v.startswith('cg:Z:')),None)
    if cigar is None:raise ValueError('Cohort alignment must retain CIGAR')
    operations=re.findall(r'(\d+)([MIDNSHP=X])',cigar)
    if f[4]=='-':operations.reverse()
    out.append('cg:Z:'+''.join(n+{'I':'D','D':'I'}.get(op,op) for n,op in operations))
    return out


def read_paf(path,reverse=False,query_names=None):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as handle:
        for line in handle:
            if not line.strip():continue
            f=line.rstrip().split('\t')
            if len(f)<12:raise ValueError('Malformed cohort PAF')
            # Discard out-of-scope records before expensive reverse CIGAR parsing.
            if query_names is not None and f[5 if reverse else 0] not in query_names:continue
            yield invert_paf(f) if reverse else f


def peer_arms(track,minimum=1000000,maximum_gap=1000000):
    """Aggregate assigned bases, preserving gaps; minor islands remain auditable."""
    arms=[]
    for block in track:
        if not block.get('chrom'):continue
        if arms and arms[-1]['chrom']==block['chrom'] and block['lo']-arms[-1]['end']<=maximum_gap:
            arms[-1]['end']=block['hi'];arms[-1]['aligned_bp']+=block['hi']-block['lo']
        else:arms.append(dict(chrom=block['chrom'],start=block['lo'],end=block['hi'],aligned_bp=block['hi']-block['lo']))
    minor=[a for a in arms if a['aligned_bp']<minimum]
    major=[]
    for arm in arms:
        if arm['aligned_bp']<minimum:continue
        if major and major[-1]['chrom']==arm['chrom'] and arm['start']-major[-1]['end']<=maximum_gap:
            major[-1]=dict(major[-1],end=arm['end'],aligned_bp=major[-1]['aligned_bp']+arm['aligned_bp'])
        else:major.append(dict(arm))
    events=[dict(left=a['chrom'],right=b['chrom'],lo=a['end'],hi=b['start'],unassigned_bp=b['start']-a['end'],
                 left_arm_bp=a['aligned_bp'],right_arm_bp=b['aligned_bp']) for a,b in zip(major,major[1:]) if a['chrom']!=b['chrom']]
    return events,minor


def discover(peer_tracks,sample,minimum=1000000,cluster_distance=250000):
    observations=[];minor=[]
    for peer,scaffolds in peer_tracks:
        for scaffold,track in scaffolds.items():
            events,islands=peer_arms(track,minimum)
            minor.extend(dict(assembly=peer['id'],scaffold=scaffold,**v) for v in islands)
            if peer['sample']==sample or not peer.get('eligible'):continue
            observations.extend(dict(peer=peer['id'],sample=peer['sample'],scaffold=scaffold,**v) for v in events)
    observations.sort(key=lambda v:(v['scaffold'],v['left'],v['right'],v['lo']))
    groups=[]
    for observation in observations:
        matches=[g for g in groups if (g['scaffold'],g['left'],g['right'])==(observation['scaffold'],observation['left'],observation['right']) and
                 observation['lo']<=g['hi']+cluster_distance and observation['hi']>=g['lo']-cluster_distance]
        if matches:
            group=matches[-1];group['lo']=min(group['lo'],observation['lo']);group['hi']=max(group['hi'],observation['hi']);group['observations'].append(observation)
        else:groups.append(dict(scaffold=observation['scaffold'],left=observation['left'],right=observation['right'],lo=observation['lo'],hi=observation['hi'],observations=[observation]))
    groups.sort(key=lambda g:(g['scaffold'],g['lo'],g['right']))
    for i,group in enumerate(groups,1):
        group['id']='C%02d'%i;group['independent_individuals']=sorted({r['sample'] for r in group['observations']})
    return groups,minor


def summarize_side(track,lo,hi,minimum=50000):
    counts=Counter()
    for b in track:
        if b.get('chrom'):counts[b['chrom']]+=max(0,min(hi,b['hi'])-max(lo,b['lo']))*b.get('aligned_bp',b['hi']-b['lo'])/max(1,b['hi']-b['lo'])
    total=sum(counts.values());identity,n=max(counts.items(),key=lambda v:v[1]) if counts else (None,0)
    return dict(chrom=identity,aligned_bp=total,dominance=n/max(1,total),qualified=total>=minimum and n/max(1,total)>=.9)


def coarse_tracks(track,width=100000):
    counts=defaultdict(Counter)
    for block in track:
        if not block.get('chrom'):continue
        for index in range(block['lo']//width,(block['hi']-1)//width+1):
            counts[index][block['chrom']]+=min(block['hi'],(index+1)*width)-max(block['lo'],index*width)
    result=[]
    for index,values in sorted(counts.items()):
        chrom,n=max(values.items(),key=lambda v:v[1]);total=sum(values.values())
        if total>=25000 and n/total>=.9:result.append(dict(lo=index*width,hi=(index+1)*width,chrom=chrom,aligned_bp=total))
    return result


def gap_spanning(path,lo,hi,mapq=20,anchor=1000):
    """A known N gap may be represented by a deletion; other large indels split support."""
    seen=set();spans=set()
    with open(path) as handle:
        for line in handle:
            if line.startswith('@'):continue
            f=line.rstrip().split('\t')
            if len(f)<11 or int(f[4])<mapq:continue
            copy=list(f);copy[4]=str(max(mapq,30))
            if not qualified_molecule(copy) or f[0] in seen:continue
            seen.add(f[0]);pos=int(f[3])-1;begin=None;end=None;chains=[]
            for n,op in re.findall(r'(\d+)([MIDNSHP=X])',f[5]):
                n=int(n)
                if op in 'M=X':
                    if begin is None:begin=pos
                    pos+=n;end=pos
                elif op in 'IDN':
                    expected=op=='D' and lo<=pos and pos+n<=hi
                    if op=='N' or (n>50 and not expected):
                        if begin is not None:chains.append((begin,end))
                        begin=end=None
                    if op in 'DN':pos+=n
            if begin is not None:chains.append((begin,end))
            if any(start<=lo-anchor and stop>=hi+anchor for start,stop in chains):spans.add(f[0])
    return len(spans)


def peer_test(peer,track,lo,hi,flank=500000):
    left=summarize_side(track,max(0,lo-flank),lo);right=summarize_side(track,hi,hi+flank)
    return dict(peer=peer['id'],sample=peer['sample'],eligible=peer.get('eligible',False),left=left,right=right,
                qualified=left['qualified'] and right['qualified'])


def summarize_peers(tests,sample,expected):
    groups=defaultdict(set)
    for t in tests:
        if t['sample']==sample or not t['eligible'] or not t['qualified']:continue
        groups[t['sample']].add((t['left']['chrom'],t['right']['chrom']))
    separate=[];continuous=[];conflicting=[]
    for individual,pairs in groups.items():
        if pairs=={tuple(expected)}:separate.append(individual)
        elif len(pairs)==1 and next(iter(pairs))[0]==next(iter(pairs))[1]:continuous.append(individual)
        else:conflicting.append(individual)
    return dict(separate=sorted(separate),continuous=sorted(continuous),conflicting=sorted(conflicting))


def support_grid(path,lo,hi,step=1000,anchor=1000,mapq=20):
    points=sorted(set([lo,hi]+list(range(lo,hi+1,step))));counts=[0]*len(points)
    import bisect
    seen=set()
    with open(path) as handle:
        for line in handle:
            if line.startswith('@'):continue
            f=line.rstrip().split('\t')
            if len(f)<11 or int(f[1])&(4|256|512|1024|2048) or int(f[4])<mapq:continue
            # Apply the same identity/CIGAR screen at both MAPQ tiers.
            copy=list(f);copy[4]=str(max(30,mapq))
            if not qualified_molecule(copy) or f[0] in seen:continue
            seen.add(f[0]);indices=set()
            for start,end in molecule_chains(f):
                indices.update(range(bisect.bisect_left(points,start+anchor),bisect.bisect_right(points,end-anchor)))
            for index in indices:counts[index]+=1
    weak=[];start=None
    for i,(point,n) in enumerate(zip(points,counts)):
        if n<2 and start is None:start=point
        if start is not None and (n>=2 or i==len(points)-1):
            end=point if n>=2 else min(hi,point+step)
            weak.append(dict(start=start,end=end,minimum=min(counts[points.index(start):i+1])));start=None
    fraction=sum(n>=2 for n in counts)/len(counts);longest=max((v['end']-v['start'] for v in weak),default=0)
    state='continuous' if min(counts)>=2 else 'mostly_continuous' if fraction>=.95 and longest<=3000 else 'patchy'
    return dict(mapq=mapq,anchor_bp=anchor,step_bp=step,minimum=min(counts),supported_fraction=fraction,
                longest_weak_run_bp=longest,state=state,weak_runs=weak,probes=[dict(position=p,molecules=n) for p,n in zip(points,counts)])


def contact_result(focal,controls):
    def ratio(row):
        denominator=math.sqrt(row.get('left_within',0)*row.get('right_within',0))
        return row.get('cross',0)/denominator if denominator else None
    r=ratio(focal)
    usable=[c for c in controls if min(c.get('left_within',0),c.get('right_within',0))>=100 and ratio(c) is not None]
    # One appropriate control population suffices; report coverage matching explicitly.
    matched=[c for c in usable if all(.1<=c.get(k,0)/max(1,focal.get(k,0))<=10 for k in ('left_within','right_within'))]
    values=sorted(ratio(c) for c in matched)
    floor=values[int(.1*(len(values)-1))] if values else None
    upper=(focal.get('cross',0)+3)/math.sqrt(focal['left_within']*focal['right_within']) if r is not None else None
    calibrated=min(focal.get('left_within',0),focal.get('right_within',0))>=100 and len(matched)>=5
    return dict(raw_counts=dict(focal),ratio=r,matched_controls=len(matched),comparison_floor=floor,
                upper_ratio=upper,calibrated=calibrated,loss=bool(calibrated and floor>0 and upper<.25*floor),control_ratios=values,
                control_ids=[c['id'] for c in matched if 'id' in c],minimum_within_pairs=100,minimum_matched_controls=5,depth_matching_factor=10)


def read_depth(path,lo,hi,width=1000,mapq=20):
    """Mean aligned primary-molecule depth; large deletions are not covered bases."""
    values=[0.0]*math.ceil((hi-lo)/width);seen=set()
    with open(path) as handle:
        for line in handle:
            if line.startswith('@'):continue
            f=line.rstrip().split('\t')
            if len(f)<11 or int(f[4])<mapq:continue
            copy=list(f);copy[4]='30'
            if not qualified_molecule(copy) or f[0] in seen:continue
            seen.add(f[0]);pos=int(f[3])-1
            for n,op in re.findall(r'(\d+)([MIDNSHP=X])',f[5]):
                n=int(n)
                if op in 'M=X':
                    start,end=max(lo,pos),min(hi,pos+n)
                    if start<end:
                        for i in range((start-lo)//width,(end-lo-1)//width+1):
                            a,b=lo+i*width,min(hi,lo+(i+1)*width)
                            values[i]+=max(0,min(end,b)-max(start,a))/(b-a)
                if op in 'MDN=X':pos+=n
    return [dict(start=lo+i*width,end=min(hi,lo+(i+1)*width),depth=v) for i,v in enumerate(values)]


def matched_control_sites(intervals,gaps,lengths,flank=250000,count=24):
    """Match separation as well as flank width; never mix distant and point assays."""
    controls={};links={}
    suspects=list(intervals.values())
    for key,target in intervals.items():
        width=target['hi']-target['lo'];gap_target='cut_bp' in target
        if gap_target:
            sites=[g for g in gaps if .5*width<=g['hi']-g['lo']<=2*width]
        else:
            sites=[dict(scaffold=sc,lo=point,hi=point+width) for sc,length in sorted(lengths.items())
                   for point in range(flank,length-width-flank,2*flank+width)]
        available=[s for s in sites if s['lo']>=flank and s['hi']+flank<=lengths[s['scaffold']] and
                   not any(v['scaffold']==s['scaffold'] and s['lo']<v['hi']+flank and s['hi']>v['lo']-flank for v in suspects) and
                   (gap_target or not any(g['scaffold']==s['scaffold'] and s['lo']-flank<g['hi'] and s['hi']+flank>g['lo'] for g in gaps))]
        available.sort(key=lambda s:(s['scaffold']!=target['scaffold'],abs(s['lo']-target['lo']),s['scaffold'],s['lo']))
        chosen=[]
        for site in available:
            if any(site['scaffold']==s['scaffold'] and abs(site['lo']-s['lo'])<2*flank for s in chosen):continue
            chosen.append(site)
            if len(chosen)>=count:break
        links[key]=[]
        for site in chosen:
            role='gap_control' if gap_target else 'continuous_control'
            cid=role+'_'+hashlib.sha256(('%s:%d:%d'%(site['scaffold'],site['lo'],site['hi'])).encode()).hexdigest()[:12]
            controls[cid]=dict(site,id=cid,role=role);links[key].append(cid)
    return controls,links


def recommend(event):
    peers=event['peers'];grid=event.get('hifi',{});options=event.get('cut_options',[])
    if len(peers['continuous'])>=2 and not peers['separate'] and not peers['conflicting']:
        return 'RETAIN','Other individuals support one chromosome on both sides.'
    eligible=[c for c in options if len(c['peers']['separate'])>=2 and not c['peers']['continuous'] and not c['peers']['conflicting'] and
              'spanning' in c and c['spanning']<2 and min(c.get('flank_molecules',[0,0]))>=2 and not c.get('native_continuity') and c.get('bridge_consistent',True) and
              sum(v['loss'] for v in c.get('contacts',[]))>=2]
    if len(eligible)==1:
        event['recommended_cut']=eligible[0]
        return 'CUT','Independent chromosome separation and contact loss in two libraries support this exact gap cut.'
    if any(c.get('spanning',0)>=2 and c['lo']<=event['lo']<=event['hi']<=c['hi'] for c in options):
        return 'RETAIN','At least two HiFi molecules directly span the verified gap containing this transition.'
    if grid.get('state') in ('continuous','mostly_continuous'):
        return 'RETAIN','Local HiFi molecules support continuity across the transition'+(' with short weak pockets.' if grid['state']=='mostly_continuous' else '.')
    if peers['conflicting'] or peers['continuous']:return 'REVIEW','Chromosome comparisons conflict; leave intact while reviewing the conflicting evidence.'
    if len(peers['separate'])>=2 or len(event.get('independent_individuals',[]))>=2:
        return 'SUSPECT','The chromosome transition is supported, but an exact unsupported join has not been established.'
    return 'REVIEW','This is a screening signal with insufficient independent confirmation.'
