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

POLICY='cohort-review-v3'


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
            arms[-1]['end']=block['hi'];arms[-1]['aligned_bp']+=block.get('aligned_bp',block['hi']-block['lo'])
        else:arms.append(dict(chrom=block['chrom'],start=block['lo'],end=block['hi'],aligned_bp=block.get('aligned_bp',block['hi']-block['lo'])))
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


def screen_transition_flanks(group,local,sample):
    tests=[peer_test(dict(id=t['peer'],sample=t['sample'],eligible=t['eligible']),t['blocks'],group['lo'],group['hi']) for t in local]
    votes=summarize_peers(tests,sample,(group['left'],group['right']))
    single=len({o['sample'] for o in group['observations']})==1
    opposed=single and len(votes['continuous'])>=2 and len(votes['continuous'])>len(votes['separate'])+len(votes['conflicting'])
    return dict(opposed=opposed,peer_tests=tests,peers=votes)


def discover(peer_tracks,sample,minimum=1000000,cluster_distance=250000):
    observations=[];minor=[]
    width=max(1,min(100000,minimum//10))
    for peer,scaffolds in peer_tracks:
        for scaffold,track in scaffolds.items():
            # Use regional dominance for chromosome-scale topology, preserving
            # exact alignments in peer tracks for the subsequent local assays.
            events,islands=peer_arms(coarse_tracks(track,width),minimum)
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
    retained=[]
    for group in groups:
        local=[dict(peer=p['id'],sample=p['sample'],eligible=p['eligible'],blocks=m.get(group['scaffold'],[])) for p,m in peer_tracks]
        regional=consensus_regions(local,sample,width)
        transitions,_=peer_arms(regional,minimum)
        matches=[t for t in transitions if (t['left'],t['right'])==(group['left'],group['right']) and t['lo']<=group['hi']+cluster_distance and t['hi']>=group['lo']-cluster_distance]
        if matches:
            transition=min(matches,key=lambda t:abs(t['lo']-group['lo'])+abs(t['hi']-group['hi']))
            group['screening_range']=dict(lo=group['lo'],hi=group['hi'])
            lo,hi=transition['lo'],transition['hi'];resolution=width
            focused=[dict(t,blocks=[dict(b,lo=max(b['lo'],max(0,lo-2*width)),hi=min(b['hi'],hi+2*width)) for b in t['blocks'] if b['hi']>max(0,lo-2*width) and b['lo']<hi+2*width]) for t in local]
            fine=consensus_regions(focused,sample,max(1,width//100))
            left=[r for r in fine if r['chrom']==group['left'] and r['lo']<=hi]
            right=[r for r in fine if r['chrom']==group['right'] and r['hi']>=lo]
            if left and right:
                left_end=max(r['hi'] for r in left)
                right_start=min(r['lo'] for r in right if r['hi']>=left_end) if any(r['hi']>=left_end for r in right) else hi
                if left_end<=right_start:lo,hi=left_end,right_start;resolution=max(1,width//100)
            group.update(lo=lo,hi=hi,localization='Independent-individual regional chromosome transition',localization_resolution_bp=resolution)
        else:
            screening=screen_transition_flanks(group,local,sample)
            same_region=any(r['lo']<=group['lo'] and r['hi']>=group['hi'] and r['chrom'] in (group['left'],group['right']) for r in regional)
            # Missing bins do not invalidate qualified same-chromosome flanks.
            # This rejects the fusion hypothesis, not every possible local error.
            if same_region or screening['opposed']:
                minor.append(dict(group,kind='rejected_peer_transition',reason='Independent local assignments support the same chromosome on both sides',peer_tests=screening['peer_tests'],peers=screening['peers']))
                continue
            group['localization']='Peer transition; cohort boundary unresolved'
        group.setdefault('localization_resolution_bp',width)
        retained.append(group)
    groups=retained
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
            counts[index][block['chrom']]+=(min(block['hi'],(index+1)*width)-max(block['lo'],index*width))*block.get('aligned_bp',block['hi']-block['lo'])/max(1,block['hi']-block['lo'])
    result=[]
    for index,values in sorted(counts.items()):
        chrom,n=max(values.items(),key=lambda v:v[1]);total=sum(values.values())
        if total>=width*.25 and n/total>=.9:result.append(dict(lo=index*width,hi=(index+1)*width,chrom=chrom,aligned_bp=n))
    return result


def consensus_regions(tracks,sample,width=10000):
    """Regional identities, counting each independent individual once.

    Conflicting sister haplotypes abstain. Require two individuals when available
    and a strict majority of informative individuals; missing regions stay blank.
    """
    eligible=[t for t in tracks if t['eligible'] and t['sample']!=sample]
    minimum=min(2,len({t['sample'] for t in eligible}))
    from statistics import median
    bins=defaultdict(lambda:defaultdict(set));coverage={}
    for track in eligible:
        for block in coarse_tracks(track['blocks'],width):
            bins[block['lo']][track['sample']].add(block['chrom'])
            key=(block['lo'],track['sample'],block['chrom'])
            coverage[key]=max(coverage.get(key,0),block['aligned_bp'])
    result=[]
    for lo,individuals in sorted(bins.items()):
        votes=Counter(next(iter(labels)) for labels in individuals.values() if len(labels)==1)
        if not votes:continue
        chrom,n=max(votes.items(),key=lambda item:item[1])
        if n<minimum or n<=sum(votes.values())/2:continue
        # Conflicting haplotypes count as uncertainty, rather than a supporting vote.
        if n<=len(individuals)/2:continue
        assigned=median(coverage[lo,individual,chrom] for individual,labels in individuals.items() if labels=={chrom})
        if result and result[-1]['hi']==lo and result[-1]['chrom']==chrom:
            result[-1]['hi']=lo+width;result[-1]['aligned_bp']+=assigned
        else:result.append(dict(lo=lo,hi=lo+width,chrom=chrom,aligned_bp=assigned))
    return result


def gap_association(event,option,sample):
    """Keep unmeasured association separate from observed chromosome contradiction."""
    lo=min(event['lo'],option['lo']);hi=max(event['hi'],option['hi'])
    evidence=[]
    for track in event['peer_tracks']:
        if not track['eligible'] or track['sample']==sample:continue
        blocks=[b for b in coarse_tracks(track['blocks'],1000) if b['lo']<hi and b['hi']>lo]
        bases=sum(max(0,min(hi,b['hi'])-max(lo,b['lo']))*b['aligned_bp']/(b['hi']-b['lo']) for b in blocks)
        order=[]
        for b in blocks:
            if not order or order[-1]!=b['chrom']:order.append(b['chrom'])
        evidence.append(dict(peer=track['peer'],sample=track['sample'],assigned_bp=bases,order=order,
            consistent=order in ([event['left']],[event['right']],[event['left'],event['right']]),observable=bases>=.5*max(1,hi-lo)))
    supporting={r['sample'] for r in evidence if r['observable'] and r['consistent']}
    opposing={r['sample'] for r in evidence if r['observable'] and not r['consistent']}
    resolution=event.get('localization_resolution_bp',0)
    nearby=option['lo']<=event['hi']+resolution and option['hi']>=event['lo']-resolution
    localized=event.get('localization')=='Independent-individual regional chromosome transition'
    if len(opposing)>=2:status='contradicted';reason='Two independent individuals show a different intervening chromosome order.'
    elif opposing:status='unresolved';reason='An intervening chromosome assignment conflicts; independent confirmation is incomplete.'
    elif len(supporting)>=2:status='confirmed';reason='Two independent individuals support a consistent intervening chromosome order.'
    elif nearby and localized:status='confirmed';reason='Gap overlaps the cohort-localized transition within its recorded measurement resolution.'
    else:status='unresolved';reason='Intervening chromosome assignment is insufficient; no contradictory order is established.'
    return dict(status=status,reason=reason,resolution_bp=resolution,within_resolution=nearby,
                supporting_individuals=sorted(supporting),opposing_individuals=sorted(opposing),measurements=evidence)


def flank_support(path,lo,hi,length,width=10000):
    """Measure adjacent neighborhoods, retaining sparse endpoint counts separately."""
    from statistics import median
    summaries=[]
    for side,start,end in [('left',max(0,lo-width),max(0,lo-1000)),('right',min(length,hi+1000),min(length,hi+width))]:
        grid=support_grid(path,start,end)
        counts=[p['molecules'] for p in grid['probes']]
        summaries.append(dict(side=side,start=start,end=end,minimum=min(counts),median=median(counts),maximum=max(counts),
            supported_fraction=grid['supported_fraction'],probes=grid['probes'],
            observable=len(counts)>=5 and median(counts)>=2 and grid['supported_fraction']>=.8))
    return summaries


def observed_flanks(option):
    return len(option.get('flank_support',[]))==2 and all(f['observable'] for f in option['flank_support'])


def cut_prerequisites(option):
    """State every unavailable or opposing prerequisite; never select a cut."""
    failures=[];peers=option['peers']
    if len(peers['separate'])<2:failures.append('Fewer than two independent individuals pass both local chromosome flank tests.')
    if peers['continuous']:failures.append('Other individuals support the same chromosome on both sides.')
    if peers['conflicting']:failures.append('Local chromosome assignments conflict between individuals or haplotypes.')
    if 'spanning' not in option:failures.append('Gap-spanning HiFi support was not measured.')
    elif option['spanning']>=2:failures.append('At least two HiFi molecules span this gap.')
    if not observed_flanks(option):failures.append('Both adjacent HiFi neighborhoods have not established sufficient read observability.')
    if option.get('native_continuity'):failures.append('A native contig path continues across the gap.')
    if option.get('association',{}).get('status')!='confirmed':failures.append('Association between this gap and the chromosome transition is unresolved or contradicted.')
    contacts=option.get('contacts',[])
    if not any(c['calibrated'] and c['loss'] for c in contacts):failures.append('No available Hi-C library establishes contact loss against suitable controls.')
    if any(c['calibrated'] and c['loss'] for c in contacts) and any(c['calibrated'] and not c['loss'] for c in contacts):
        failures.append('Measured Hi-C libraries disagree about contact loss; review the library disagreement.')
    return failures


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
    weak=[];weak_points=[];start=None
    for i,(point,n) in enumerate(zip(points,counts)):
        if n<2 and start is None:start=point
        if start is not None and (n>=2 or i==len(points)-1):
            end=point if n>=2 else min(hi,point+step)
            first=points.index(start)
            if end>start:
                weak.append(dict(start=start,end=end,minimum=min(counts[first:i+1]),boundary_censored=(first==0 or i==len(points)-1)))
            else:weak_points.append(dict(position=start,molecules=counts[first],boundary_censored=True))
            start=None
    fraction=sum(n>=2 for n in counts)/len(counts);longest=max((v['end']-v['start'] for v in weak),default=0)
    state='continuous' if min(counts)>=2 else 'mostly_continuous' if fraction>=.95 and longest<=3000 else 'patchy'
    return dict(mapq=mapq,anchor_bp=anchor,step_bp=step,minimum=min(counts),supported_fraction=fraction,
                longest_weak_run_bp=longest,state=state,weak_runs=weak,weak_points=weak_points,probes=[dict(position=p,molecules=n) for p,n in zip(points,counts)])


def contact_result(focal,controls):
    def ratio(row):
        denominator=math.sqrt(row.get('left_within',0)*row.get('right_within',0))
        return row.get('cross',0)/denominator if denominator else None
    r=ratio(focal)
    audit=[];matched=[]
    for control in controls:
        reason=control.get('context_rejection','')
        if not reason:
            if min(control.get('left_within',0),control.get('right_within',0))<100:reason='insufficient_within_side_contacts'
            elif ratio(control) is None:reason='unobservable_contact_ratio'
            elif not all(.1<=control.get(k,0)/max(1,focal.get(k,0))<=10 for k in ('left_within','right_within')):reason='flank_depth_mismatch'
        audit.append(dict(control,ratio=ratio(control),accepted=not reason,rejection=reason or 'accepted'))
        if not reason:matched.append(control)
    values=sorted(ratio(c) for c in matched)
    floor=values[int(.1*(len(values)-1))] if values else None
    upper=(focal.get('cross',0)+3)/math.sqrt(focal['left_within']*focal['right_within']) if r is not None else None
    focal_observable=min(focal.get('left_within',0),focal.get('right_within',0))>=100
    calibrated=focal_observable and len(matched)>=5
    failure=[]
    if not focal_observable:failure.append('fewer_than_100_within_side_pairs')
    if len(matched)<5:failure.append('fewer_than_5_matched_controls')
    return dict(raw_counts=dict(focal),ratio=r,matched_controls=len(matched),comparison_floor=floor,
                upper_ratio=upper,calibrated=calibrated,loss=bool(calibrated and floor>0 and upper<.25*floor),control_ratios=values,
                control_ids=[c['id'] for c in matched if 'id' in c],control_audit=audit,calibration_failure=failure,
                tested_controls=len(controls),minimum_within_pairs=100,minimum_matched_controls=5,depth_matching_factor=10)



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


def matched_control_sites(intervals,gaps,lengths,flank=250000,count=96):
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
        # Represent multiple scaffolds before spending the pool on nearby gaps.
        by_scaffold=defaultdict(list)
        for site in available:by_scaffold[site['scaffold']].append(site)
        chosen=[];offsets={sc:0 for sc in by_scaffold}
        while len(chosen)<count:
            added=False
            for scaffold,sites in by_scaffold.items():
                while offsets[scaffold]<len(sites):
                    site=sites[offsets[scaffold]];offsets[scaffold]+=1
                    if any(site['scaffold']==v['scaffold'] and abs(site['lo']-v['lo'])<2*flank for v in chosen):continue
                    chosen.append(site);added=True;break
                if len(chosen)>=count:break
            if not added:break
        links[key]=[]
        for site in chosen:
            role='gap_control' if gap_target else 'continuous_control'
            cid=role+'_'+hashlib.sha256(('%s:%d:%d'%(site['scaffold'],site['lo'],site['hi'])).encode()).hexdigest()[:12]
            controls[cid]=dict(site,id=cid,role=role);links[key].append(cid)
    return controls,links


def recommend(event):
    peers=event['peers'];grid=event.get('hifi',{});options=event.get('cut_options',[])
    event.pop('recommended_cut',None)
    for option in options:option['unmet_cut_prerequisites']=cut_prerequisites(option)
    supported=grid.get('state') in ('continuous','mostly_continuous')
    strict=event.get('hifi_strict',{})
    contacts=event.get('contacts',[])
    losses=[c['library'] for c in contacts if c['calibrated'] and c['loss']]
    nonlosses=[c['library'] for c in contacts if c['calibrated'] and not c['loss']]
    event['evidence_conflicts']=['HiFi continuity versus measured Hi-C contact loss in '+library for library in losses] if supported else []
    if losses and nonlosses:event['evidence_conflicts'].append('Hi-C libraries disagree: contact loss in '+', '.join(losses)+'; no demonstrated loss in '+', '.join(nonlosses))
    if supported and strict.get('state') not in ('continuous','mostly_continuous'):
        event['evidence_conflicts'].append('MAPQ-20 continuity is not confirmed by the stricter MAPQ-30 read check')
    if event['evidence_conflicts']:
        return 'REVIEW','Keep intact pending review: '+ '; '.join(event['evidence_conflicts'])+'.'
    if len(peers['continuous'])>=2 and not peers['separate'] and not peers['conflicting']:
        return 'RETAIN','Other individuals support one chromosome on both sides.'
    eligible=[c for c in options if not c['unmet_cut_prerequisites']]
    if len(eligible)==1:
        event['recommended_cut']=eligible[0]
        return 'CUT','Independent chromosome separation, observable HiFi flanks without gap-spanning support, and contact loss in the available Hi-C evidence support this exact gap cut.'
    if any(c.get('spanning',0)>=2 and c['lo']<=event['lo']<=event['hi']<=c['hi'] for c in options):
        return 'RETAIN','At least two HiFi molecules directly span the verified gap containing this transition.'
    if any(len(c['peers']['separate'])>=2 and not c['peers']['continuous'] and not c['peers']['conflicting'] and
           c.get('spanning')==0 and observed_flanks(c) and c.get('association',{}).get('status')=='confirmed' and not c.get('native_continuity') for c in options):
        return 'SUSPECT','A nearby chromosome-consistent gap has observable read support on both sides but no spanning HiFi molecule; its contact evidence needs review.'
    if supported:
        return 'RETAIN','Local HiFi molecules support continuity across the transition'+(' with short weak pockets.' if grid['state']=='mostly_continuous' else '.')
    if peers['conflicting'] or peers['continuous']:return 'REVIEW','Chromosome comparisons conflict; leave intact while reviewing the conflicting evidence.'
    if len(peers['separate'])>=2 or len(event.get('independent_individuals',[]))>=2:
        if options:return 'SUSPECT','The chromosome transition is supported; the gap evidence does not yet meet the stated cut recommendation requirements.'
        return 'SUSPECT','The chromosome transition is supported; an exact cut site has not yet been established.'
    return 'REVIEW','This is a screening signal with insufficient independent confirmation.'
