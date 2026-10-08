"""Chromosome-block and distance-matched evidence for repeat-obscured gap joins."""
from collections import defaultdict
from chimera_controls import flank_placement,informative_hifi,contact_ratio,sam_measure,project_query
import math
from chimera_tracks import chromosome_at

OFFSETS=(100000,250000,500000)

def chromosome_blocks(hits,interval,peer,gaps,anchor=50000):
    labels=peer.get('chromosome_labels',{});trials=[]
    for offset in OFFSETS:
        lo=interval['lo']-interval['start'];hi=interval['hi']-interval['start']
        bounds=(lo-offset-anchor,lo-offset,hi+offset,hi+offset+anchor)
        if bounds[0]<0 or bounds[3]>interval['end']-interval['start']:continue
        left=flank_placement(hits,bounds[0],bounds[1]);right=flank_placement(hits,bounds[2],bounds[3])
        lc=chromosome_at(labels,left[5],project_query(left,(bounds[0]+bounds[1])//2)) if left else None
        rc=chromosome_at(labels,right[5],project_query(right,(bounds[2]+bounds[3])//2)) if right else None
        trials.append(dict(offset_bp=offset,left_chrom=lc,right_chrom=rc,
            left_target=left[5] if left else None,right_target=right[5] if right else None,
            usable=bool(lc and rc),left_edge=interval['lo']-offset,right_edge=interval['hi']+offset))
    usable=[t for t in trials if t['usable']]
    pairs={(t['left_chrom'],t['right_chrom']) for t in usable}
    qualified=len(usable)>=2 and len(pairs)==1
    pair=next(iter(pairs)) if qualified else None
    localized=False
    if qualified and pair[0]!=pair[1]:
        nearest=min(usable,key=lambda t:t['offset_bp'])
        between=[g for g in gaps if g['scaffold']==interval['scaffold'] and nearest['left_edge']<=g['lo']<g['hi']<=nearest['right_edge']]
        localized=len(between)==1 and between[0]['lo']==interval['lo'] and between[0]['hi']==interval['hi']
    return dict(peer=peer['id'],sample=peer['sample'],auto_evidence=peer.get('auto_evidence') is True,
        chromosome_pair=list(pair) if pair else None,qualified=qualified,unique_gap_localization=localized,trials=trials)


def summarize_blocks(rows,sample):
    by_sample=defaultdict(list)
    for row in rows:
        if row['sample']!=sample and row['auto_evidence']:by_sample[row['sample']].append(row)
    agreed={}
    observed=set()
    for individual,values in by_sample.items():
        # Missing haplotypes are not contradictions; conflicting usable haplotypes are.
        qualified=[v for v in values if v['qualified']]
        pairs={tuple(v['chromosome_pair']) for v in qualified}
        observed.update(pairs)
        if len(pairs)==1:
            pair=next(iter(pairs))
            if pair[0]!=pair[1] and all(v['unique_gap_localization'] for v in qualified):agreed[individual]=pair
    counts=defaultdict(list)
    for individual,pair in agreed.items():counts[pair].append(individual)
    unanimous=len(counts)==1 and len(observed)==1
    supporting=next(iter(counts.values())) if unanimous else []
    pair=next(iter(counts)) if unanimous else None
    return dict(localized=unanimous and len(supporting)>=3,independent_individuals=supporting,
        chromosome_pair=list(pair) if pair else None,contradictory_pairs=len(observed)>1,peer_assays=rows)


def farther_hifi(path,interval):
    result={}
    for offset in OFFSETS:
        if interval['lo']-offset-1000<0 or interval['hi']+offset+1000>interval['end']:continue
        # These measure flank observability, never bridges across the enlarged interval.
        value=sam_measure(path,interval['lo']-offset,interval['hi']+offset)
        result[str(offset)]=dict(informative=informative_hifi(value),raw=value)
    return result


def farther_regions(intervals,lengths):
    result={}
    for key,interval in intervals.items():
        for offset in OFFSETS:
            lo,hi=interval['lo']-offset,interval['hi']+offset
            if lo<250000 or hi+250000>lengths[interval['scaffold']]:continue
            result[key+'@'+str(offset)]=dict(scaffold=interval['scaffold'],lo=lo,hi=hi,base=key,offset_bp=offset)
    return result


def farther_contacts(candidate,controls,counts,libraries,regions):
    trials=[]
    for offset in OFFSETS:
        key=candidate['key']+'@'+str(offset)
        if key not in regions:continue
        for lib in sorted(libraries):
            focal=counts[lib,key];matched=[]
            for c in controls:
                ck=c['key']+'@'+str(offset)
                if ck not in regions or not c.get('farther_hifi',{}).get(str(offset),{}).get('informative'):continue
                if c['role']=='continuous_control' and c['hifi']['spanning']<2:continue
                if c['role']=='gap_control' and c.get('peer_continuous_individuals',0)<3:continue
                value=counts[lib,ck]
                if min(value.get('left_within',0),value.get('right_within',0))<100:continue
                if not all(.5<=value.get(side,0)/max(1,focal.get(side,0))<=2 for side in ('left_within','right_within')):continue
                matched.append((c,contact_ratio(value)))
            populations={role:sum(c['role']==role for c,r in matched) for role in ('continuous_control','gap_control')}
            informative=(candidate.get('farther_hifi',{}).get(str(offset),{}).get('informative') is True and
                min(focal.get('left_within',0),focal.get('right_within',0))>=100 and all(n>=5 for n in populations.values()))
            minimum=min((r for c,r in matched),default=0)
            upper=(focal.get('cross',0)+3)/max(1,math.sqrt(focal.get('left_within',0)*focal.get('right_within',0)))
            trials.append(dict(offset_bp=offset,library=lib,informative=informative,control_populations=populations,
                matched_control_ids=[c['key'] for c,r in matched],minimum_control_ratio=minimum,upper_count_allowance_ratio=upper,
                support_loss=informative and minimum>0 and upper<.1*minimum,raw_counts=dict(focal)))
    informative_offsets=[o for o in OFFSETS if len([t for t in trials if t['offset_bp']==o])==len(libraries) and
        all(t['informative'] and t['support_loss'] for t in trials if t['offset_bp']==o)]
    contradiction=any(t['informative'] and not t['support_loss'] for t in trials)
    return dict(supported_offsets=informative_offsets,pass_all=len(libraries)>=2 and len(informative_offsets)>=2 and not contradiction,
        contradictory_informative_trial=contradiction,trials=trials)
