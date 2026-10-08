"""GitHub-readable evidence presentation; no cutting decisions are inferred."""
import json
import re


def peer_observability(t):
    if t['relationship']!='uninformative':
        return 'Both sides meet chromosome-assignment requirements'
    reasons=[]
    for side in ('left','right'):
        value=t[side]
        if value.get('qualified'):
            continue
        if not value.get('aligned_bp'):
            reason='no qualifying aligned bases after filtering; raw matches may exist'
        elif not value.get('chrom'):
            reason='aligned sequence has ambiguous or mixed chromosome assignment'
        else:
            reason='assignment to '+str(value['chrom'])+' has insufficient qualifying bases, dominance or informative bins'
        reasons.append(side+': '+reason)
    return '; '.join(reasons) or 'Chromosome relationship could not be classified'


def decision_context(m, sample):
    """Describe measured chromosome relationships; review priority never selects a cut."""
    groups = {}
    for t in m.get('chromosome_tracks', []):
        if t.get('auto_evidence') and t['sample'] != sample:
            groups.setdefault(t['sample'], set())
            if t['relationship'] in ('different_chromosomes', 'same_chromosome'):
                groups[t['sample']].add((t['relationship'], t['left']['chrom'], t['right']['chrom']))
    separate=[];same=[];unclear=[];conflict=[];pairs=set()
    for individual, values in sorted(groups.items()):
        if len(values)>1:
            conflict.append(individual)
        elif not values:
            unclear.append(individual)
        else:
            relationship,left,right=next(iter(values))
            text=individual+': '+str(left)+' → '+str(right)
            (separate if relationship=='different_chromosomes' else same).append(text)
            pairs.add(str(left)+' → '+str(right))
    partial=[]
    for t in m.get('chromosome_tracks',[]):
        if t.get('auto_evidence') and t['sample']!=sample and t['relationship']=='uninformative' and t['left'].get('chrom') and t['right'].get('chrom'):
            partial.append(t.get('peer',t['sample'])+': '+t['left']['chrom']+' → '+t['right']['chrom'])
    partial_individuals={t['sample'] for t in m.get('chromosome_tracks',[]) if t.get('auto_evidence') and t['sample']!=sample and t['relationship']=='uninformative' and t['left'].get('chrom') and t['right'].get('chrom') and t['left']['chrom']!=t['right']['chrom']}
    for_cut='Separate chromosomes: '+'; '.join(separate) if separate else ('Below-threshold separate-chromosome observations in '+str(len(partial_individuals))+' independent individuals; local confirmation incomplete' if partial_individuals else 'No qualified local evidence supporting a break')
    against='Same chromosome: '+'; '.join(same) if same else 'No opposing chromosome evidence observed'
    spans=m.get('hifi_spanning_molecules')
    if spans and m.get('hifi_informative'):
        against+='; '+str(spans)+' qualified immediate spanning molecules'
    if m.get('native_continuity'):
        against+='; continuity within one native primary contig (assembly context, not independent validation)'
    limits=[]
    if partial:limits.append('Below-threshold chromosome assignments (not qualified votes): '+'; '.join(partial))
    if unclear:limits.append('No informative two-sided chromosome assignment: '+', '.join(unclear))
    if conflict:limits.append('Discordant haplotype assignments within individual: '+', '.join(conflict))
    if not m.get('hifi_informative'):limits.append('Immediate HiFi assay not informative; zero spanning reads is inconclusive')
    if not m.get('hic_informative'):limits.append('Hi-C assay not calibrated/informative for a cut decision')
    if m.get('local_path_support')=='supported_grid':
        against+='; local HiFi continuity supported across the sampled transition'
    if same and not separate and not conflict:
        priority='Evidence favors retaining this sampled boundary'
    elif m.get('verified_gap') and len(separate)>=2 and not same and not conflict and not (spans and m.get('hifi_informative')) and not m.get('native_continuity') and m.get('local_path_support')!='supported_grid':
        priority='Prioritize gap-cut review'
    elif separate:
        priority='Investigate chromosome transition'+('; exact cut not localized' if not m.get('verified_gap') else '; mixed/limited evidence')
    elif partial_individuals:
        priority='Review chromosome transition; below-threshold observations need stronger local confirmation'
    else:
        priority='Insufficient evidence to propose a break'
    sisters=[t for t in m.get('chromosome_tracks',[]) if t['sample']==sample]
    sister_context='; '.join(t.get('peer',sample)+': '+str(t['left'].get('chrom') or '?')+' → '+str(t['right'].get('chrom') or '?')+'; '+str(t.get('placement_relationship','placement not assessed'))+'; '+('qualified chromosome assignments' if t['relationship']!='uninformative' else 'local chromosome confirmation incomplete') for t in sisters) or 'Sister haplotype comparison unavailable'
    partial_groups={}
    for t in m.get('chromosome_tracks',[]):
        if t.get('auto_evidence') and t['sample']!=sample and t['relationship']=='uninformative' and t['left'].get('chrom') and t['right'].get('chrom'):
            pair=t['left']['chrom']+' → '+t['right']['chrom']
            partial_groups.setdefault(pair,[]).append(t['sample'])
    partial_summary='; '.join(pair+': '+str(len(set(samples)))+' individuals / '+str(len(samples))+' haplotypes (below threshold)' for pair,samples in sorted(partial_groups.items()))
    return dict(chromosome_context='; '.join(sorted(pairs)) or 'Unresolved from qualified local alignments',
        evidence_for_cut=for_cut,evidence_against_cut=against,review_priority=priority,
        evidence_limits='; '.join(limits) or 'See local measurements and controls',
        local_confirmation=str(len(separate))+' independent individuals separate; '+str(len(same))+' same chromosome; '+str(len(conflict))+' discordant',
        partial_confirmation=partial_summary or 'No additional below-threshold chromosome pairs',sister_context=sister_context)


def assess_transitions(rows, measurements, provenance, candidates):
    """Persist event membership and bridge observations separately from cut selection."""
    intervals=provenance.get('intervals', {})
    for r in rows:
        c=candidates[r['source_candidate_id']]
        detected=(c.get('left_chrom'),c.get('right_chrom'))
        r['detected_transition']=' → '.join(detected) if all(x and x!='.' for x in detected) else ''
        if r['detected_transition'] and r['chromosome_context']=='Unresolved from qualified local alignments':
            r['chromosome_context']=r['detected_transition']
        r.update(assessment_status='review_required',transition_id=r['id'],preferred_candidate='',bridge_status='not_assessed',related_candidate='')
        if r['review_priority']=='Evidence favors retaining this sampled boundary':
            r['assessment_status']='retention_supported'
            r['action']='RETAIN'
            r['cut_bp']=''
            r['localization_status']='assessed_gap' if measurements.get(r['source_candidate_id'],{}).get('verified_gap') else 'unlocalized'
    for gap in rows:
        gm=measurements.get(gap['source_candidate_id'],{})
        if not gm.get('verified_gap') or gap['assessment_status']=='retention_supported':continue
        sources=set(gm.get('source_candidates',[]))
        linked=[r for r in rows if r is not gap and r['scaffold']==gap['scaffold'] and
                measurements.get(r['source_candidate_id'],{}).get('packet_interval_id') in sources]
        # Multiple source transitions remain separate: a nearby gap is not a unique localization.
        if len(linked)!=1:continue
        source=linked[0];pair=source['detected_transition'].split(' → ')
        if len(pair)!=2 or pair[0]==pair[1]:continue
        if gap['chromosome_context']!=' → '.join(pair):continue
        lo=min(int(source['review_end']),int(gap['review_start']))
        hi=max(int(source['review_end']),int(gap['review_start']))
        observations={};reversal=False
        for r in (source,gap):
            m=measurements.get(r['source_candidate_id'],{})
            interval=intervals.get(m.get('packet_interval_id'),{})
            if 'start' not in interval:continue
            for track in m.get('chromosome_tracks',[]):
                if not track.get('auto_evidence') or track['sample']==provenance.get('sample'):continue
                tiles=observations.setdefault(track['peer'],{})
                for b in track.get('bridge_segments',track.get('bins',[])):
                    a=interval['start']+b['lo'];z=interval['start']+b['hi']
                    if a<hi and z>lo and b.get('chrom'):
                        tiles[(a,z)]=b['chrom']
        sufficient=[];discordant_bp={};coverage_bp={}
        for peer,tiles in observations.items():
            ordered=sorted(tiles.items());seen_right=False;covered=[];bad_bp=0
            for (a,z),chrom in ordered:
                if chrom not in pair or (seen_right and chrom==pair[0]):bad_bp+=max(0,min(hi,z)-max(lo,a))
                if chrom==pair[1]:seen_right=True
                covered.append((max(lo,a),min(hi,z)))
            end=lo;bp=0
            for a,z in covered:
                bp+=max(0,z-max(a,end));end=max(end,z)
            discordant_bp[peer]=bad_bp
            coverage_bp[peer]=bp
            if bad_bp>=10000:reversal=True
            if ordered and bp>=.5*max(1,hi-lo) and bad_bp==0:sufficient.append(peer)
        status='reversal_or_other_chromosome_observed' if reversal else 'no_reversal_observed_adequate_coverage' if sufficient else 'insufficient_bridge_observability'
        intervening=[r['id'] for r in rows if r is not source and r is not gap and r['scaffold']==source['scaffold'] and r['review_start']!='' and lo<int(r['review_start'])<hi]
        if intervening:status='intervening_transition_requires_review'
        audit=dict(source=source['id'],gap=gap['id'],start=lo,end=hi,status=status,adequate_peers=sufficient,discordant_bp=discordant_bp,coverage_bp=coverage_bp,intervening_candidates=intervening,
                   coverage_fraction_required=.5,discordant_bp_threshold=10000,
                   peer_tiles={k:[dict(start=a,end=z,chrom=c) for (a,z),c in sorted(v.items())] for k,v in observations.items()})
        gm['transition_bridge_assessment']=audit
        source['related_candidate']=gap['id'];gap['related_candidate']=source['id']
        source['bridge_status']=gap['bridge_status']=status
        gap['detected_transition']=source['detected_transition']
        if status not in ('reversal_or_other_chromosome_observed','intervening_transition_requires_review'):
            source['transition_id']=gap['transition_id']=source['id']
        sibling_gaps=[r for r in rows if measurements.get(r['source_candidate_id'],{}).get('verified_gap') and
                      measurements.get(source['source_candidate_id'],{}).get('packet_interval_id') in
                      measurements.get(r['source_candidate_id'],{}).get('source_candidates',[]) and
                      r['scaffold']==source['scaffold'] and r['chromosome_context']==gap['chromosome_context']]
        if status=='no_reversal_observed_adequate_coverage' and len(sibling_gaps)==1:
            source['assessment_status']='supporting_measurement'
            source['transition_id']=gap['transition_id']=source['id']
            source['preferred_candidate']=gap['preferred_candidate']=gap['id']
            gap['detected_transition']=source['detected_transition']
    return rows



def render(*args, **kwargs):
    from chimera_report import render as write_report
    return write_report(*args, **kwargs)
