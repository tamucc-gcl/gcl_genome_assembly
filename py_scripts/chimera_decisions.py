"""Human review advice derived from evidence; never authorizes a sequence change."""
from collections import defaultdict


def bp(value):
    return format(int(value), ',')


def pair(row):
    return tuple((row.get('detected_transition') or row.get('chromosome_context','')).split(' → '))


def peer_counts(row, measurement, sample):
    expected=pair(row);qualified=defaultdict(set);partial=defaultdict(set)
    for track in measurement.get('chromosome_tracks',[]):
        if not track.get('auto_evidence') or track['sample']==sample:continue
        identity=(track['left'].get('chrom'),track['right'].get('chrom'))
        if not all(identity):continue
        target=partial if track['relationship']=='uninformative' else qualified
        target[track['sample']].add(identity)
    matches=[];same=[];conflicts=[]
    for individual,identities in qualified.items():
        if len(identities)!=1:conflicts.append(individual)
        elif expected in identities:matches.append(individual)
        elif next(iter(identities))[0]==next(iter(identities))[1]:same.append(individual)
        else:conflicts.append(individual)
    weaker=[s for s,identities in partial.items() if expected in identities and s not in qualified]
    weak_same=sorted(s for s,identities in partial.items() if any(a==b for a,b in identities))
    weak_other=sorted(s for s,identities in partial.items() if any(a!=b and (a,b)!=expected for a,b in identities))
    return dict(matches=sorted(matches),same=sorted(same),conflicts=sorted(conflicts),partial=sorted(weaker),weak_same=weak_same,weak_other=weak_other)


def advise(row,m,sample):
    """Explicit advice is distinct from the user's CUT/RETAIN disposition."""
    counts=peer_counts(row,m,sample)
    retained=row['assessment_status']=='retention_supported'
    positive_reads=(m.get('hifi_informative') and m.get('hifi_spanning_molecules',0)>=2) or m.get('local_path_support')=='supported_grid'
    if retained:
        code='KEEP_REJECTED_SIGNAL';advice='Keep — fusion signal rejected'
        reason='Local comparisons place both sides on the same chromosome.'
    elif counts['same'] or counts['conflicts'] or positive_reads:
        code='KEEP_CONFLICTING_EVIDENCE';advice='Keep for now — evidence conflicts'
        reason='Evidence for continuity or conflicting chromosome assignments must be resolved before cutting.'
    elif m.get('verified_gap') and len(counts['matches'])>=2 and not m.get('native_continuity') and row.get('cut_bp')!='':
        code='CONSIDER_GAP_CUT';advice='Consider cutting at this gap'
        reason='Multiple independent individuals support the detected chromosome pair, and an exact scaffolding gap is available.'
    elif not m.get('verified_gap'):
        code='KEEP_PENDING_LOCALIZATION';advice='Keep for now — locate the break first'
        reason='The chromosome transition is suspicious, but the evidence does not identify an exact cut.'
    else:
        code='KEEP_INSUFFICIENT_EVIDENCE';advice='Keep for now — cut evidence is insufficient'
        reason='An available gap alone does not establish a misjoin.'
    return dict(suggested_action=code,recommendation=advice,recommendation_reason=reason),counts


def focus(group,measurements,sample):
    ordered=sorted(group,key=lambda r:(advise(r,measurements.get(r['source_candidate_id'],{}),sample)[0]['suggested_action']!='CONSIDER_GAP_CUT',r.get('cut_bp')=='',r['id']))
    return ordered[0]


def location(row):
    if row.get('cut_bp')!='':return 'Gap cut: '+bp(row['cut_bp'])+' bp'
    if row.get('review_start') not in ('',None) and row.get('review_end') not in ('',None):
        return bp(row['review_start'])+'–'+bp(row['review_end'])+' bp; exact cut unknown'
    return 'Location unavailable'


def assessment_label(registry,groups):
    if groups:return str(len(groups))+' potential chimera'+('s' if len(groups)!=1 else '')
    if registry['assessment_status']=='unavailable':return 'Not assessed — reference alignment unavailable'
    if registry['assessment_status']=='chromosome scope unresolved':return 'No chromosome-scale scaffolds identified'
    if registry['assessment_status']=='assessed':return 'No potential chimeras detected'
    return 'Not assessed — inputs incomplete'


def evidence(row,m,sample):
    counts=peer_counts(row,m,sample);support=[];against=[];limits=[]
    if counts['matches']:
        support.append(str(len(counts['matches']))+' independent individual'+('s place' if len(counts['matches'])!=1 else ' places')+
                       ' the two sides on '+row['detected_transition']+'. ('+', '.join(counts['matches'])+'.)')
    if counts['partial']:
        support.append(str(len(counts['partial']))+(' additional' if counts['matches'] else '')+' individual'+('s show' if len(counts['partial'])!=1 else ' shows')+
                       ' the same chromosome pair, with too little local alignment for confirmation. ('+', '.join(counts['partial'])+'.)')
    if not support:support.append('The reference alignment changes chromosome identity here; independent local confirmation is missing.')
    if counts['same']:against.append(str(len(counts['same']))+' independent individuals place both sides on one chromosome. ('+', '.join(counts['same'])+'.)')
    if counts['conflicts']:against.append('Chromosome assignments conflict in '+', '.join(counts['conflicts'])+'.')
    if counts['weak_same']:against.append('Lower-confidence alignments place both sides on one chromosome in '+', '.join(counts['weak_same'])+'. Those alignments do not pass the local confirmation test.')
    if counts['weak_other']:limits.append('Lower-confidence alignments suggest a different chromosome pair in '+', '.join(counts['weak_other'])+'; those assignments do not pass the local confirmation test.')
    if m.get('native_continuity'):
        against.append('Both sides are within one original hifiasm contig. This supports the assembler’s choice, but does not independently validate it.')
    elif m.get('verified_gap'):
        support.append('The gap separates two original contigs; cutting would undo a scaffolding join.')
    if m.get('hifi_informative') and m.get('hifi_spanning_molecules',0)>=2:
        against.append(str(m['hifi_spanning_molecules'])+' HiFi molecules span the tested interval with adequate flank coverage.')
    if m.get('local_path_support')=='supported_grid':against.append('HiFi reads support continuity at the sampled positions across the region.')
    if not against:against.append('No decisive evidence for keeping this join was measured. This is not proof that it is wrong.')
    if row.get('cut_bp')=='':limits.append('Only a transition range is available. Its midpoint is not a justified cut.')
    bridge=m.get('transition_bridge_assessment',{})
    if bridge and bridge.get('status')!='no_reversal_observed_adequate_coverage':
        length=bridge['end']-bridge['start'];best=max(bridge.get('coverage_bp',{}).values(),default=0)
        limits.append('The nearby gap has not been shown to be the chromosome-change point: at best %.1f%% of the %s bp intervening region is chromosome-assigned in an independent peer.'%(100*best/max(1,length),bp(length)))
    if not m.get('hifi_informative') and 'hifi_spanning_molecules' in m:
        limits.append('HiFi coverage at one or both immediate flanks is too low for the spanning-read test.')
    libraries=m.get('libraries',[])
    if libraries and not any(t.get('informative') for t in libraries):
        limits.append('The Hi-C comparisons lack enough matched controls and/or usable flank contacts to classify this join reliably.')
    sisters=[t for t in m.get('chromosome_tracks',[]) if t['sample']==sample]
    if sisters and all(t['relationship']=='uninformative' for t in sisters):
        limits.append('The sister-haplotype comparison could not resolve this local join. Matching chromosome combinations alone do not prove the same junction.')
    return support,against,limits


def decision_card(group,measurements,provenance,cohort=False):
    sample=provenance.get('sample',group[0]['assembly'].rsplit('_hap',1)[0]);row=focus(group,measurements,sample)
    m=measurements.get(row['source_candidate_id'],{});advice,_=advise(row,m,sample)
    ids=' / '.join(r['id'] for r in group)
    label=(row['assembly']+' · ' if cohort else '')+ids+' — '+row['detected_transition']
    md=['### '+label,'', '**Suggested action: '+advice['recommendation']+'.** '+advice['recommendation_reason'],
        '', '**Where:** '+row['scaffold']+'; '+location(row)+'.','']
    if m.get('verified_gap'):
        md+=['**Cut option:** use row '+row['id']+' at **'+bp(row['cut_bp'])+' bp**, at the end of the '+bp(int(row['gap_end'])-int(row['gap_start']))+' bp N gap. This is a physical cut option, not a proven biological breakpoint.','']
    support,against,limits=evidence(row,m,sample)
    for title,items in [('Evidence for cutting',support),('Evidence for keeping',against),('What remains uncertain',limits)]:
        if items:md+=['**'+title+'**','']+['- '+item for item in items]+['']
    if len(group)>1:
        md+=['**Related measurement:** '+', '.join(r['id']+' tests '+location(r) for r in group if r is not row)+'. These measurements describe one possible fusion; they are not instructions to make two cuts.','']
    raw=m.get('hifi_raw',{});libs=m.get('libraries',[])
    if raw:
        md+=['**HiFi:** '+str(m.get('hifi_spanning_molecules','Not measured'))+' molecules span the assessed interval; '+str(raw.get('left_molecules','?'))+' / '+str(raw.get('right_molecules','?'))+' molecules cover the left / right flanks. '+
             ('The flank coverage test passes.' if m.get('hifi_informative') else 'The flank coverage test fails.')+
             (' The '+bp(int(row['review_end'])-int(row['review_start']))+' bp interval may exceed read lengths; zero spanning molecules alone does not support cutting.' if not m.get('verified_gap') else ' Zero spanning molecules with inadequate flanks is inconclusive.' if not m.get('hifi_informative') else ''),'']
    if libs:
        md+=['**Hi-C:** '+'; '.join(t['library']+': '+str(t['raw_counts'].get('cross',0))+' cross-join pairs, '+str(t['raw_counts'].get('left_within',0))+' / '+str(t['raw_counts'].get('right_within',0))+' within-flank pairs (left / right)' for t in libs)+'. '+
             ('Calibrated comparisons are available; see the assay table.' if any(t.get('informative') for t in libs) else 'These are measured counts, but the cut test is not calibrated.'),'']
    if cohort:md+=['[Plots and full evidence for this event]('+row['report_path']+')','']
    return md
