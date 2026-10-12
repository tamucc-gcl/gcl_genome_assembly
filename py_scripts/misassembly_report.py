#!/usr/bin/env python3
"""One decision record and one readable report, backed by retained evidence."""
import argparse
import json
import shutil
from pathlib import Path
from misassembly_core import write_table,POLICY,table,cut_prerequisites

FIELDS=['selected','review_disposition','id','assembly','scaffold','chromosomes','event_id','event_assessment','event_reason','recommendation','option_advice','option_association','location_kind','location_rank','assay_probe_bp','cut_bp','review_start','review_end','action','gap_start','gap_end',
        'localization_status','reviewer','reason','coordinate_stage','assessment_sha256','decision_source','evidence_packet_id']


def bp(value):return format(int(value),',')


def mdtable(headers,rows):
    return ['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(str(v).replace('|','&#124;').replace('\n',' ') for v in row)+' |' for row in rows]


def advice(event):
    return {'CUT':'Recommend cutting at the specified gap','RETAIN':'Keep intact — local continuity supported','SUSPECT':'Suspect misjoin — review the evidence','REVIEW':'Review — keep intact pending a decision'}[event['recommendation']]



CONTROL_REASONS={
    'insufficient_independent_same_chromosome_context':'chromosome context incomplete',
    'hifi_mapping_unavailable':'HiFi unavailable',
    'continuous_control_without_hifi_continuity':'control read continuity unsupported',
    'insufficient_within_side_contacts':'too few within-side contacts',
    'unobservable_contact_ratio':'within-side support unavailable',
    'flank_depth_mismatch':'flank coverage does not match',
}


def control_summary(result):
    from collections import Counter
    reasons=Counter(c['rejection'] for c in result.get('control_audit',[]) if not c['accepted'])
    detail='; '.join(str(n)+' '+CONTROL_REASONS.get(reason,reason) for reason,n in sorted(reasons.items()))
    focal=result['raw_counts']
    if min(focal.get('left_within',0),focal.get('right_within',0))<100:
        detail=('candidate has fewer than 100 within-side pairs on a flank; '+detail).rstrip('; ')
    return str(result.get('tested_controls',len(result.get('control_audit',[]))))+' tested; '+str(result['matched_controls'])+' matched'+('; '+detail if detail else '')


def result_summary(data):
    active=sum(e['recommendation']!='RETAIN' for e in data['events'])
    retained=len(data['events'])-active
    if data['status']=='no_chromosome_scope':return 'Chromosome-scale classification unresolved' if data.get('scope_unresolved') else 'No chromosome-scale scaffolds identified'
    if data['status']!='assessed':return 'Not assessed — no informative cohort comparisons'
    if data['events']:return str(len(data['events']))+' chromosome transition(s) detected; '+str(active)+' need review; '+str(retained)+' provisionally keep intact'
    return 'No potential chimeras detected'

def option_advice(event,option):
    if event.get('recommended_cut',{}).get('id')==option['id']:return 'Recommended cut'
    association=option.get('association',{}).get('status','unresolved')
    if association!='confirmed':return 'Not suggested: nearby gap is not established as the transition location'
    if event['recommendation']=='RETAIN':return 'Not suggested: local evidence favors keeping intact'
    return 'Available at the transition; no cut recommendation'


def sister_context(event,data,assemblies):
    """Match detected chromosome pairs, without claiming identical junctions."""
    pair=sorted((event['left'],event['right']))
    return [(other,found) for other in assemblies if other['sample']==data['sample'] and other['assembly']!=data['assembly']
            for found in other['events'] if sorted((found['left'],found['right']))==pair]


def gap_summary(event):
    options=event.get('cut_options',[])
    if not options:return 'No available gap cut'
    return '; '.join(c['id']+': '+bp(c['cut_bp'])+' bp ('+
                    ('recommended' if event.get('recommended_cut',{}).get('id')==c['id'] else 'available; not recommended' if c.get('association',{}).get('status')=='confirmed' else 'nearby only; not associated')+')' for c in options)


def event_card(event,assembly,cohort=False,counterparts=()):
    md=['### '+(assembly+' · ' if cohort else '')+event['id']+' — '+event['left']+' → '+event['right'],'',
        '**'+advice(event)+'.** '+event['reason'],'',
        '**Location:** '+event['scaffold']+', transition range '+bp(event['lo'])+'–'+bp(event['hi'])+' bp.','',
        '**Evidence for breaking**','']
    if event.get('localization'):md+=['Boundary localization: '+event['localization']+'.','']
    peers=event['peers'];grid=event.get('hifi',{})
    broad=event.get('independent_individuals',[])
    if broad:md+=['- Alignments from '+str(len(broad))+' independent individuals identify different chromosome arms along this assembly: '+', '.join(broad)+'. This does not mean those individuals carry the fusion. The local flank tests below measure a separate, narrower question.']
    if event.get('localization_warning'):md+=['- '+event['localization_warning']]
    if event.get('screening_range'):
        bracket=event['screening_range'];md+=['- Original chromosome-arm bracket: '+bp(bracket['lo'])+'–'+bp(bracket['hi'])+' bp.']
    if peers['separate']:md+=['- '+str(len(peers['separate']))+' independent individuals place the flanks on different chromosomes: '+', '.join(peers['separate'])+'.']
    else:md+=['- Cohort-wide alignment detects a chromosome change; qualified local confirmation is incomplete.']
    if grid.get('state')=='patchy' and grid.get('weak_runs'):md+=['- The local read scan finds '+str(len(grid['weak_runs']))+' weak pocket(s); the longest is '+bp(grid['longest_weak_run_bp'])+' bp. Weak mapping alone does not establish a misjoin.']
    options=event.get('cut_options',[])
    unsupported_gaps=[c for c in options if c.get('association',{}).get('status')=='confirmed' and c.get('spanning')==0]
    md+=['','**Neighboring-sequence support**' if unsupported_gaps else '**Evidence for retaining**','']
    if grid.get('probes'):md+=['- '+format(100*grid['supported_fraction'],'.1f')+'% of sampled positions have at least two HiFi molecules crossing the position with 1 kb on each side. Minimum: '+str(grid['minimum'])+'; MAPQ ≥20. Small indels are tolerated.']
    else:md+=['- HiFi support was not measured.']
    if unsupported_gaps:md+=['- These nearby supported positions do not establish support across '+', '.join(c['id'] for c in unsupported_gaps)+': zero HiFi molecules span the actual gap.']
    if peers['continuous']:md+=['- Same-chromosome context in '+', '.join(peers['continuous'])+'.']
    if event.get('native_continuity'):md+=['- Both sides occur within one native hifiasm contig. This is assembly context, not independent read confirmation.']
    if event.get('evidence_conflicts'):md+=['','**Conflicting evidence:** '+ '; '.join(event['evidence_conflicts'])+'. Retention of the sequence pending review does not resolve that conflict.','']
    if peers['conflicting']:md+=['- Conflicting chromosome assignments in '+', '.join(peers['conflicting'])+'.']
    sisters=[t for t in event.get('peer_tests',[]) if t.get('sister')]
    if counterparts:
        md+=['','**Other haplotype — detected chromosome pair:** '+ '; '.join(
            '[ '+other['assembly']+' '+found['id']+' ]('+('assemblies/'+other['assembly']+'/report.md' if cohort else '../'+other['assembly']+'/report.md')+'#'+found['id'].lower()+'): '+found['left']+' → '+found['right']+' on '+found['scaffold'] for other,found in counterparts)+'. The same chromosome pair is detected; identical junctions have not been established. This is within-individual context, not an independent vote.','']
    if sisters:md+=['','**Other-haplotype local flank measurement:** '+ '; '.join(t['peer']+': '+str(t['left']['chrom'] or '?')+' → '+str(t['right']['chrom'] or '?')+' (qualified)' if t['qualified'] else t['peer']+': incomplete; no conclusion about the counterpart transition' for t in sisters)+'. Detailed assignments remain in the peer table.','']
    md+=['','**Cut options**','']
    if options:
        md+=mdtable(['Option','Gap / cut position','Individuals separating flanks','HiFi crossing / immediate left / right','Hi-C contact loss','Use'],[
            [c['id'],bp(c['lo'])+'–'+bp(c['hi'])+'; cut '+bp(c['cut_bp']),len(c['peers']['separate']),str(c.get('spanning','Not measured'))+' / '+' / '.join(map(str,c.get('flank_molecules',['?','?']))),
             ', '.join(t['library'] for t in c.get('contacts',[]) if t['loss']) or 'Not established',
             option_advice(event,c)] for c in options])
        md+=['','An available gap is not necessarily a suggested cut. Its association with the chromosome transition is assessed separately. Do not cut every option. All options begin unselected.','']
        for option in options:
            if option.get('association',{}).get('status')!='confirmed':
                distance=max(0,event['lo']-option['hi'],option['lo']-event['hi'])
                md+=['**'+option['id']+': nearby gap, not established as the transition location** ('+bp(distance)+' bp from the assessed range). It is not a suggested cut for this event.','']
            if option.get('origin'):md+=['**'+option['id']+' gap source:** '+('Current scaffolding AGP.' if option['origin']=='current_agp' else 'Literal N gap inherited inside a scaffold component.')+'','']
            if option.get('contacts'):
                md+=['**'+option['id']+' — contact evidence at this exact gap**','']+mdtable(['Library','Cross / within-left / within-right','Matched controls','Finding'],[
                    [t['library'],str(t['raw_counts']['cross'])+' / '+str(t['raw_counts']['left_within'])+' / '+str(t['raw_counts']['right_within']),control_summary(t),
                     'Supports cutting: contact loss' if t['loss'] else 'No demonstrated loss' if t['calibrated'] else 'Unable to judge contact loss'] for t in option['contacts']])+['']
            failures=option.get('unmet_cut_prerequisites',[])
            if failures:md+=['**Why this gap has no cut recommendation:**','']+['- '+reason for reason in failures]+['']
            association=option.get('association',{})
            if association:md+=['**Gap association: '+association['status']+'.** '+association['reason'],'']
            if option.get('flank_support'):
                md+=['**Adjacent HiFi support (10 kb neighborhoods):**','']+mdtable(['Side','Minimum / median / maximum molecules','Supported probes','Observability'],[[f['side'],str(f['minimum'])+' / '+str(f['median'])+' / '+str(f['maximum']),format(100*f['supported_fraction'],'.0f')+'%','Observed' if f['observable'] else 'Insufficient'] for f in option['flank_support']])+['','Immediate left/right counts in the option table are single probes; neighborhood counts assess whether a sparse endpoint represents a wider lack of reads.','']
            if option.get('native_continuity'):md+=['Both sides of '+option['id']+' map within native contig '+option.get('native_left','?')+'; this argues against treating it as a simple scaffolding join.','']
    else:md+=['No scaffold gap in or within 250 kb of this transition is available. An internal sequence cut needs an exact, reviewed coordinate; the range midpoint is not a suggested cut.','']
    if event.get('contacts'):
        md+=['**Hi-C measurements at the transition**','']+mdtable(['Library','Cross / within-left / within-right pairs','Matched controls','Interpretation'],[
            [t['library'],str(t['raw_counts']['cross'])+' / '+str(t['raw_counts']['left_within'])+' / '+str(t['raw_counts']['right_within']),control_summary(t),
             'Contact loss' if t['loss'] else 'No clear contact loss' if t['calibrated'] else 'Too few usable contacts or matched controls'] for t in event['contacts']])
        md+=['','Hi-C contributed to scaffolding and is corroborating evidence. A usable comparison needs at least 100 within-side pairs on each flank and five control boundaries whose flank counts each match within a factor of ten. Without these controls, low cross-boundary counts cannot be judged against an expected baseline: this is neither support for retaining nor a demonstrated contact loss. One library can provide this evidence; additional libraries are checked for disagreement. The control table lists measurements and exclusions. Measured zero counts are shown as 0; unavailable measurements are not zeros.','']
    if event.get('hifi_strict',{}).get('probes'):
        strict=event['hifi_strict']
        finding='Local continuity is supported.' if strict['state'] in ('continuous','mostly_continuous') else 'Local continuity is not established at this confidence level; review possible mapping ambiguity.'
        md+=['**Stricter read check:** '+format(100*strict['supported_fraction'],'.1f')+'% of positions supported at MAPQ ≥30; minimum '+str(strict['minimum'])+' molecules. '+finding,'']
        if strict.get('weak_runs'):
            md+=mdtable(['Weak interval in stricter check (bp)','Minimum spanning molecules'],[[bp(w['start'])+'–'+bp(w['end']),w['minimum']] for w in strict['weak_runs']])+['','These intervals identify weak read support, not selected cuts.','']
    if event.get('hifi_all_primary',{}).get('probes'):
        md+=['**Mapping ambiguity check:** '+format(100*event['hifi_all_primary']['supported_fraction'],'.1f')+'% of positions are crossed when low-confidence primary placements are also included. This helps distinguish absent reads from uncertain placement; these reads do not authorize cutting or establish continuity.','']
    if grid.get('weak_points'):
        md+=['**Unsupported positions at the assay edge:** '+', '.join(bp(p['position'])+' bp ('+str(p['molecules'])+' molecules)' for p in grid['weak_points'])+'. Their extent is not bounded by this assay; they are not zero-width break intervals.','']
    if grid.get('weak_runs') and event['recommendation']!='RETAIN':
        md+=['**Read support localization**','']+mdtable(['Weak interval (bp)','Minimum spanning molecules'],[[bp(w['start'])+'–'+bp(w['end']),w['minimum']] for w in grid['weak_runs']])+['','These are mapping support intervals, not selected cuts. An internal cut still requires a reviewed coordinate within the evidence interval.','']
    alternatives=event.get('alternatives',[])
    if alternatives:
        md+=['**Ranked alternative locations to review**','',
             'These are measured support-loss regions, not proposed cut coordinates. Rank prioritizes review; it is not a probability of a misjoin. The probe is the exact site used for bridging-read and Hi-C assays.','']
        md+=mdtable(['Rank / location','Support-loss range, bp','HiFi minimum: MAPQ20 / MAPQ30 / all primary','Separating / same-chromosome individuals','Assessment'],[
            [str(t['rank'])+' / '+t['id'],bp(t['lo'])+'–'+bp(t['hi']),str(t['hifi']['minimum'])+' / '+str(t['hifi_strict']['minimum'])+' / '+str(t['hifi_all_primary']['minimum']),str(len(t['peers']['separate']))+' / '+str(len(t['peers']['continuous'])),t['assessment']] for t in alternatives])+['']
        for t in alternatives:
            md+=['**'+t['id']+'**, assay probe '+bp(t['probe_bp'])+' bp: '+t['reason']+'.','',
                 'Association with the original chromosome transition: **'+t['association']['status']+'**. '+t['association']['reason'],'']
            if t.get('bridge_counts') is not None:
                counts=t['bridge_counts']
                md+=['Actual reads bridging 1 kb on each side of this probe: **'+str(counts['20'])+' MAPQ ≥20; '+str(counts['30'])+' MAPQ ≥30; '+str(counts['0'])+' including low-confidence primary placements**. '+str(t.get('bridges_with_competing_placements',0))+' have competing placements recorded in this window.','']
            if t.get('contacts'):
                md+=mdtable(['Library','Cross / within-left / within-right','Matched controls','Finding'],[
                    [c['library'],str(c['raw_counts']['cross'])+' / '+str(c['raw_counts']['left_within'])+' / '+str(c['raw_counts']['right_within']),c['matched_controls'],'Contact loss' if c['loss'] else 'No demonstrated loss' if c['calibrated'] else 'Uncalibrated'] for c in t['contacts']])+['']
        if event.get('unassessed_alternative_count'):md+=['Additional support-loss regions not individually assayed: '+str(event['unassessed_alternative_count'])+'. Review remains incomplete for these regions.','']
        md+=['Actual bridging reads and competing placements observed within the assessment window are retained in the bridging-read table. MAPQ is a placement-confidence filter; low-confidence bridges do not establish unique continuity. Placements outside this window have not been exhaustively enumerated.','']
    if cohort:md+=['[Evidence plots and peer-by-peer comparisons](assemblies/'+assembly+'/report.md#'+event['id'].lower()+')','']
    return md


def decisions(data):
    rows=[]
    for event in data['events']:
        targets=(event['cut_options'] or [dict(id=event['id'],lo=event['lo'],hi=event['hi'])])+event.get('alternatives',[])
        for target in targets:
            alternative=target.get('kind')=='regional_read_loss'
            gap='cut_bp' in target and not alternative
            rows.append(dict(selected='NO',review_disposition='PENDING',id=target['id'],assembly=data['assembly'],scaffold=event['scaffold'],
                chromosomes=event['left']+' → '+event['right'],event_id=event['id'],event_assessment=event['recommendation'],event_reason=event['reason'],recommendation=advice(event),
                option_advice=target['assessment']+'; assay probe is not a cut' if alternative else option_advice(event,target) if gap else 'No available gap cut; an internal cut requires a reviewed coordinate',option_association=target.get('association',{}).get('status','unresolved') if gap or alternative else 'no_gap',cut_bp=target.get('cut_bp',''),review_start=target['lo'] if alternative else event['lo'],review_end=target['hi'] if alternative else event['hi'],
                location_kind='alternative_region' if alternative else 'gap_option' if gap else 'original_transition',location_rank=target.get('rank',''),assay_probe_bp=target.get('probe_bp',''),
                action='UNJOIN_UNSUPPORTED' if gap else 'UNRESOLVED',gap_start=target['lo'] if gap else '',gap_end=target['hi'] if gap else '',localization_status='available_gap' if gap else 'unlocalized',
                reviewer='',reason='',coordinate_stage='pre_finishing',assessment_sha256=data['assessment_sha256'],decision_source='review',evidence_packet_id=data['assembly']))
    return rows


INSTRUCTIONS='''# Review and apply cuts

Copy `review.tsv` outside generated output before editing. Recommendations are advice; all selections start at NO.

`event_id`, `event_assessment`, `event_reason` and `recommendation` describe the detected transition. `option_advice` and `option_association` describe the particular gap cut, separately. A nearby unassociated gap is not a suggested cut just because it is listed.

| Decision | review_disposition | selected | Required |
| --- | --- | --- | --- |
| Cut | CUT | YES | reviewer, reason, exact cut_bp and cutting action |
| Keep | RETAIN | NO | reviewer and reason |
| Decide later | DEFER | NO | reviewer and reason |
| Unreviewed | PENDING | NO | Nothing |

For an approved gap cut, keep its supplied cut_bp, gap_start, gap_end and action=UNJOIN_UNSUPPORTED. Select one option per intended cut, not every nearby option. A retained event's gap options remain unselected.

For an internal cut, supply an exact cut_bp, action=BREAK_PROBABLE_MISJOIN, localization_status=localized and blank gap fields. Do not use a broad range midpoint without evidence.

Rows with location_kind=alternative_region describe ranked support-loss regions. assay_probe_bp is where evidence was measured; it is not a cut recommendation. These rows deliberately have blank cut_bp. Select a cut only after reviewing the region and supplying an exact coordinate.

To add an undetected cut, append a row with a new unique id. Copy assembly, assessment_sha256 and coordinate_stage from assembly_registry.tsv. Enter the original scaffold, exact coordinates, reviewer, reason, decision_source=review and evidence_packet_id=assembly. Set CUT and YES. Multiple cuts use separate rows, all in original coordinates.

Coordinates are zero-based on the original pre-finishing scaffold. cut_bp divides [0,cut_bp) and [cut_bp,length). All sequence, including gap Ns, is preserved.

Run the same full pipeline with `--chimera_break /absolute/path/reviewed.tsv -resume`. The cutter validates checksums, scaffold identity, literal gaps, conflicts and minimum piece lengths. It emits corrected assemblies, a cut audit, coordinate map and exact reconstruction verification. Chromosome reassignment and finishing use the corrected sequences. Automated cutting remains deferred.
'''


def main():
    p=argparse.ArgumentParser();p.add_argument('--packet',action='append',required=True);p.add_argument('--actions');p.add_argument('--verification',action='append',default=[]);a=p.parse_args()
    packet_data=[(Path(value),json.loads((Path(value)/'evidence.json').read_text())) for value in sorted(a.packet)]
    all_data=[data for packet,data in packet_data]
    if len({d['assembly'] for d in all_data})!=len(all_data):raise ValueError('Duplicate assembly packet')
    rows=[];registry=[];manifest={};versions=[]
    applied=[json.loads(Path(v).read_text()) for v in a.verification]
    reviewed=table(a.actions) if a.actions else []
    for packet,data in packet_data:
        assembly=data['assembly']
        for event in data['events']:
            for option in event.get('cut_options',[]):option['unmet_cut_prerequisites']=cut_prerequisites(option)
        target=Path('assemblies')/assembly;target.parent.mkdir(exist_ok=True)
        shutil.copytree(packet,target);rows+=decisions(data)
        (target/'evidence.json').write_text(json.dumps(data,indent=2),encoding='utf-8')
        versions.extend(dict(process='MISASSEMBLY_ASSESS',tool=tool,version=version) for tool,version in data.get('tools',{}).items())
        mapping=data['mapping']
        if mapping.get('bam'):manifest[assembly]=dict(bam=mapping['bam'],index=mapping['index'],fasta=data.get('assessment_fasta'),provenance=str((packet/'evidence.json').resolve()))
        result=result_summary(data)
        registry.append(dict(assembly=assembly,assessment_sha256=data['assessment_sha256'],assessment_fasta=data.get('assessment_fasta',''),coordinate_stage='pre_finishing',assessment_status=data['status'],detected_transitions=len(data['events']),review_candidates=sum(e['recommendation']!='RETAIN' for e in data['events']),retained_boundaries=sum(e['recommendation']=='RETAIN' for e in data['events']),rejected_screens=sum(v.get('kind')=='rejected_peer_transition' for v in data.get('minor_signals',[])),result=result))
        md=['# '+assembly+' — misassembly review','','[All assemblies](../../README.md) · [Manual decisions](../../cut-instructions.md)','', '**'+result+'.**','']
        if data.get('assessment_fasta'):md+=['Coordinates and read mappings refer to this original assessment FASTA: `'+data['assessment_fasta']+'`.','']
        if data.get('screening_coverage'):
            coverage=data['screening_coverage'];screen=[]
            for scaffold in data['scope']:
                useful=[r for r in coverage if r['scaffold']==scaffold and r['eligible'] and r['sample']!=data['sample']]
                individuals={r['sample'] for r in useful if r['assigned_bp']>=data['minimum_arm_bp']}
                fractions=[r['assigned_fraction'] for r in useful if r['assigned_fraction'] is not None]
                screen.append([scaffold,len(individuals),format(100*max(fractions),'.1f')+'%' if fractions else 'Unavailable'])
            md+=['<details>','<summary>What was screened</summary>','']+mdtable(['Chromosome-scale scaffold','Individuals with chromosome assignments','Best assigned fraction'],screen)+['','An absence of detected transitions applies to the usable comparisons above. Unmapped sequence is not evidence of continuity.','</details>','']
        for event in sorted(data['events'],key=lambda e:e['recommendation']=='RETAIN'):
            md+=['<a name="'+event['id'].lower()+'"></a>','']+event_card(event,assembly,counterparts=sister_context(event,data,all_data))
            if (target/'plots'/(event['id']+'.png')).exists():md+=['![Chromosome identity, HiFi support, Hi-C support and telomeres](plots/'+event['id']+'.png)','','Blue and orange distinguish the focal chromosome pair consistently across haplotypes; other chromosomes are muted. Pale backgrounds show the independent-individual consensus assignment; ambiguous regions remain neutral. Yellow marks the original transition bracket. Purple regions and R labels show the six highest-ranked alternative assay locations; the table retains all assessed locations. Dashed red lines mark available gap options, not selected cuts. Spanning molecules are HiFi reads. The Hi-C curves show binned support per library, not calibrated cut recommendations. Telomere counts are in the bottom panel. Gaps in curves mean unavailable measurements.','']
            md+=['<details>','<summary>Peer chromosome assignments</summary>','']+mdtable(['Assembly','Left → right','Local assignment','Role'],[
                [t['peer'],str(t['left']['chrom'] or '?')+' → '+str(t['right']['chrom'] or '?'),'Confirmed' if t['qualified'] else 'Partial or ambiguous',
                 'Sister haplotype' if t['sample']==data['sample'] else 'Independent individual' if t['eligible'] else 'Context only'] for t in event['peer_tests']])+['','</details>','']
        rejected=[v for v in data.get('minor_signals',[]) if v.get('kind')=='rejected_peer_transition']
        if rejected:
            md+=['<details>','<summary>Rejected chromosome-fusion screens ('+str(len(rejected))+')</summary>','']+mdtable(
                ['Scaffold / interval','Screened chromosomes','Same-chromosome individuals','Reason'],
                [[v['scaffold']+': '+bp(v['lo'])+'–'+bp(v['hi']),v['left']+' → '+v['right'],', '.join(v.get('peers',{}).get('continuous',[])),v['reason']] for v in rejected])+[
                '', 'These rejected fusion hypotheses do not establish that every local adjacency is correct. Their exact measurements and original detector observations remain in evidence.json.','</details>','']
        scaffold_figures=sorted((target/'plots').glob('*.contacts.png'))
        if scaffold_figures:
            md+=['## Scaffold contact maps','']
            for image in scaffold_figures:md+=['!['+image.stem+'](plots/'+image.name+')','']
        records=[('Alternative locations','alternative_locations.tsv'),('Actual bridging reads','bridging_reads.tsv'),('Read support','read_support.tsv'),('Hi-C counts','contact_support.tsv'),('Hi-C controls and exclusions','control_support.tsv'),('Chromosome assignments','chromosome_blocks.tsv'),('Telomere counts','telomeres.tsv'),('Full measurements and provenance','evidence.json')]
        links=['['+label+']('+name+')' for label,name in records if (target/name).is_file()]
        md+=['## Evidence records','', ' · '.join(links)+'.', '',
             'Small chromosome-labelled islands and rejected single-peer fusion hypotheses are retained in the machine-readable screening audit. They are not counted as candidates needing review. Missing homology is distinguished from evidence for a different placement.']
        (target/'report.md').write_text('\n'.join(md)+'\n',encoding='utf-8')
    if reviewed:
        sources={d['assembly']:d['assessment_sha256'] for d in all_data}
        identities=set()
        for row in reviewed:
            identity=(row.get('assembly'),row.get('id'))
            if identity in identities:raise ValueError('Duplicate manual review identity')
            identities.add(identity)
            if row.get('assembly') not in sources or row.get('assessment_sha256')!=sources[row['assembly']] or row.get('coordinate_stage')!='pre_finishing':raise ValueError('Manual review refers to a different source assembly')
            disposition=row.get('review_disposition')
            if row.get('selected') not in ('YES','NO') or disposition not in ('PENDING','CUT','RETAIN','DEFER') or (row['selected']=='YES')!=(disposition=='CUT'):raise ValueError('Manual selection and disposition disagree')
            if disposition!='PENDING' and (not row.get('reviewer','').strip() or not row.get('reason','').strip()):raise ValueError('Completed manual review needs reviewer and reason')
        by_key={(r['assembly'],r['id']):r for r in reviewed}
        rows=[by_key.pop((r['assembly'],r['id']),r) for r in rows]+list(by_key.values())
    write_table('review.tsv',rows,FIELDS);write_table('assembly_registry.tsv',registry,list(registry[0]))
    Path('mapping_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    write_table('versions.tsv',versions,['process','tool','version'])
    all_events=[(d,e) for d in all_data for e in d['events']]
    events=sorted(all_events,key=lambda pair:pair[1]['recommendation']=='RETAIN')
    counts={k:sum(e['recommendation']==k for _,e in all_events) for k in ('CUT','RETAIN','SUSPECT','REVIEW')}
    md=['# Chimera and misassembly review','',str(len(all_data))+' assemblies included; '+str(sum(d['status']=='assessed' for d in all_data))+' have usable cohort comparisons. **'+str(len(events))+' chromosome transitions detected: '+str(counts['CUT'])+' recommended cuts; '+str(counts['SUSPECT'])+' suspect misjoins; '+str(counts['REVIEW'])+' conflicting or incomplete assessments; '+str(counts['RETAIN'])+' provisionally keep intact.**','',
        'Detected chromosome transitions remain listed even when local evidence favors keeping the sequence intact. Local continuity does not prove biological chromosome fusion. Cut options begin unselected. [Editable decisions](review.tsv) · [Cut instructions](cut-instructions.md)','',
        '## Assembly overview','']+mdtable(['Assembly','Result','Evidence'],[[r['assembly'],r['result'],'[Report](assemblies/'+r['assembly']+'/report.md)'] for r in registry])
    if reviewed:
        md+=['','## Submitted manual decisions','',str(sum(r.get('selected')=='YES' for r in reviewed))+' cuts selected; '+str(sum(len(v['selected_actions']) for v in applied))+' actions verified by exact sequence reconstruction. These evidence measurements describe the original pre-finishing assemblies. Corrected assemblies are reassigned and finished downstream.','']
    if events:
        md+=['','## Detected chromosome transitions and suggested actions','']+mdtable(['Assembly / event','Chromosomes','Suggested action','Transition range','Gap cut options'],[
            ['['+d['assembly']+' '+e['id']+'](assemblies/'+d['assembly']+'/report.md#'+e['id'].lower()+')',e['left']+' → '+e['right'],advice(e),
             e['scaffold']+': '+bp(e['lo'])+'–'+bp(e['hi'])+' bp',gap_summary(e)] for d,e in events])
        md+=['','## Evidence for each decision','']
        for d,e in events:md+=event_card(e,d['assembly'],cohort=True,counterparts=sister_context(e,d,all_data))
    md+=['## Scope and reproducibility','',
         'Discovery compares all chromosome-scale scaffolds with the other assemblies, including screening the selected naming reference. The assigned-sequence arm threshold is recorded in each evidence packet; short islands remain in the audit. Individuals, not haplotypes, supply independent votes.', '',
         'HiFi continuity uses a grid of 1 kb probes with 1 kb flanks at MAPQ 20 and 30. Indels up to 50 bp are tolerated; larger indels and skipped regions interrupt support. Known scaffold N gaps are assessed separately. Retention can be supported by continuous or ≥95% supported regions with weak pockets no longer than 3 kb. These are explicit review heuristics, not proof of evolutionary fusion.', '',
         '[Assembly registry](assembly_registry.tsv) · [Retained mapping manifest](mapping_manifest.json). Policy: '+POLICY+'.']
    Path('README.md').write_text('\n'.join(md)+'\n',encoding='utf-8');Path('cut-instructions.md').write_text(INSTRUCTIONS,encoding='utf-8')
    Path('coordinate_status.md').write_text('# Misassembly coordinates\n\nAll discovery, evidence and manual actions refer to source-bound pre-finishing scaffolds. Cuts are verified by exact sequence reconstruction.\n',encoding='utf-8')


if __name__=='__main__':main()
