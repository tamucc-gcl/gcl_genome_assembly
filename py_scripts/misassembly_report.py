#!/usr/bin/env python3
"""One decision record and one readable report, backed by retained evidence."""
import argparse
import json
import shutil
from pathlib import Path
from misassembly_core import write_table,POLICY,table

FIELDS=['selected','review_disposition','id','assembly','scaffold','chromosomes','recommendation','cut_bp','review_start','review_end','action','gap_start','gap_end',
        'localization_status','reviewer','reason','coordinate_stage','assessment_sha256','decision_source','evidence_packet_id']


def bp(value):return format(int(value),',')


def mdtable(headers,rows):
    return ['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(str(v).replace('|','&#124;').replace('\n',' ') for v in row)+' |' for row in rows]


def advice(event):
    return {'CUT':'Recommend cutting at the specified gap','RETAIN':'Retain — evidence supports continuity','SUSPECT':'Suspect misjoin — review cut options / localization','REVIEW':'Review — evidence is incomplete or conflicting'}[event['recommendation']]


def event_card(event,assembly,cohort=False):
    md=['### '+(assembly+' · ' if cohort else '')+event['id']+' — '+event['left']+' → '+event['right'],'',
        '**'+advice(event)+'.** '+event['reason'],'',
        '**Location:** '+event['scaffold']+', transition range '+bp(event['lo'])+'–'+bp(event['hi'])+' bp.','',
        '**Evidence for breaking**','']
    peers=event['peers'];grid=event.get('hifi',{})
    if peers['separate']:md+=['- '+str(len(peers['separate']))+' independent individuals place the flanks on different chromosomes: '+', '.join(peers['separate'])+'.']
    else:md+=['- Cohort-wide alignment detects a chromosome change; qualified local confirmation is incomplete.']
    if grid.get('state')=='patchy' and grid.get('weak_runs'):md+=['- The local read scan finds '+str(len(grid['weak_runs']))+' weak pocket(s); the longest is '+bp(grid['longest_weak_run_bp'])+' bp. Weak mapping alone does not establish a misjoin.']
    md+=['','**Evidence for retaining**','']
    if grid.get('probes'):md+=['- '+format(100*grid['supported_fraction'],'.1f')+'% of sampled positions have at least two HiFi molecules crossing the position with 1 kb on each side. Minimum: '+str(grid['minimum'])+'; MAPQ ≥20. Small indels are tolerated.']
    else:md+=['- HiFi support was not measured.']
    if peers['continuous']:md+=['- Same-chromosome context in '+', '.join(peers['continuous'])+'.']
    if event.get('native_continuity'):md+=['- Both sides occur within one native hifiasm contig. This is assembly context, not independent read confirmation.']
    if peers['conflicting']:md+=['- Conflicting chromosome assignments in '+', '.join(peers['conflicting'])+'.']
    sisters=[t for t in event.get('peer_tests',[]) if t.get('sister')]
    if sisters:md+=['','**Other haplotype:** '+ '; '.join(t['peer']+': '+str(t['left']['chrom'] or '?')+' → '+str(t['right']['chrom'] or '?')+(' (qualified local assignments)' if t['qualified'] else ' (incomplete local assignments)') for t in sisters)+'. This is within-individual context, not an independent vote.','']
    md+=['','**Cut options**','']
    options=event.get('cut_options',[])
    if options:
        md+=mdtable(['Option','Gap / cut position','Individuals separating flanks','HiFi crossing / left / right','Hi-C contact loss','Use'],[
            [c['id'],bp(c['lo'])+'–'+bp(c['hi'])+'; cut '+bp(c['cut_bp']),len(c['peers']['separate']),str(c.get('spanning','Not measured'))+' / '+' / '.join(map(str,c.get('flank_molecules',['?','?']))),
             ', '.join(t['library'] for t in c.get('contacts',[]) if t['loss']) or 'Not established',
             'Recommended cut' if event.get('recommended_cut',{}).get('id')==c['id'] else 'Not recommended' if event['recommendation']=='RETAIN' else 'Manual review option'] for c in options])
        md+=['','A gap outside the transition range is a nearby alternative, not a proven localization. Do not cut every option. All options begin unselected.','']
        for option in options:
            if option.get('contacts'):
                md+=['**'+option['id']+' — contact evidence at this exact gap**','']+mdtable(['Library','Cross / within-left / within-right','Matched controls','Finding'],[
                    [t['library'],str(t['raw_counts']['cross'])+' / '+str(t['raw_counts']['left_within'])+' / '+str(t['raw_counts']['right_within']),t['matched_controls'],
                     'Supports cutting: contact loss' if t['loss'] else 'No demonstrated loss' if t['calibrated'] else 'Comparison not calibrated'] for t in option['contacts']])+['']
            if not option.get('bridge_consistent',True):md+=['The chromosome assignments between '+option['id']+' and the transition do not establish this gap as its cut site.','']
            if option.get('native_continuity'):md+=['Both sides of '+option['id']+' map within native contig '+option.get('native_left','?')+'; this argues against treating it as a simple scaffolding join.','']
    else:md+=['No scaffold gap in or within 250 kb of this transition is available. An internal sequence cut needs an exact, reviewed coordinate; the range midpoint is not a suggested cut.','']
    if event.get('contacts'):
        md+=['**Hi-C measurements at the transition**','']+mdtable(['Library','Cross / within-left / within-right pairs','Matched controls','Interpretation'],[
            [t['library'],str(t['raw_counts']['cross'])+' / '+str(t['raw_counts']['left_within'])+' / '+str(t['raw_counts']['right_within']),t['matched_controls'],
             'Contact loss' if t['loss'] else 'No clear contact loss' if t['calibrated'] else 'Too few usable contacts or matched controls'] for t in event['contacts']])
        md+=['','Hi-C contributed to scaffolding and is corroborating evidence. Measured zero counts are shown as 0; unavailable measurements are not zeros.','']
    if event.get('hifi_strict',{}).get('probes'):
        strict=event['hifi_strict'];md+=['**Stricter read check:** '+format(100*strict['supported_fraction'],'.1f')+'% of positions supported at MAPQ ≥30. A difference between tiers can indicate mapping ambiguity; it is not automatically a break.','']
    if event.get('hifi_all_primary',{}).get('probes'):
        md+=['**Mapping ambiguity check:** '+format(100*event['hifi_all_primary']['supported_fraction'],'.1f')+'% of positions are crossed when low-confidence primary placements are also included. This helps distinguish absent reads from uncertain placement; these reads do not authorize cutting or establish continuity.','']
    if grid.get('weak_runs') and event['recommendation']!='RETAIN':
        md+=['**Read support localization**','']+mdtable(['Weak interval (bp)','Minimum spanning molecules'],[[bp(w['start'])+'–'+bp(w['end']),w['minimum']] for w in grid['weak_runs']])+['','These are mapping support intervals, not selected cuts. An internal cut still requires a reviewed coordinate within the evidence interval.','']
    if cohort:md+=['[Evidence plots and peer-by-peer comparisons](assemblies/'+assembly+'/report.md#'+event['id'].lower()+')','']
    return md


def decisions(data):
    rows=[]
    for event in data['events']:
        targets=event['cut_options'] or [dict(id=event['id'],lo=event['lo'],hi=event['hi'])]
        for target in targets:
            gap='cut_bp' in target
            rows.append(dict(selected='NO',review_disposition='PENDING',id=target['id'],assembly=data['assembly'],scaffold=event['scaffold'],
                chromosomes=event['left']+' → '+event['right'],recommendation=advice(event) if not gap or event.get('recommended_cut',{}).get('id')==target['id'] or event['recommendation']=='RETAIN' else 'Manual review option; no recommended cut here',cut_bp=target.get('cut_bp',''),review_start=event['lo'],review_end=event['hi'],
                action='UNJOIN_UNSUPPORTED' if gap else 'UNRESOLVED',gap_start=target['lo'] if gap else '',gap_end=target['hi'] if gap else '',localization_status='available_gap' if gap else 'unlocalized',
                reviewer='',reason='',coordinate_stage='pre_finishing',assessment_sha256=data['assessment_sha256'],decision_source='review',evidence_packet_id=data['assembly']))
    return rows


INSTRUCTIONS='''# Review and apply cuts

Copy `review.tsv` outside generated output before editing. Recommendations are advice; all selections start at NO.

| Decision | review_disposition | selected | Required |
| --- | --- | --- | --- |
| Cut | CUT | YES | reviewer, reason, exact cut_bp and cutting action |
| Keep | RETAIN | NO | reviewer and reason |
| Decide later | DEFER | NO | reviewer and reason |
| Unreviewed | PENDING | NO | Nothing |

For an approved gap cut, keep its supplied cut_bp, gap_start, gap_end and action=UNJOIN_UNSUPPORTED. Select one option per intended cut, not every nearby option. A retained event's gap options remain unselected.

For an internal cut, supply an exact cut_bp, action=BREAK_PROBABLE_MISJOIN, localization_status=localized and blank gap fields. Do not use a broad range midpoint without evidence.

To add an undetected cut, append a row with a new unique id. Copy assembly, assessment_sha256 and coordinate_stage from assembly_registry.tsv. Enter the original scaffold, exact coordinates, reviewer, reason, decision_source=review and evidence_packet_id=assembly. Set CUT and YES. Multiple cuts use separate rows, all in original coordinates.

Coordinates are zero-based on the original pre-finishing scaffold. cut_bp divides [0,cut_bp) and [cut_bp,length). All sequence, including gap Ns, is preserved.

Run the same full pipeline with `--chimera_break /absolute/path/reviewed.tsv -resume`. The cutter validates checksums, scaffold identity, literal gaps, conflicts and minimum piece lengths. It emits corrected assemblies, a cut audit, coordinate map and exact reconstruction verification. Chromosome reassignment and finishing use the corrected sequences. Automated cutting remains deferred.
'''


def main():
    p=argparse.ArgumentParser();p.add_argument('--packet',action='append',required=True);p.add_argument('--actions');p.add_argument('--verification',action='append',default=[]);a=p.parse_args()
    all_data=[];rows=[];registry=[];manifest={};versions=[]
    applied=[json.loads(Path(v).read_text()) for v in a.verification]
    reviewed=table(a.actions) if a.actions else []
    for value in sorted(a.packet):
        packet=Path(value);data=json.loads((packet/'evidence.json').read_text());assembly=data['assembly']
        if any(d['assembly']==assembly for d in all_data):raise ValueError('Duplicate assembly packet')
        target=Path('assemblies')/assembly;target.parent.mkdir(exist_ok=True)
        shutil.copytree(packet,target);all_data.append(data);rows+=decisions(data)
        versions.extend(dict(process='MISASSEMBLY_ASSESS',tool=tool,version=version) for tool,version in data.get('tools',{}).items())
        mapping=data['mapping']
        if mapping.get('bam'):manifest[assembly]=dict(bam=mapping['bam'],index=mapping['index'],fasta=data.get('assessment_fasta'),provenance=str((packet/'evidence.json').resolve()))
        result=str(len(data['events']))+' potential chimera(s)' if data['events'] else 'No potential chimeras detected' if data['status']=='assessed' else 'No chromosome-scale scaffolds identified' if data['status']=='no_chromosome_scope' else 'Not assessed — no informative cohort comparisons'
        if data['status']=='no_chromosome_scope' and data.get('scope_unresolved'):result='Chromosome-scale classification unresolved'
        registry.append(dict(assembly=assembly,assessment_sha256=data['assessment_sha256'],assessment_fasta=data.get('assessment_fasta',''),coordinate_stage='pre_finishing',assessment_status=data['status'],potential_chimeras=len(data['events']),result=result))
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
        for event in data['events']:
            md+=['<a name="'+event['id'].lower()+'"></a>','']+event_card(event,assembly)
            if (target/'plots'/(event['id']+'.png')).exists():md+=['![Chromosome identity and local HiFi support](plots/'+event['id']+'.png)','','Shading marks the chromosome transition range. Dashed red lines mark available gap cut options. Blank chromosome tracks mean no usable assignment; gaps in read curves mean unmeasured positions.','']
            md+=['<details>','<summary>Peer chromosome assignments</summary>','']+mdtable(['Assembly','Left → right','Local assignment','Role'],[
                [t['peer'],str(t['left']['chrom'] or '?')+' → '+str(t['right']['chrom'] or '?'),'Confirmed' if t['qualified'] else 'Partial or ambiguous',
                 'Sister haplotype' if t['sample']==data['sample'] else 'Independent individual' if t['eligible'] else 'Context only'] for t in event['peer_tests']])+['','</details>','']
        scaffold_figures=sorted((target/'plots').glob('*.scaffold.png'))
        if scaffold_figures:
            md+=['## Scaffold contact maps and telomeres','']
            for image in scaffold_figures:md+=['!['+image.stem+'](plots/'+image.name+')','']
        md+=['## Evidence records','', '[Read support](read_support.tsv) · [Hi-C counts](contact_support.tsv) · [Chromosome assignments](chromosome_blocks.tsv) · [Telomere counts](telomeres.tsv) · [Full measurements and provenance](evidence.json).', '',
             'Small chromosome-labelled islands below the chromosome-arm threshold are retained in the machine-readable screening audit. They are not counted as chromosome-sized fusions. Missing homology is distinguished from evidence for a different placement.']
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
    events=[(d,e) for d in all_data for e in d['events']]
    counts={k:sum(e['recommendation']==k for _,e in events) for k in ('CUT','RETAIN','SUSPECT','REVIEW')}
    md=['# Chimera and misassembly review','',str(len(all_data))+' assemblies included; '+str(sum(d['status']=='assessed' for d in all_data))+' have usable cohort comparisons. **'+str(counts['CUT'])+' recommended cuts; '+str(counts['RETAIN'])+' supported retentions; '+str(counts['SUSPECT'])+' suspect transitions; '+str(counts['REVIEW'])+' other reviews.**','',
        'No cuts are selected automatically. [Editable decisions](review.tsv) · [Cut instructions](cut-instructions.md)','',
        '## Assembly overview','']+mdtable(['Assembly','Result','Evidence'],[[r['assembly'],r['result'],'[Report](assemblies/'+r['assembly']+'/report.md)'] for r in registry])
    if reviewed:
        md+=['','## Submitted manual decisions','',str(sum(r.get('selected')=='YES' for r in reviewed))+' cuts selected; '+str(sum(len(v['selected_actions']) for v in applied))+' actions verified by exact sequence reconstruction. These evidence measurements describe the original pre-finishing assemblies. Corrected assemblies are reassigned and finished downstream.','']
    if events:
        md+=['','## Suggested actions','']+mdtable(['Assembly / event','Chromosomes','Suggested action','Cut or transition range'],[
            ['['+d['assembly']+' '+e['id']+'](assemblies/'+d['assembly']+'/report.md#'+e['id'].lower()+')',e['left']+' → '+e['right'],advice(e),
             bp(e['recommended_cut']['cut_bp'])+' bp on '+e['scaffold'] if e.get('recommended_cut') else e['scaffold']+': '+bp(e['lo'])+'–'+bp(e['hi'])+' bp'] for d,e in events])
        md+=['','## Evidence for each decision','']
        for d,e in events:md+=event_card(e,d['assembly'],cohort=True)
    md+=['## Scope and reproducibility','',
         'Discovery compares all chromosome-scale scaffolds with the other assemblies, including screening the selected naming reference. The assigned-sequence arm threshold is recorded in each evidence packet; short islands remain in the audit. Individuals, not haplotypes, supply independent votes.', '',
         'HiFi continuity uses a grid of 1 kb probes with 1 kb flanks at MAPQ 20 and 30. Indels up to 50 bp are tolerated; larger indels and skipped regions interrupt support. Known scaffold N gaps are assessed separately. Retention can be supported by continuous or ≥95% supported regions with weak pockets no longer than 3 kb. These are explicit review heuristics, not proof of evolutionary fusion.', '',
         '[Assembly registry](assembly_registry.tsv) · [Retained mapping manifest](mapping_manifest.json). Policy: '+POLICY+'.']
    Path('README.md').write_text('\n'.join(md)+'\n',encoding='utf-8');Path('cut-instructions.md').write_text(INSTRUCTIONS,encoding='utf-8')
    Path('coordinate_status.md').write_text('# Misassembly coordinates\n\nAll discovery, evidence and manual actions refer to source-bound pre-finishing scaffolds. Cuts are verified by exact sequence reconstruction.\n',encoding='utf-8')


if __name__=='__main__':main()
