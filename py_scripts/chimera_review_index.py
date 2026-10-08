"""Cohort transition families, cross-haplotype context and manual decision interface."""
import argparse
import csv
import json
from pathlib import Path
from chimera_report import families, family_row, table


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--packet',action='append',required=True);args=parser.parse_args()
    rows=[];registry=[];seen=set();fields=None;measurements={};provenance={}
    for value in args.packet:
        packet=Path(value);reg=json.loads((packet/'registry.json').read_text(encoding='utf-8'));registry.append(reg)
        measurements[reg['assembly']]=json.loads((packet/'measurements.json').read_text(encoding='utf-8'))
        provenance[reg['assembly']]=json.loads((packet/'assessment_provenance.json').read_text(encoding='utf-8'))
        with (packet/'review.tsv').open(encoding='utf-8') as handle:
            reader=csv.DictReader(handle,delimiter='\t')
            if fields is None:fields=reader.fieldnames
            if fields!=reader.fieldnames:raise ValueError('Review headers differ')
            for row in reader:
                identity=(row['assembly'],row['id'])
                if identity in seen:raise ValueError('Duplicate candidate identity')
                seen.add(identity);rows.append(row)
    registry.sort(key=lambda r:r['assembly']);rows.sort(key=lambda r:(r['assembly'],r['id']))
    if len({r['assembly'] for r in registry})!=len(registry):raise ValueError('Duplicate assessed assembly')
    for name,values,headers in [('chimera_review.tsv',rows,fields),('assembly_registry.tsv',registry,list(registry[0]))]:
        with open(name,'w',newline='',encoding='utf-8') as handle:
            writer=csv.DictWriter(handle,fieldnames=headers,delimiter='\t');writer.writeheader();writer.writerows(values)
    groups=[]
    for reg in registry:groups+=families([r for r in rows if r['assembly']==reg['assembly']])
    retained=[r for r in rows if r['assessment_status']=='retention_supported']
    md=['# Chimera Detection and Review','','## Assessment summary','',
        str(len(registry))+' assemblies assessed; '+str(len(rows))+' candidate boundaries; '+str(len(groups))+' transition families require review; '+str(len(retained))+' boundaries have evidence against the chromosome-fusion hypothesis.',
        'No cuts selected or applied. Human review dispositions start at PENDING. Related measurements are grouped even when a preferred cut cannot yet be established.','',
        '[Editable decision file](chimera_review.tsv) · [Assessment registry](assembly_registry.tsv) · [Manual cut instructions](cut-instructions.md)','']
    md+=table(['Assembly','Measurements','Assessment / reason','Evidence'],[
        [r['assembly'],r['candidate_count'],r['assessment_status']+': '+r['assessment_reason'],'[Report]('+r['assembly']+'.review/report.md)'] for r in registry])
    if groups:
        md+=['','## Transition families requiring review','']+table(['Assembly','Measurements','Scaffold','Detected chromosomes','Measured intervals, bp','Available gap cut, bp','Independent chromosome observations','Native context'],[family_row(g,True) for g in groups])
        md+=['','A gap-end coordinate is a physically available cut location, not a uniquely established biological breakpoint. Separate chromosome assignments support a fusion hypothesis; they do not by themselves authorize cutting. Below-threshold observations and native path continuity are shown in the detailed report.']
    # Corresponding chromosome combinations are not automatically homologous junctions.
    combinations={}
    for group in groups:
        first=group[0];pair=(first['detected_transition'] or first['chromosome_context']).split(' → ')
        if len(pair)==2 and pair[0]!=pair[1]:combinations.setdefault(tuple(sorted(pair)),[]).append(group)
    if combinations:
        md+=['','## Chromosome combinations across haplotypes and individuals','']
        values=[]
        for pair,matched in sorted(combinations.items()):
            if len(matched)<2:continue
            for group in matched:
                first=group[0];sample=provenance[first['assembly']].get('sample',first['assembly'])
                sister=[]
                for row in group:
                    for track in measurements[row['assembly']].get(row['source_candidate_id'],{}).get('chromosome_tracks',[]):
                        if track['sample']==sample:
                            target=track.get('left_target');right=track.get('right_target')
                            lo=track.get('left_target_bp');hi=track.get('right_target_bp')
                            partners=[r for r in rows if r['assembly']==track['peer'] and r['scaffold']==target and target==right and lo is not None and hi is not None and r['review_start']!='' and min(lo,hi)<=int(r['review_start'])<=int(r['review_end'])<=max(lo,hi)]
                            sister.append(track['peer']+': '+('local anchors bracket '+', '.join(r['id'] for r in partners) if partners else 'corresponding local boundary not established'))
                values.append([' + '.join(pair),first['assembly'],first['scaffold'],' / '.join(r['id'] for r in group),'; '.join(dict.fromkeys(sister)) or 'No sister placement assay'])
        if values:md+=table(['Chromosome combination','Assembly','Scaffold','Measurements','Sister-haplotype correspondence'],values)
        md+=['','Matching chromosome combinations across haplotypes are contextual evidence. Homologous local placement, order and orientation must establish whether the same junction is present. Two haplotypes from one individual are never counted as independent biological replication.']
    if retained:
        md+=['','## Boundaries with evidence against a chromosome fusion','']+table(['Assembly','Measurement','Initial detection','Local context','Assessed interval','Interpretation'],[
            [r['assembly'],'['+r['id']+']('+r['report_path']+')',r['detected_transition'],r['chromosome_context'],r['review_range'],'Opposes proposed chromosome transition; exact adjacency not validated'] for r in retained])
    md+=['','## Review and apply decisions','',
        'C01, C02, etc. identify measurements within an assembly. A transition family may contain several measurements. Review each family once; select an exact cut row only after deciding its location. All measurements remain in the TSV and linked evidence.',
        'Record review_disposition=PENDING, RETAIN, CUT or DEFER. selected=YES authorizes a cut only with the required action, coordinates, reviewer and reason. A disposition of RETAIN or DEFER leaves selected=NO.',
        'Zero detected candidates means no qualifying chromosome-scale composite entered the local assays; it does not validate every assembly join. See each assembly’s assessment scope and reason.']
    Path('README.md').write_text('\n'.join(md)+'\n',encoding='utf-8')
    Path('cut-instructions.md').write_text(INSTRUCTIONS,encoding='utf-8')


INSTRUCTIONS='''# Specify manual cuts

Copy chimera_review.tsv outside generated outputs before editing. Coordinates use the original pre-finishing FASTA and original scaffold IDs, not final chromosome names. cut_bp divides [0,cut_bp) and [cut_bp,length). All bases are retained.

## Detected boundaries

Review each transition family, including supporting measurements and opposing evidence. Record review_disposition as PENDING, RETAIN, CUT or DEFER. For RETAIN/DEFER, keep selected=NO and add reviewer/reason. For an approved cut, set selected=YES, review_disposition=CUT, reviewer and reason. Other measurements in the family remain NO. No cuts are inferred from priority or review disposition alone.

UNJOIN_UNSUPPORTED requires verified all-N gap_start/gap_end and cut_bp within its edges. BREAK_PROBABLE_MISJOIN requires an exact cut_bp and localization_status=localized; leave gap fields blank. RETAIN is not a cutting action. If overriding a retention assessment, explicitly change the action and provide the required coordinates and rationale. Keep decision_source=review and preserve source identity fields.

## Add an undetected cut

Append a row with a new unique id. Copy assembly, assessment_sha256 and coordinate_stage from assembly_registry.tsv, including for assemblies with zero candidates. Enter the original scaffold ID, exact cut_bp and appropriate action/gap fields. Supply reviewer, reason, decision_source=review, evidence_packet_id, selected=YES and review_disposition=CUT. Multiple cuts use separate rows in original coordinates. User-added rows do not themselves generate new evidence assays.

## Apply and verify

Repeat the original pipeline command with the same inputs and add:

```bash
--chimera_break /absolute/path/reviewed-cuts.tsv -resume
```

The cut process validates source checksums, scaffold identities, coordinates, literal gaps, conflicts and minimum piece lengths. It writes a cut audit, coordinate lift and exact sequence reconstruction verification. Downstream chromosome assignment and finishing rerun for affected inputs; original evidence remains in pre-finishing coordinates. An all-N gap cut at its end retains the gap Ns on the left piece. Automatic cutting remains deferred.
'''


if __name__=='__main__':main()
