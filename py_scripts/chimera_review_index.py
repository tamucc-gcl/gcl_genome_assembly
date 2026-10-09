"""Cohort transition families, cross-haplotype context and manual decision interface."""
import argparse
import csv
import json
from pathlib import Path
from chimera_report import families, table


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
    from chimera_decisions import assessment_label, advise, decision_card, focus, location
    cut_options=[]
    for group in groups:
        sample=provenance[group[0]['assembly']].get('sample')
        row=focus(group,measurements[group[0]['assembly']],sample)
        if advise(row,measurements[row['assembly']].get(row['source_candidate_id'],{}),sample)[0]['suggested_action']=='CONSIDER_GAP_CUT':cut_options.append(group)
    pending=[g for g in groups if g not in cut_options]
    md=['# Chimera review','',
        '**'+str(len(groups))+' potential chimeras in '+str(len({g[0]['assembly'] for g in groups}))+' assemblies.** '+str(len(cut_options))+' gap cut'+('s' if len(cut_options)!=1 else '')+' suggested for consideration; '+str(len(pending))+' other event'+('s' if len(pending)!=1 else '')+' should remain intact pending review.',
        '', 'No cuts have been selected or applied. Suggestions are for human review.', '',
        '[Editable decision file](chimera_review.tsv) · [How to approve or add a cut](cut-instructions.md)','',
        '## Assembly overview','']
    md+=table(['Assembly','Result','Evidence'],[
        [reg['assembly'],assessment_label(reg,families([r for r in rows if r['assembly']==reg['assembly']])),
         '[Report]('+reg['assembly']+'.review/report.md)'] for reg in registry])
    if not groups:md+=['','0 candidate boundaries require a chromosome-fusion decision.','']
    if groups:
        md+=['','## Suggested actions','']
        values=[]
        for group in cut_options+pending:
            assembly=group[0]['assembly'];sample=provenance[assembly].get('sample');row=focus(group,measurements[assembly],sample)
            advice,_=advise(row,measurements[assembly].get(row['source_candidate_id'],{}),sample)
            values.append([assembly,'['+' / '.join(r['id'] for r in group)+']('+row['report_path']+')',row['scaffold']+' · '+row['detected_transition'],advice['recommendation'],location(row)])
        md+=table(['Assembly','Event','Scaffold / chromosomes','Suggested action','Cut or transition range'],values)
        md+=['', '“Consider cutting” identifies a concrete manual cut option. “Keep for now” means leave the sequence intact because evidence conflicts or an exact cut is not established; it does not mean the fusion is biologically confirmed.', '',
             'Coordinates refer to the original scaffolds before finishing. Related C01 / C02 measurements describe one possible fusion, not two cuts.','',
             '## Evidence for each decision','']
        for group in cut_options+pending:
            md+=decision_card(group,measurements[group[0]['assembly']],provenance[group[0]['assembly']],cohort=True)
    combinations={}
    for group in groups:
        identity=(group[0]['detected_transition'] or '').split(' → ')
        if len(identity)==2:combinations.setdefault(tuple(sorted(identity)),[]).append(group)
    shared=[(identity,matched) for identity,matched in combinations.items() if len(matched)>1]
    if shared:
        md+=['## What the two haplotypes tell us','']
        for identity,matched in shared:
            md+=['- **'+' + '.join(identity)+'** occurs in '+', '.join(g[0]['assembly']+' ('+' / '.join(r['id'] for r in g)+')' for g in matched)+'.']
        md+=['', 'The same chromosome combination in both haplotypes is a reason to consider a real fusion. It does not establish that the exact same join is present: the local sister-haplotype comparisons in this run must resolve that separately. Two haplotypes of one individual are not independent individuals.','']
    if retained:
        md+=['<details>','<summary>Screened out: '+str(len(retained))+' reference-alignment signals — no cuts suggested</summary>','',
             'Keep these joins. Local comparisons place both sides on the same chromosome, so the initial reference-label change was rejected as a chromosome-fusion signal. These are not additional potential chimeras.','']
        md+=table(['Assembly','Record','Why excluded'],[
            [r['assembly'],'['+r['id']+']('+r['report_path']+')',r['chromosome_context']+' in local comparisons'] for r in retained])
        md+=['','The detector deliberately screens reference-label changes before comparing other assemblies. Short reference-labelled segments can create false positives; the local comparison stage removes them from the action list while preserving the audit.','</details>','']
    md+=['## Record and apply decisions','',
        'Edit [chimera_review.tsv](chimera_review.tsv). Its recommendation columns reproduce the advice above; they never select a cut. To approve a cut, set review_disposition=CUT and selected=YES, and supply reviewer and reason. Keep other rows unselected.', '',
        'For an event with only a range, an exact coordinate must be supplied before cutting. See [cut instructions](cut-instructions.md) for the required action and for adding an undetected cut.','',
        '<details>','<summary>Scope and audit files</summary>','',
        '“No potential chimeras detected” applies to this chromosome-fusion screen, not every possible assembly error. Assemblies without chromosome-scale scaffolds or reference alignments are identified separately above.', '',
        '[Assessment registry](assembly_registry.tsv). Full measurements and source checksums are linked from each assembly report.','', '</details>']
    Path('README.md').write_text('\n'.join(md)+'\n',encoding='utf-8')
    Path('cut-instructions.md').write_text(INSTRUCTIONS,encoding='utf-8')


INSTRUCTIONS='''# Specify manual cuts

Copy chimera_review.tsv outside generated outputs before editing. Coordinates use the original pre-finishing FASTA and original scaffold IDs, not final chromosome names. cut_bp divides [0,cut_bp) and [cut_bp,length). All bases are retained.

## Choose what to do

| Your decision | review_disposition | selected | Also supply |
| --- | --- | --- | --- |
| Make the cut | CUT | YES | reviewer, reason and an exact cut_bp |
| Keep the join | RETAIN | NO | reviewer and reason |
| Decide later | DEFER | NO | reviewer and reason |
| Not reviewed | PENDING | NO | Nothing |

The recommendation columns contain pipeline advice. Changing them does not approve a cut. You must set CUT and YES explicitly.

For a suggested gap cut, keep the supplied scaffold, cut_bp, action, gap_start and gap_end. Set CUT, YES, your reviewer name and your reason. Leave the related transition measurement unselected: it is not a second cut.

For a transition with only a range, do not copy the range midpoint into cut_bp. After choosing and justifying an exact coordinate, enter it as cut_bp, set action=BREAK_PROBABLE_MISJOIN and localization_status=localized, and clear gap_start/gap_end. Then record CUT, YES, reviewer and reason.

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
