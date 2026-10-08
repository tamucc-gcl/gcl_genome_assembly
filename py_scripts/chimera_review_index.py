"""Merge assessment packets into a cohort Markdown report and editable cut file."""
import argparse,csv,json
from pathlib import Path

def main():
    p=argparse.ArgumentParser();p.add_argument('--packet',action='append',required=True);a=p.parse_args()
    rows=[];fields=None;seen=set();registry=[]
    for packet in a.packet:
        packet=Path(packet);registry.append(json.loads((packet/'registry.json').read_text()))
        with (packet/'review.tsv').open(encoding='utf-8') as h:
            reader=csv.DictReader(h,delimiter='\t')
            if fields is None:fields=reader.fieldnames
            if reader.fieldnames!=fields:raise ValueError('Review headers differ')
            for row in reader:
                identity=(row['assembly'],row['id'])
                if identity in seen:raise ValueError('Duplicate candidate identity')
                seen.add(identity);rows.append(row)
    rows.sort(key=lambda r:(r['assembly'],r['scaffold'],r['id']))
    registry.sort(key=lambda r:r['assembly'])
    if len({r['assembly'] for r in registry})!=len(registry):raise ValueError('Duplicate assessed assembly')
    with open('chimera_review.tsv','w',newline='',encoding='utf-8') as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
    with open('assembly_registry.tsv','w',newline='',encoding='utf-8') as h:
        w=csv.DictWriter(h,fieldnames=['assembly','assessment_sha256','coordinate_stage','candidate_count'],delimiter='\t');w.writeheader();w.writerows(registry)
    md=['# Chimera Detection and Review','','## Assessment summary','',
        str(len(registry))+' assemblies assessed; '+str(len(rows))+' candidate boundaries in '+str(len({r['assembly'] for r in rows}))+' assemblies; '+str(sum(r['assessment_status']=='review_required' for r in rows))+' boundaries require review; '+str(sum(r['assessment_status']=='retention_supported' for r in rows))+' have evidence favoring retention. Evidence generation applies no cuts. All review selections start at NO.','',
        '[Editable review TSV](chimera_review.tsv) · [All assessment identities](assembly_registry.tsv) · [Cut instructions](cut-instructions.md)','','| Assembly | Candidate boundaries | Evidence |','| --- | ---: | --- |']
    for r in registry:md.append('| '+r['assembly']+' | '+str(r['candidate_count'])+' | [Report]('+r['assembly']+'.review/report.md) |')
    md+=['','## Candidate boundaries','','| Candidate | Assembly | Scaffold | Chromosomes left → right | Region to review, bp | Exact cut, bp | Review priority | For cutting | Against cutting | Related measurement / bridge |','| --- | --- | --- | --- | --- | ---: | --- | --- | --- | --- |']
    for r in rows:
        if r['assessment_status']!='review_required':continue
        vals=['['+r['id']+']('+r['report_path']+')',r['assembly'],r['scaffold'],r['chromosome_context'],r['review_range'],r['cut_bp'] or 'Not assigned',r['review_priority'],r['evidence_for_cut'],r['evidence_against_cut'],(r['related_candidate'] or 'None')+' / '+r['bridge_status']]
        md.append('| '+' | '.join(str(v).replace('|','&#124;').replace('\n',' ') for v in vals)+' |')
    for status,title in [('retention_supported','Boundaries with evidence favoring retention'),('supporting_measurement','Supporting transition measurements')]:
        md+=['','## '+title,'','| Assembly | Measurement | Transition | Preferred candidate | Bridge assessment |','| --- | --- | --- | --- | --- |']
        for r in rows:
            if r['assessment_status']==status:md.append('| '+' | '.join([r['assembly'],'['+r['id']+']('+r['report_path']+')',r['chromosome_context'],r['preferred_candidate'],r['bridge_status']])+' |')
    md+=['','IDs C01, C02, etc. are scoped to an assembly. Separate-chromosome assignments support reviewing a fusion; same-chromosome assignments oppose a fusion at that sampled boundary. No informative assignment means the measurement cannot decide; discordant within-individual assignments require investigation. Prioritize gap-cut review means a verified gap and at least two informative independent individuals with separate-chromosome assignments, without opposing or discordant chromosome assignments. This is a review aid, not a calibrated cut classifier.','','No detected candidate is not a guarantee of structural correctness. Review ranges are measured chromosome-transition intervals, not permission to cut at their midpoint. An exact coordinate is needed only when requesting a cut. No opposing evidence observed is not equivalent to positive evidence for cutting. Detailed reports distinguish independent chromosome context, local support and measurement limits.']
    Path('README.md').write_text('\n'.join(md)+'\n',encoding='utf-8')
    Path('cut-instructions.md').write_text("""# Specify manual cuts

Copy chimera_review.tsv outside the generated outputs before editing. All coordinates use the original pre-finishing FASTA, not final chromosome names. Coordinates are zero-based; cut_bp divides [0,cut_bp) and [cut_bp,length).

## Detected boundaries

Set selected=YES only for a reviewed cut. Supply reviewer and reason; preserve assembly, assessment_sha256, coordinate_stage, id and evidence identifiers. UNJOIN_UNSUPPORTED requires a literal all-N gap, gap_start/gap_end and a cut within its edges. BREAK_PROBABLE_MISJOIN requires an exact cut_bp and localization_status=localized; leave gap bounds blank. Keep decision_source=review. Keep retained or unresolved rows at NO and record the review rationale. Proposed coordinates alone are not permission to cut.

## Add an undetected cut

Append a row to the same TSV with a new unique id. Copy assembly, assessment_sha256 and coordinate_stage from assembly_registry.tsv, including for assemblies with zero candidates. Enter an exact original scaffold ID and cut_bp; choose the appropriate action and gap fields as above. Supply reviewer, reason, decision_source=review and an evidence_packet_id identifying your supporting evidence. report_path can link your review note. User-added cuts do not require a detector candidate, but do not automatically trigger new local evidence assays. Multiple cuts use separate rows in original coordinates.

## Apply and verify

Repeat the original pipeline command with the same inputs and add:

```bash
--chimera_break /absolute/path/reviewed-cuts.tsv -resume
```

Source checksums, scaffold identities, gap bounds, coordinates, conflicting cuts and minimum piece lengths are validated. Every base is preserved; cutting at a gap's end retains its Ns on the left. Cutting writes an audit, exact reconstruction verification and coordinate lift, then updates chromosome assignment for affected cohorts. Original evidence BAMs remain in original assessment coordinates.
""",encoding='utf-8')

if __name__=='__main__':main()
