"""Create source-bound manual review records and compact Markdown evidence packets."""
import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path
from chimera_markdown import render, decision_context, assess_transitions
from chimera_decisions import advise

FIELDS=['selected','id','assembly','coordinate_stage','assessment_sha256','scaffold','action','cut_bp','gap_start','gap_end','review_start','review_end','review_range','localization_explanation',
        'decision_source','evidence_packet_id','reviewer','reason','localization_status','evidence_summary','report_path','source_candidate_id','chromosome_context','evidence_for_cut','evidence_against_cut','review_priority','evidence_limits',
        'detected_transition','assessment_status','transition_id','preferred_candidate','bridge_status','related_candidate','review_disposition','structural_context','local_confirmation','partial_confirmation','sister_context',
        'suggested_action','recommendation','recommendation_reason']


def write_table(path,rows):
    with path.open('w',newline='',encoding='utf-8') as handle:
        writer=csv.DictWriter(handle,fieldnames=FIELDS,delimiter='\t');writer.writeheader();writer.writerows(rows)


def generate(assembly,calls,context,out,coordinate_audit=None,supplement=None,name_map=None,supplement_status='not recorded'):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);context=Path(context)
    def read(name):
        path=context/name
        return json.loads(path.read_text(encoding='utf-8')) if path.exists() else {}
    provenance=read('provenance.json');measurements=read('decision_measurements.json')
    provenance['supplementary_status']=supplement_status
    if coordinate_audit:
        with open(coordinate_audit,encoding='utf-8') as handle:
            provenance['coordinate_audit']={r['metric']:r['value'] for r in csv.DictReader(handle,delimiter='\t')}
    if name_map and Path(name_map).name!='NO_HARMONIZE':
        with open(name_map,encoding='utf-8') as handle:
            scope=list(csv.DictReader((line for line in handle if not line.startswith('#')),delimiter='\t'))
        provenance['chromosome_scope']=dict(eligible_scaffolds=sum(r.get('chromosome_member')=='yes' for r in scope),
            unresolved=any(r.get('chromosome_member')=='unresolved' for r in scope),
            reasons=sorted({r.get('chromosome_scope_reason','not recorded') for r in scope}))
    with open(calls,encoding='utf-8') as handle:
        records=list(csv.DictReader((line for line in handle if not line.startswith('#')),delimiter='\t'))
    candidates={hashlib.sha256(json.dumps(r,sort_keys=True).encode()).hexdigest()[:20]:r for r in records}
    for key,m in measurements.items():
        if m.get('review_only'):candidates[key]=dict(scaffold=m['scaffold'],assembly_sha256=m['assessment_sha256'])
    intervals=provenance.get('intervals',{})
    def location(item):
        key,c=item;m=measurements.get(key,{});interval=intervals.get(m.get('packet_interval_id'),{})
        digits=''.join(x for x in c['scaffold'] if x.isdigit())
        return (int(digits) if digits else 0,c['scaffold'],int(m.get('cut_bp') or c.get('transition_lo') or interval.get('lo') or 0),key)
    rows=[]
    for number,(key,c) in enumerate(sorted(candidates.items(),key=location),1):
        m=measurements.get(key,{});gap=m.get('verified_gap') is True
        m['detector_context']={k:c[k] for k in ('arm_pattern','left_arm_aligned_bp','right_arm_aligned_bp','left_arm_fraction','right_arm_fraction') if k in c and c[k]!='.'}
        interval=intervals.get(m.get('packet_interval_id'),{})
        lo=m.get('gap_start') if gap else c.get('transition_lo',interval.get('lo'))
        hi=m.get('gap_end') if gap else c.get('transition_hi',interval.get('hi'))
        lo='' if lo in (None,'.') else lo;hi='' if hi in (None,'.') else hi
        row=dict.fromkeys(FIELDS,'')
        row.update(selected='NO',id='C%02d'%number,assembly=assembly,coordinate_stage='pre_finishing',
            assessment_sha256=c.get('assembly_sha256',m.get('assessment_sha256',provenance.get('sha256',''))),
            scaffold=c['scaffold'],action='UNJOIN_UNSUPPORTED' if gap else 'UNRESOLVED',
            cut_bp=m.get('cut_bp','') if gap else '',gap_start=m.get('gap_start',''),gap_end=m.get('gap_end',''),
            review_start=lo,review_end=hi,review_range=('%s–%s'%(lo,hi)) if lo!='' and hi!='' else 'Interval unavailable',
            localization_status='available_gap' if gap else 'unlocalized',decision_source='review',review_disposition='PENDING',
            evidence_packet_id=assembly+'.sequence_context',source_candidate_id=key)
        row['report_path']=assembly+'.review/report.md#candidate-'+row['id'].lower()
        row.update(decision_context(m,provenance.get('sample',assembly.rsplit('_hap',1)[0])))
        row['structural_context']=('Internal to one native primary contig' if m.get('native_continuity') else
            'Between distinct native contigs' if m.get('native_left') not in (None,'.') and m.get('native_right') not in (None,'.') and m['native_left']!=m['native_right'] else 'Native placement unresolved')
        row['localization_explanation']=('Verified all-N gap; gap-end coordinate is available, not an established biological breakpoint' if gap else
            'Measured chromosome-transition interval; no unique failed seam has been established')
        rows.append(row)
    assess_transitions(rows,measurements,provenance,candidates)
    for row in rows:
        row.update(advise(row,measurements.get(row['source_candidate_id'],{}),provenance.get('sample',assembly.rsplit('_hap',1)[0]))[0])
        row['evidence_summary']=row['review_priority']+'; detected '+(row['detected_transition'] or row['chromosome_context'])+'; '+row['local_confirmation']
    if supplement:
        dest=out/'scaffold_evidence';dest.mkdir(exist_ok=True)
        for path in Path(supplement).iterdir():
            if path.is_file() and path.name.startswith(assembly+'.') and path.suffix in ('.png','.tsv'):
                shutil.copyfile(path,dest/path.name)
    write_table(out/'review.tsv',rows)
    (out/'measurements.json').write_text(json.dumps(measurements,indent=2),encoding='utf-8')
    (out/'assessment_provenance.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
    render(assembly,rows,measurements,provenance,[],context,out)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    for key in ('assembly','calls','context','out'):parser.add_argument('--'+key,required=True)
    parser.add_argument('--coordinate-audit');parser.add_argument('--supplement');parser.add_argument('--name-map')
    parser.add_argument('--supplement-status',choices=('enabled','disabled'),default='enabled')
    args=parser.parse_args();generate(args.assembly,args.calls,args.context,args.out,args.coordinate_audit,args.supplement,args.name_map,args.supplement_status)
