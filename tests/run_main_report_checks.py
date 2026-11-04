"""Real R renderer checks, including first-pass/reviewed and disabled sections."""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import re
import subprocess


def run(rscript,root,case,actual=None):
    repo=Path(__file__).resolve().parents[1]
    out=root/case;out.mkdir(parents=True,exist_ok=True)
    empty=out/'NO_REPORT_INPUT';empty.write_text('not produced\n')
    names=[]
    if actual:
        eligibility=actual/'assembly_eligibility.json'
        data=json.loads(eligibility.read_text())
        names=sorted((actual/'assembly/final').glob('*.name_map.tsv'))
        qc=actual/'qc/assembly/assembly_qc_metrics.csv'
        decisions=actual/'assembly/chimeras/review.tsv'
        verifications=list((actual/'assembly/chimeras').glob('*verification.json'))
        artifacts=[dict(type='interactive_qc',id='run',path='reports/assembly_qc_report.html')]
        for path in (actual/'pangenome').rglob('*.md'):
            relative=path.relative_to(actual)
            note=next((line for line in path.read_text().splitlines() if line.startswith('Status:')),'') if path.name=='source_report.md' else ''
            artifacts.append(dict(type='pangenome',id=relative.parts[1],path=relative.as_posix(),note=note))
        for path in (actual/'organelle').glob('*mito_stats.tsv'):
            artifacts.append(dict(type='organelle_stats',id=path.name,path=path.relative_to(actual).as_posix(),content=path.read_text()))
        for path in (actual/'organelle').glob('*mito_status.tsv'):
            artifacts.append(dict(type='organelle_recovery',id=path.name.removesuffix('_mito_status.tsv'),organelle='mitochondrion',
                                  path=path.relative_to(actual).as_posix(),content=path.read_text()))
        expected=[]
        for sample in sorted({a['sample'] for a in data['assemblies']}):
            expected.append(dict(sample=sample,status='accepted',expected=sum(a['sample']==sample for a in data['assemblies']),reason='Accepted input'))
        painting=[]
    else:
        records=[];painting=[]
        for hap in (['primary'] if case=='primary' else ['hap1','hap2']):
            aid='sample_'+hap
            records.append(dict(id=aid,sample='sample',taxid='1',haplotype=hap,completion='finalized',fasta=aid+'.fasta',graph_member=False))
            nm=out/(aid+'.name_map.tsv');nm.write_text('old_name\tnew_name\tlength\tclass\nold\tchr1_1\t1000\tscaffold\n')
            names.append(nm)
            tsv=out/(aid+'.painting.tsv')
            tsv.write_text('sequence\tlength\tstart_0based\tend_exclusive\tchromosome\tbasis\nchr1_1\t1000\t0\t1000\tchr1\tprimary_mapping_chain\n')
            painting.append(tsv)
        data=dict(assemblies=records,linear_references=[dict(id='sample_hap1',taxid='1')])
        eligibility=out/'eligibility.json';eligibility.write_text(json.dumps(data))
        expected=[dict(sample='sample',status='accepted',expected=len(records),reason='Accepted input'),
                  dict(sample='skipped',status='skipped',expected=0,reason='No assembly reads')]
        qc=out/'qc.csv'
        with qc.open('w',newline='') as handle:
            w=csv.writer(handle);w.writerow(['sample_id','stage','analysis','metric','hap1','hap2','both'])
            for analysis,metric,a,b in [('busco','total_busco',100,100),('busco','complete',95,96),('merqury','qv',50,51),
                                       ('merqury','kmer_completeness',80,81),('quast','N50',1000,1000)]:
                w.writerow(['sample','final',analysis,metric,a,b,'NA'])
        decisions=out/'review.tsv'
        columns=['assembly','id','event_id','selected','review_disposition','cut_bp','assessment_sha256','chromosomes','scaffold','reason']
        decision_rows=[['sample_hap1','C01','C01','NO','PENDING','','hash','chr1 → chr2','old',''],
                       ['sample_hap1','C01G01','C01','YES','CUT','500','hash','chr1 → chr2','old','Reviewed cut']]
        with decisions.open('w',newline='') as handle:
            w=csv.writer(handle,delimiter='\t');w.writerow(columns);w.writerows(decision_rows)
        v=out/'verification.json';v.write_text(json.dumps(dict(exact_parent_reconstruction=True,selected_actions=[dict(zip(columns,decision_rows[1]))])))
        verifications=[v];artifacts=[]
        if case=='primary':decisions=empty;verifications=[empty]
        if case=='no_recovery':artifacts=[dict(type='organelle_recovery',id='sample',organelle='mitochondrion',path='organelle/status.tsv',content='sample_id\tstatus\nsample\tno_recovery\n')]
        if case=='disabled':qc=empty;decisions=empty;painting=[];verifications=[empty]
        if case=='first_pass':
            with decisions.open('w',newline='') as handle:
                w=csv.writer(handle,delimiter='\t');w.writerow(columns);w.writerow(decision_rows[0])
            verifications=[empty]
        if case=='wrong_fingerprint':
            payload=json.loads(v.read_text());payload['selected_actions'][0]['assessment_sha256']='other'
            v.write_text(json.dumps(payload))
        if case=='manual_additions':
            decision_rows=[['sample_hap1','M01','','YES','CUT','300','hash','user-added','old','Manual addition'],
                           ['sample_hap1','M02','','YES','CUT','700','hash','user-added','old','Second manual addition']]
            with decisions.open('w',newline='') as handle:
                w=csv.writer(handle,delimiter='\t');w.writerow(columns);w.writerows(decision_rows)
            v.write_text(json.dumps(dict(exact_parent_reconstruction=True,selected_actions=[dict(zip(columns,row)) for row in decision_rows])))
    spec=dict(eligibility=str(eligibility),expected=expected,name_maps=list(map(str,names)),paintings=list(map(str,painting)),
        qc=str(qc),decisions=str(decisions),verifications=list(map(str,verifications)),artifacts=artifacts,versions=str(empty),
        run=dict(run_name='renderer_validation',nextflow='23.10.1',profile='local',top_scaffolds=20,chimera_enabled=case!='disabled'))
    path=out/'spec.json';path.write_text(json.dumps(spec))
    env=dict(os.environ,TZ='UTC');env['PATH']=str(Path(rscript).parent)+os.pathsep+env['PATH']
    result=subprocess.run([rscript,str(repo/'r_scripts/generate_main_report.R'),str(path)],cwd=out,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (out/'renderer.log').write_text(result.stdout)
    assert result.returncode==0,result.stdout
    text=(out/'assembly_report.md').read_text()
    index=(out/'reports/main/evidence_index.md').read_text()
    for target in {item['path'] for item in artifacts}:
        assert index.count('../../'+target)==1, 'Artifact listed in multiple evidence categories: '+target
    assert not re.search(r'!\[[^\]]*[Qq]uality',text)
    if case=='disabled':
        assert '## Assembly quality' not in text and '## Pangenome' not in text and '## Chimera' not in text
    else:
        assert '| Individual | BUSCO complete (%) |' in text
        with qc.open() as handle: source=list(csv.DictReader(handle))
        measured=list(csv.DictReader((out/'reports/main/quality_values.tsv').open(),delimiter='\t'))
        for row in measured:
            match=next(r for r in source if r['sample_id']==row['sample_id'] and r['stage']=='final' and
                       r['analysis']==row['analysis'] and r['metric']==row['metric'])
            value=float(match[row['qc_haplotype']])
            if row['metric']=='complete':
                total=next(r for r in source if r['sample_id']==row['sample_id'] and r['stage']=='final' and
                           r['analysis']=='busco' and r['metric']=='total_busco')
                value=100*value/float(total[row['qc_haplotype']])
            if row['metric']=='N50':value/=1e6
            assert math.isclose(float(row['measurement']),value,rel_tol=1e-12), 'Report changed the metric: '+str(row)
        if case=='reviewed':
            assert 'Cut applied and verified' in text and 'Decisions needing attention' not in text
        if case=='first_pass':assert 'Review pending' in text and 'Decisions needing attention' in text
        if case=='wrong_fingerprint':assert 'Cut selected; verification needed' in text and 'Cut applied and verified' not in text
        if case=='manual_additions':
            events=list(csv.DictReader((out/'reports/main/assembly_decisions.tsv').open(),delimiter='\t'))
            assert {r['event_id'] for r in events}=={'M01','M02'}
            assert all(r['status']=='Cut applied and verified' for r in events)
        if case=='no_recovery':assert '## Organelle results' in text and '| sample | mitochondrion | no recovery |' in text
        if case=='primary':assert '| sample | 95.0 | 50.0 | 80.0 |' in text
    if actual:
        events=list(csv.DictReader((out/'reports/main/assembly_decisions.tsv').open(),delimiter='\t'))
        assert sum(r['status']=='Cut applied and verified' for r in events)==4
        assert sum(r['status']=='Retained by review' for r in events)==2
    print('PASS',case,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('rscript');p.add_argument('--actual-root',type=Path)
    a=p.parse_args();root=Path(__file__).resolve().parents[1]/'comparisons/main-report-validation'
    for case in ['reviewed','first_pass','wrong_fingerprint','manual_additions','no_recovery','primary','disabled']:run(a.rscript,root,case)
    if a.actual_root:run(a.rscript,root,'accepted_run',a.actual_root.resolve())
