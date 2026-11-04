"""Real R renderer checks, including first-pass/reviewed and disabled sections."""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys


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
        for path in (actual/'pangenome').rglob('sharing_summary.tsv'):
            artifacts.append(dict(type='pangenome_summary',id=path.relative_to(actual).parts[1],
                                  path=path.relative_to(actual).as_posix(),content=path.read_text()))
        for path in (actual/'pangenome').rglob('source_sequences.tsv'):
            artifacts.append(dict(type='graph_review_rows',id=path.relative_to(actual).parts[1],path=path.relative_to(actual).as_posix(),content=path.read_text()))
        for path in (actual/'pangenome').rglob('source_audit.json'):
            artifacts.append(dict(type='graph_review_audit',id=path.relative_to(actual).parts[1],path=path.relative_to(actual).as_posix(),content=path.read_text()))
        for path in (actual/'contact_maps').glob('*.png'):
            aid=re.sub(r'_final_.*','',path.name)
            artifacts.append(dict(type='contact_map',id=aid,path=path.relative_to(actual).as_posix()))
        for path in (actual/'organelle').glob('*_mito_circular.png'):
            artifacts.append(dict(type='organelle',id=path.name.removesuffix('_mito_circular.png'),path=path.relative_to(actual).as_posix()))
        for path in (actual/'snail_plots').glob('*_final_snail.svg'):
            artifacts.append(dict(type='snail',id=path.name.removesuffix('_final_snail.svg'),path=path.relative_to(actual).as_posix()))
        for path in (actual/'organelle').glob('*mito_stats.tsv'):
            artifacts.append(dict(type='organelle_stats',id=path.name,path=path.relative_to(actual).as_posix(),content=path.read_text()))
        for path in (actual/'organelle').glob('*mito_status.tsv'):
            artifacts.append(dict(type='organelle_recovery',id=path.name.removesuffix('_mito_status.tsv'),organelle='mitochondrion',
                                  path=path.relative_to(actual).as_posix(),content=path.read_text()))
        registry_path=actual/'assembly/chimeras/assembly_registry.tsv'
        artifacts.append(dict(type='chimera_registry',id='run',path=registry_path.relative_to(actual).as_posix(),content=registry_path.read_text()))
        for path in (actual/'provenance').iterdir():
            if path.is_file():artifacts.append(dict(type='provenance',id=path.name,path=path.relative_to(actual).as_posix()))
        audit_dir=out/'reports/organelle_annotations';audit_dir.mkdir(parents=True,exist_ok=True)
        for path in (actual/'organelle').glob('*.gb'):
            aid=path.name.removesuffix('_mitogenome.gb')
            prefix=audit_dir/path.stem
            subprocess.run([sys.executable,str(repo/'py_scripts/annotation_counts.py'),'--genbank',str(path),
                            '--sample',aid,'--prefix',str(prefix)],check=True)
            for kind,suffix in [('annotation_counts','.counts.tsv'),('annotation_loci','.loci.tsv')]:
                product=Path(str(prefix)+suffix)
                artifacts.append(dict(type=kind,id=aid,path=product.relative_to(out).as_posix(),content=product.read_text()))
        registry=actual/'reports/main/artifact_registry.tsv'
        if registry.exists():
            existing={a['path'] for a in artifacts}
            for item in csv.DictReader(registry.open(),delimiter='\t'):
                if item['path'] not in existing:
                    artifacts.append(item)
        expected=[]
        for sample in sorted({a['sample'] for a in data['assemblies']}):
            expected.append(dict(sample=sample,status='accepted',expected=sum(a['sample']==sample for a in data['assemblies']),reason='Accepted input'))
        # Recover skipped inputs from the retained report for this read-only rendering check.
        original_report=(actual/'assembly_report.md').read_text()
        for sample,status,observed,count,reason in re.findall(r'^\| ([^|]+) \| (skipped|rejected|incomplete) \| (\d+) / (\d+) \| ([^|]+) \|$',original_report,re.M):
            expected.append(dict(sample=sample,status=status,expected=int(count),reason=reason))
        painting=sorted((actual/'reports/chromosome_painting').rglob('*.painting.tsv'))
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
                                       ('merqury','kmer_completeness',80,81),('quast','N50',1000,1000),
                                       ('mapped_reads','coverage',30,31),('hic_contact','retention_pct',0,20),('hic_contact','trans_to_cis',0,.5)]:
                w.writerow(['sample','final',analysis,metric,a,b,97 if metric=='kmer_completeness' and case!='missing_combined' else 'NA'])
        decisions=out/'review.tsv'
        columns=['assembly','id','event_id','selected','review_disposition','cut_bp','assessment_sha256','chromosomes','scaffold','reason']
        decision_rows=[['sample_hap1','C01','C01','NO','PENDING','','hash','chr1 → chr2','old',''],
                       ['sample_hap1','C01G01','C01','YES','CUT','500','hash','chr1 → chr2','old','Reviewed cut']]
        with decisions.open('w',newline='') as handle:
            w=csv.writer(handle,delimiter='\t');w.writerow(columns);w.writerows(decision_rows)
        v=out/'verification.json';v.write_text(json.dumps(dict(exact_parent_reconstruction=True,selected_actions=[dict(zip(columns,decision_rows[1]))])))
        verifications=[v];artifacts=[]
        if case=='unplaced':
            for nm in names:
                nm.write_text(nm.read_text().replace('chr1_1','unplaced_1'))
            for tsv in painting:
                tsv.write_text(tsv.read_text().replace('chr1_1','unplaced_1'))
        if case=='primary':decisions=empty;verifications=[empty]
        if case=='no_recovery':artifacts=[dict(type='organelle_recovery',id='sample',organelle='mitochondrion',path='organelle/status.tsv',content='sample_id\tstatus\nsample\tno_recovery\n')]
        if case=='disabled':qc=empty;decisions=empty;painting=[];verifications=[empty]
        if case in ('first_pass','deferred'):
            if case=='deferred':
                decision_rows[0][4]='DEFER';decision_rows[0][-1]='Waiting for additional evidence'
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
    if actual:
        spec['run']=json.loads((actual/'reports/main/report_manifest.json').read_text())['run']
        controls=json.loads((actual/'provenance/run_provenance.json').read_text())['context']['parameters']
        spec['run'].update(contract_enabled=True,comparisons_enabled=controls['run_post_assembly'] and controls['run_pairwise_alignments'],
            comparison_mode=controls['pairwise_alignment_mode'],graph_enabled=controls['run_pangenome'],graph_validate_only=controls['pangenome_validate_only'],
            contact_maps_enabled=controls['run_final_contact_maps'],contact_resolutions=controls['hic_plot_resolutions'],
            samples=[dict(sample=s,hifi=True,hic=True,shortread=False) for s in sorted({r['sample'] for r in data['assemblies']})])
    if case=='missing_outputs':spec['run'].update(contract_enabled=True,qc_mode='none',comparisons_enabled=True,comparison_mode='within_sample')
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
    headings=re.findall(r'^## (.+)$',text,re.M)
    assert headings.index('Final assemblies, chromosome name maps and mitogenomes') < headings.index('Assembly quality') if case!='disabled' else True
    assert '## Contents' in text
    length_rows=list(csv.DictReader((out/'reports/main/assembly_lengths.tsv').open(),delimiter='\t'))
    source_names={nm.name.removesuffix('.name_map.tsv'):list(csv.DictReader(nm.open(),delimiter='\t')) for nm in names}
    for row in length_rows:
        source=source_names[row['assembly']]
        total=sum(int(r['length']) for r in source)
        assigned=sum(int(r['length']) for r in source if r['new_name'].startswith('chr'))
        assert int(row['assembly_bp'])==total and int(row['chromosome_assigned_bp'])==assigned
        assert math.isclose(float(row['chromosome_assigned_percent']),100*assigned/total)
    if case=='disabled':
        assert '## Assembly quality' not in text and '## Pangenome' not in text and '## Chimera' not in text
    else:
        assert '| Individual | BUSCO complete (%) |' in text
        assert 'Assembly length (Mb)' in text and 'Chromosome-assigned (%)' in text and '**hap1 / hap2**' in text
        if not actual and case!='primary':
            assert ('Not measured (80.0 / 81.0)' if case=='missing_combined' else '97.0 (80.0 / 81.0)') in text
        with qc.open() as handle: source=list(csv.DictReader(handle))
        combined=list(csv.DictReader((out/'reports/main/combined_kmer_values.tsv').open(),delimiter='\t'))
        for row in combined:
            source_row=next(r for r in source if r['sample_id']==row['sample'] and r['stage']=='final' and r['metric']=='kmer_completeness')
            if source_row['both']!='NA':assert math.isclose(float(row['combined']),float(source_row['both']))
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
            assert 'Cut applied and verified' in text and '## Needs attention' not in text
        if case=='first_pass':assert 'Review pending' in text and '## Needs attention' in text
        if case=='deferred':assert 'Review deferred' in text and '## Needs attention' in text
        if case=='wrong_fingerprint':assert 'Cut selected; verification needed' in text and 'Cut applied and verified' not in text
        if case=='manual_additions':
            events=list(csv.DictReader((out/'reports/main/assembly_decisions.tsv').open(),delimiter='\t'))
            assert {r['event_id'] for r in events}=={'M01','M02'}
            assert all(r['status']=='Cut applied and verified' for r in events)
        if case=='no_recovery':assert '## Organelle results' in text and '| sample | mitochondrion | no recovery |' in text
        if case=='primary':assert '| sample | 95.0 | 50.0 | 80.0 |' in text
        read_rows=list(csv.DictReader((out/'reports/main/final_read_support.tsv').open(),delimiter='\t'))
        for row in read_rows:
            match=next(r for r in source if r['sample_id']==row['sample_id'] and r['stage']=='final' and r['analysis']==row['analysis'] and r['metric']==row['metric'])
            assert math.isclose(float(row['value']),float(match[row['qc_haplotype']]))
        assert '## Final read support' in text
    if case=='missing_outputs':
        completeness=list(csv.DictReader((out/'reports/main/report_completeness.tsv').open(),delimiter='\t'))
        assert any(r['Feature']=='dotplot' and r['Status']=='MISSING' for r in completeness)
        assert 'Expected output missing: dotplot' in text
    if painting:
        assert 'rather than a combined cross-assembly consensus' in text
        assert 'not independent alignment evidence' in text
        for ref in data['linear_references']:
            assert 'Painting reference: '+ref['id']+'.' in text
        scope=list(csv.DictReader((out/'reports/main/display_scope.tsv').open(),delimiter='\t'))
        source_names={nm.name.removesuffix('.name_map.tsv'):list(csv.DictReader(nm.open(),delimiter='\t')) for nm in names}
        for row in scope:
            source=source_names[row['assembly']]
            unplaced=[r for r in source if not r['new_name'].startswith('chr')]
            assert int(row['unplaced_pieces'])==len(unplaced)
            assert int(row['unplaced_bp'])==sum(int(r['length']) for r in unplaced)
        assert 'separate scale' in text and 'not additive' in text
        if case=='unplaced':assert all(int(r['shown_scaffolds'])==1 and int(r['unplaced_pieces'])==1 for r in scope)
    if actual:
        events=list(csv.DictReader((out/'reports/main/assembly_decisions.tsv').open(),delimiter='\t'))
        assert sum(r['status']=='Cut applied and verified' for r in events)==4
        assert sum(r['status']=='Retained by review' for r in events)==2
        assert 'scaffold_5:35,426,000' in text
        assert 'Distinct graph sequence (Mb)' in text and 'Complete pangenome report' in text
        assert 'Status: **REVIEW**' not in text
        hub=(out/'reports/main/pangenome_373251.md').read_text()
        assert 'Technical validation and provenance' in hub and 'Status: **REVIEW**' in hub
        contacts=index.split('### Hi-C contact maps',1)[1].split('<a id=',1)[0]
        assert contacts.count('| Sde-CTlk_104_hap1 |')==1
        assert '98.6 (78.5 / 78.4)' in text
        assert '| Individual | hap1 | hap2 | Mitogenome |' in text
        assert '| Individual | hap1 snail | hap2 snail | Mitogenome map |' in text
        assert text.count('alt="Sde-CTlk_104_hap1 final snail plot"')==1
        assert '[GenBank](organelle/Sde-CTlk_104_mitogenome.gb)' in text
        assert index.index('### Dotplots') < index.index('### Ribbon plots') < index.index('### Hi-C contact maps')
        counts=list(csv.DictReader((out/'reports/main/annotation_counts.tsv').open(),delimiter='\t'))
        assert len(counts)==5 and all(int(r['trna_loci'])==22 and int(r['protein_coding_loci'])==13 and int(r['rrna_loci'])==2 for r in counts)
        assert 'Spratelloides delicatulus' in text
        assert (out/'reports/main/assembly_stage_changes.md').is_file()
        completeness=list(csv.DictReader((out/'reports/main/report_completeness.tsv').open(),delimiter='\t'))
        assert not any(r['Status']=='MISSING' for r in completeness),[r for r in completeness if r['Status']=='MISSING']
        attention=list(csv.DictReader((out/'reports/main/needs_attention.tsv').open(),delimiter='\t'))
        assert len(attention)==5
        assert all('CPla' not in r['Item'] for r in attention)
        assert any('Sde-CTlk_104_hap2 / chr9_1'==r['Item'] and '48.9%' in r['Issue'] for r in attention)
        assert all('80%' in r['Issue'] and 'source_report.md' in r['Details'] for r in attention)
        assert all(r['Area']=='Pangenome graph' for r in attention)
        assert 'sequence omitted from the graph' in text and 'do not identify new assembly misjoins or recommend assembly edits' in text
        assert 'Sde-CPla_115_hap1' in (out/'reports/main/assessment_coverage.md').read_text()
        changes=list(csv.DictReader((out/'reports/main/stage_changes.tsv').open(),delimiter='\t'))
        with qc.open() as handle: stage_source=list(csv.DictReader(handle))
        for row in changes:
            assert math.isclose(float(row['change']),float(row['value'])-float(row['previous_value']),abs_tol=1e-9)
            for stage_key,value_key in [('stage','value'),('previous_stage','previous_value')]:
                source_row=next(r for r in stage_source if r['sample_id']==row['sample_id'] and r['stage']==row[stage_key] and
                                r['analysis']==row['analysis'] and r['metric']==row['previous_metric' if stage_key=='previous_stage' else 'metric'])
                expected=float(source_row[row['haplotype']])
                if row['measure']=='BUSCO (pp)':
                    total=next(r for r in stage_source if r['sample_id']==row['sample_id'] and r['stage']==row[stage_key] and
                               r['analysis']=='busco' and r['metric']=='total_busco')
                    expected=100*expected/float(total[row['haplotype']])
                if row['measure'] in ('Length (Mb)','N50 (Mb)'):expected/=1e6
                assert math.isclose(float(row[value_key]),expected,abs_tol=1e-9)
    print('PASS',case,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('rscript');p.add_argument('--actual-root',type=Path)
    a=p.parse_args();root=Path(__file__).resolve().parents[1]/'comparisons/main-report-validation'
    for case in ['reviewed','first_pass','deferred','wrong_fingerprint','manual_additions','no_recovery','primary','unplaced','missing_combined','missing_outputs','disabled']:run(a.rscript,root,case)
    if a.actual_root:run(a.rscript,root,'accepted_run',a.actual_root.resolve())
