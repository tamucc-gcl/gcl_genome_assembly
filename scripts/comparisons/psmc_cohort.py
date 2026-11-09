"""Prepare and review a normal all-individual main.nf PSMC acceptance run."""
import argparse
import csv
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
from prepare_integration_params import prepare as baseline_params

REPO=Path(__file__).resolve().parents[2]


def prepare(baseline,out):
    baseline=Path(baseline).resolve();out=Path(out).resolve();out.mkdir(parents=True,exist_ok=True)
    params=baseline_params(baseline,REPO/'nextflow.config',out/'baseline-params.json')
    # Selection and mode belonged to earlier one-individual experiments. Reset
    # them explicitly so the normal workflow plans the entire current cohort.
    params.update(run_psmc=True,psmc_manifest=None,psmc_mode='auto',psmc_assembly_preset='asm20',
        psmc_assembly_min_alignment_bp=10000,psmc_mutation_rate=None,psmc_generation_time=None,
        run_final_hifi_mapping=True)
    (out/'baseline-params.json').write_text(json.dumps(params,indent=2)+'\n')
    eligibility=json.loads((baseline/'assembly_eligibility.json').read_text())
    expected=sorted({r['sample'] for r in eligibility['assemblies']})
    (out/'psmc_expected.json').write_text(json.dumps(dict(samples=expected),indent=2)+'\n')
    print('PSMC cohort selection: '+', '.join(expected))


def audit(root):
    root=Path(root);psmc=root/'results/post_assembly/psmc';plan=json.loads((psmc/'psmc_plan.json').read_text())
    expected=json.loads((root/'psmc_expected.json').read_text())['samples']
    observed=[j['sample'] for j in plan['jobs']]+[j['sample'] for j in plan['skipped']]
    issues=[]
    if sorted(observed)!=expected:issues.append('Planned/skipped individuals differ from complete expected cohort')
    with (psmc/'summary.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    if sorted(r['sample'] for r in rows)!=expected:issues.append('PSMC result index omits or duplicates cohort members')
    for job in plan['jobs']:
        base=psmc/str(job['taxid'])/job['sample'];q=json.loads((base/'input/input_qc.json').read_text())
        if q['reference_sha256']!=job['reference']['sha256']:issues.append(job['sample']+': reference identity mismatch')
        if job['mode']=='assemblies':
            settings=json.loads((base/'genotypes/dipcall_settings.json').read_text())
            params=json.loads((root/'results/provenance/run_provenance.json').read_text())['context']['parameters']
            if settings['assembly_preset']!=params['psmc_assembly_preset']:issues.append(job['sample']+': preset differs from recorded run parameters')
            if settings['confidence_thresholds']['mask_min_alignment_bp']!=int(params['psmc_assembly_min_alignment_bp']) or settings['confidence_thresholds']['sam_filter_min_alignment_bp']!=int(params['psmc_assembly_min_alignment_bp']):
                issues.append(job['sample']+': inconsistent mask/call alignment lengths')
    result=dict(status='FAIL' if issues else 'PASS',expected_samples=expected,planned=len(plan['jobs']),
        skipped=plan['skipped'],results=rows,issues=issues)
    (root/'psmc_cohort_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    if issues:raise ValueError('; '.join(issues))
    return result


def context(root,length_run,hifi_run):
    root=Path(root).resolve();length_run=Path(length_run).resolve();hifi_run=Path(hifi_run).resolve()
    spec=json.loads((length_run/'comparison.json').read_text());dest=root/'psmc_checks';dest.mkdir(exist_ok=True)
    value=dict(repository=str(REPO),sample=spec['sample'],source=spec['source'],
        run_directories=dict(cohort=str(root),control_20kb=spec['run_directories']['asm20_20kb'],hifi=str(hifi_run)))
    (dest/'comparison.json').write_text(json.dumps(value,indent=2)+'\n')
    # Both pairwise audits reuse the already completed own-reference HiFi BAM.
    return value


def environment(root,sample):
    root=Path(root)
    with (root/'trace.tsv').open() as f:trace=list(csv.DictReader(f,delimiter='\t'))
    rows=[r for r in trace if r['name'].split(' (')[0].endswith(':PSMC_INPUT') and r['name'].endswith('('+sample+')')
        and r['status'] in ['COMPLETED','CACHED'] and r['exit']=='0']
    if len(rows)!=1:raise ValueError('Expected one successful PSMC_INPUT for '+sample)
    contexts=json.loads((root/'execution_environments.json').read_text())['tasks']
    selected=[r for r in contexts if r['workdir']==rows[0]['workdir']]
    prefixes={p for r in selected for p in r.get('conda_prefixes',[])}
    if len(prefixes)!=1:raise ValueError('Expected recorded PSMC analysis environment for '+sample)
    return next(iter(prefixes))


def collect(root):
    root=Path(root);out=root/'psmc_checks';out.mkdir(exist_ok=True);issues=[]
    validation=root/'psmc_cohort_validation.json'
    if validation.is_file():
        value=json.loads(validation.read_text());shutil.copyfile(validation,out/validation.name)
        if value['status']!='PASS':issues.extend(value['issues'])
    else:value=dict(results=[]);issues.append('Cohort validation unavailable')
    plan_path=root/'results/post_assembly/psmc/psmc_plan.json'
    if plan_path.exists():
        plan=json.loads(plan_path.read_text())
        for job in plan['jobs']:
            report=out/('support_'+job['sample'])
            if not (report/'exit_status.txt').is_file() or (report/'exit_status.txt').read_text().strip()!='0':issues.append(job['sample']+': direct read-support audit failed or missing')
    for label in ['20kb_vs_cohort','hifi_vs_cohort']:
        report=out/label
        if not (report/'exit_status.txt').is_file() or (report/'exit_status.txt').read_text().strip()!='0':issues.append(label+': pairwise audit failed or missing')
    context_path=out/'comparison.json'
    if context_path.exists():
        spec=json.loads(context_path.read_text());curves=[]
        for label,directory in spec['run_directories'].items():
            source_plan=json.loads((Path(directory)/'results/post_assembly/psmc/psmc_plan.json').read_text())
            job=next(j for j in source_plan['jobs'] if j['sample']==spec['sample'])
            path=Path(directory)/'results/post_assembly/psmc'/str(job['taxid'])/job['sample']/'report/curves.tsv'
            if path.is_file():
                with path.open() as f:curves.extend(dict(run=label,**r) for r in csv.DictReader(f,delimiter='\t'))
        if curves:
            with (out/'method_curves.tsv').open('w',newline='') as f:
                writer=csv.DictWriter(f,fieldnames=list(curves[0]),delimiter='\t');writer.writeheader();writer.writerows(curves)
            source=json.loads((Path(spec['source'])/'manifest.json').read_text())
            try:subprocess.run([str(Path(source['plot_env'])/'bin/Rscript'),str(REPO/'r_scripts/plot_psmc_method_comparison.R'),str(out)],check=True)
            except (OSError,subprocess.CalledProcessError):issues.append('Method comparison rendering failed')
    md=['# PSMC cohort acceptance checks','',f'Completion: **{"INCOMPLETE" if issues else "PASS"}**.','',
        'The normal main.nf run plans all individuals. Its PSMC reports and combined histories are in `results/post_assembly/psmc/` in the main results/review archive. These additional read checks do not modify inputs or infer biological truth.','',
        '| Individual | Method | Outcome | Callable Mb | Heterozygous SNPs | Read-support evidence |','| --- | --- | --- | ---: | ---: | --- |']
    for r in value['results']:
        called=float(r.get('callable_bp') or 0)/1e6
        link=f"[Report](support_{r['sample']}/README.md)" if (out/('support_'+r['sample'])/'README.md').is_file() else 'Not assessed'
        md.append(f"| {r['sample']} | {r['mode']} | {r['status']} | {called:.2f} | {r.get('heterozygous_snps','')} | {link} |")
    md+=['','- [Focused 20 kb versus cohort reliability](20kb_vs_cohort/README.md)',
        '- [HiFi versus cohort genotypes](hifi_vs_cohort/README.md)','',
        'The focused comparison isolates sites admitted between 20 kb and the cohort threshold. HiFi and assembly routes can have similar aggregate heterozygosity but different spatial patterns and histories; genotype agreement alone does not establish demographic equivalence. A matched-callable-scope fit comparison remains a separate sensitivity analysis if their histories disagree.']
    if (out/'method_histories.png').exists():md+=['','![Method comparison on the pilot individual](method_histories.png)','',
        '[Observed and bootstrap curve evidence](method_curves.tsv). Each method uses its own callable scope; these are sensitivity comparisons, not matched-scope demographic validation.']
    if issues:md+=['','Incomplete items: '+ '; '.join(issues)]
    (out/'README.md').write_text('\n'.join(md)+'\n',encoding='utf-8')
    (out/'validation.json').write_text(json.dumps(dict(status='INCOMPLETE' if issues else 'PASS',issues=issues),indent=2)+'\n')
    with tarfile.open(str(root)+'.psmc-checks.review.tar.gz','w:gz') as t:t.add(out,arcname='psmc_checks')
    if issues:raise ValueError('; '.join(issues))


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('prepare');q.add_argument('--baseline',required=True);q.add_argument('--out',required=True)
    for name in ['audit','collect','environment','context']:
        q=sub.add_parser(name);q.add_argument('--root',required=True)
        if name=='environment':q.add_argument('--sample',required=True)
        if name=='context':q.add_argument('--length-run',required=True);q.add_argument('--hifi-run',required=True)
    a=p.parse_args()
    if a.command=='prepare':prepare(a.baseline,a.out)
    elif a.command=='audit':audit(a.root)
    elif a.command=='collect':collect(a.root)
    elif a.command=='context':context(a.root,a.length_run,a.hifi_run)
    else:print(environment(a.root,a.sample))


if __name__=='__main__':main()
