"""Aggregate explicit current-run PSMC reports; no genotype or fitting work."""
import argparse
import csv
import json
from pathlib import Path


def index(plan_path,result_directories,out):
    plan=json.loads(Path(plan_path).read_text(encoding='utf-8'));rows=[];by_id={};curves=[]
    for directory in result_directories:
        p=Path(directory)/'summary.json';r=json.loads(p.read_text(encoding='utf-8'))
        if r['sample'] in by_id:raise ValueError('Duplicate PSMC result')
        by_id[r['sample']]=r
        if r['status']=='completed':
            with (Path(directory)/'curves.tsv').open(encoding='utf-8') as f:
                for c in csv.DictReader(f,delimiter='\t'):
                    c.update(sample=r['sample'],taxid=r['taxid'],reference=r['reference'],mode=r['mode'])
                    curves.append(c)
    expected={j['sample'] for j in plan['jobs']}
    if set(by_id)!=expected:raise ValueError('Missing/unexpected PSMC individual reports')
    for r in plan['skipped']:
        rows.append(dict(sample=r['sample'],taxid=r['taxid'],mode='not_run',status='skipped',reason=r['reason'],reference='',callable_bp='',heterozygous_snps='',scaling='',report=''))
    for sample,r in sorted(by_id.items()):
        rows.append({k:r.get(k,'') for k in ['sample','taxid','mode','status','reason','reference','callable_bp','heterozygous_snps',
            'retained_callable_bp','retained_scaffolds','retained_contigs','missing_bin_fraction','scaling','report']})
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    fields=['sample','taxid','mode','status','reason','reference','callable_bp','heterozygous_snps',
        'retained_callable_bp','retained_scaffolds','retained_contigs','missing_bin_fraction','scaling','report']
    with (out/'summary.tsv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
    counts={status:sum(r['status']==status for r in rows) for status in ['completed','not_fitted','skipped']}
    md=['# PSMC demographic history','',
        f"{counts['completed']} completed; {counts['not_fitted']} not fitted; {counts['skipped']} skipped.",
        '' if rows else 'No finalized individuals were available for PSMC; no histories were inferred.',
        '', 'One diploid individual per analysis. Long-scaffold scope may include sex-linked sequence.','',
        '| Individual | Own reference | Input mode | Result | Callable bp before conversion | Callable bp retained for fitting | Retained scaffolds | Scaling | Evidence |',
        '| --- | --- | --- | --- | ---: | ---: | ---: | --- | --- |']
    for r in rows:
        link=f"[Report]({r['taxid']}/{r['sample']}/report/report.md)" if r['report'] else str(r['reason']).replace('|','\\|')
        md.append(f"| {r['sample']} | {r['reference']} | {r['mode']} | {r['status']} | {r['callable_bp']} | {r.get('retained_callable_bp','')} | {r.get('retained_scaffolds','')} | {r['scaling']} | {link} |")
    if curves:
        with (out/'combined_curves.tsv').open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=list(curves[0]),delimiter='\t');w.writeheader();w.writerows(curves)
        md+=['','## Individuals together','','![Observed histories for all fitted individuals](combined_history.png)','',
            '[PDF](combined_history.pdf) · [Combined curve evidence](combined_curves.tsv)','',
            'Observed fits use logarithmic time and size axes so extreme curves do not hide the others. Species and scaling groups are shown separately; each individual uses its own reference. These are exploratory results; check fitted scope before interpretation.','',
            '![Individual histories and block bootstraps](individual_histories.png)','',
            '[Panel figure PDF](individual_histories.pdf)','',
            'Companion panels also use logarithmic axes, with separate panel ranges to keep bootstrap extremes visible.','']
    md+=['','[Structured summary](summary.tsv) · [Input plan](psmc_plan.json)','',
         'Curves are model-based histories. Assembly phasing, callability, heterozygote detection and selected time intervals affect inference. Bootstrap curves measure block-resampling variation, not uncertainty in mutation rate, generation time or model assumptions.']
    (out/'README.md').write_text('\n'.join(md)+'\n',encoding='utf-8')
    (out/'psmc_plan.json').write_text(json.dumps(plan,indent=2)+'\n',encoding='utf-8')



def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=['index'])
    p.add_argument('--plan',required=True)
    p.add_argument('--directories',nargs='*',default=[])
    p.add_argument('--out',required=True)
    a=p.parse_args();index(a.plan,a.directories,a.out)


if __name__=='__main__':main()
