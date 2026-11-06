"""Aggregate explicit current-run PSMC reports; no genotype or fitting work."""
import argparse
import csv
import json
from pathlib import Path


def index(plan_path,result_directories,out):
    plan=json.loads(Path(plan_path).read_text());rows=[];by_id={}
    for directory in result_directories:
        p=Path(directory)/'summary.json';r=json.loads(p.read_text())
        if r['sample'] in by_id:raise ValueError('Duplicate PSMC result')
        by_id[r['sample']]=r
    expected={j['sample'] for j in plan['jobs']}
    if set(by_id)!=expected:raise ValueError('Missing/unexpected PSMC individual reports')
    for r in plan['skipped']:
        rows.append(dict(sample=r['sample'],taxid=r['taxid'],mode='not_run',status='skipped',reason=r['reason'],reference='',callable_bp='',heterozygous_snps='',scaling='',report=''))
    for sample,r in sorted(by_id.items()):
        rows.append({k:r.get(k,'') for k in ['sample','taxid','mode','status','reason','reference','callable_bp','heterozygous_snps','scaling','report']})
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    fields=['sample','taxid','mode','status','reason','reference','callable_bp','heterozygous_snps','scaling','report']
    with (out/'summary.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
    counts={status:sum(r['status']==status for r in rows) for status in ['completed','not_fitted','skipped']}
    md=['# PSMC demographic history','',
        f"{counts['completed']} completed; {counts['not_fitted']} not fitted; {counts['skipped']} skipped.",
        '' if rows else 'No finalized individuals were available for PSMC; no histories were inferred.',
        '', 'One diploid individual per analysis. Long-scaffold scope may include sex-linked sequence.','',
        '| Individual | Input mode | Result | Callable bp | Heterozygous SNPs | Scaling | Evidence |',
        '| --- | --- | --- | ---: | ---: | --- | --- |']
    for r in rows:
        link=f"[Report]({r['taxid']}/{r['sample']}/report/report.md)" if r['report'] else str(r['reason']).replace('|','\\|')
        md.append(f"| {r['sample']} | {r['mode']} | {r['status']} | {r['callable_bp']} | {r['heterozygous_snps']} | {r['scaling']} | {link} |")
    md+=['','[Structured summary](summary.tsv) Â· [Input plan](psmc_plan.json)','',
         'Curves are model-based histories. Assembly phasing, callability, heterozygote detection and selected time intervals affect inference. Bootstrap curves measure block-resampling variation, not uncertainty in mutation rate, generation time or model assumptions.']
    (out/'README.md').write_text('\n'.join(md)+'\n')
    (out/'psmc_plan.json').write_text(json.dumps(plan,indent=2)+'\n')



def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=['index'])
    p.add_argument('--plan',required=True)
    p.add_argument('--directories',nargs='*',default=[])
    p.add_argument('--out',required=True)
    a=p.parse_args();index(a.plan,a.directories,a.out)


if __name__=='__main__':main()
