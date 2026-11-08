"""Collect length trials and their read-support audits without large intermediates."""
import argparse
import csv
import json
from pathlib import Path
import shutil
import subprocess
import tarfile


def collect(root):
    root=Path(root).resolve();spec=json.loads((root/'comparison.json').read_text());review=root/'review'
    review.mkdir(exist_ok=True);rows=[];curves=[];issues=[]
    for label,directory in spec['run_directories'].items():
        run=Path(directory)
        if not (run/'exit_status.txt').is_file() or (run/'exit_status.txt').read_text().strip()!='0':
            issues.append(label+': workflow incomplete or failed');continue
        plan=json.loads((run/'results/post_assembly/psmc/psmc_plan.json').read_text())
        job=next(j for j in plan['jobs'] if j['sample']==spec['sample'])
        base=run/'results/post_assembly/psmc'/str(job['taxid'])/job['sample']
        q=json.loads((base/'input/input_qc.json').read_text())
        rows.append(dict(run=label,callable_bp=q['callable_bp'],callable_percent=100*q['callable_fraction'],
            heterozygous_snps=q['heterozygous_snps'],heterozygosity_percent=100*q['heterozygosity'],
            retained_callable_bp=q['psmc_retained_callable_bp'],scaffolds=len(q['retained_contigs']),status=q['status']))
        dest=review/'runs'/label
        for name in ['report','fits','reference/reference.json','input/input_qc.json','input/mask_audit.tsv',
                     'input/contig_qc.tsv','input/conversion_audit.tsv','input/provenance.json',
                     'genotypes/dipcall_settings.json','genotypes/dipcall.makefile.diff','genotypes/dipcall.log']:
            src=base/name;target=dest/name
            if not src.exists():continue
            target.parent.mkdir(parents=True,exist_ok=True)
            if src.is_dir():shutil.copytree(src,target,dirs_exist_ok=True)
            else:shutil.copyfile(src,target)
        for name in ['trace.tsv','params.json']:
            if (run/name).is_file():shutil.copyfile(run/name,dest/name)
        if (base/'report/curves.tsv').is_file():
            with (base/'report/curves.tsv').open() as f:
                curves.extend(dict(run=label,**r) for r in csv.DictReader(f,delimiter='\t'))
    for trial in spec['trials']:
        label=trial['label'];audit=root/('reliability_'+label)
        if not (audit/'exit_status.txt').is_file() or (audit/'exit_status.txt').read_text().strip()!='0':
            issues.append(label+': reliability audit incomplete or failed')
        if audit.is_dir():shutil.copytree(audit,review/'reliability'/label,dirs_exist_ok=True)
    fields=['run','callable_bp','callable_percent','heterozygous_snps','heterozygosity_percent','retained_callable_bp','scaffolds','status']
    with (review/'summary.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
    if curves:
        with (review/'curves.tsv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(curves[0]),delimiter='\t');w.writeheader();w.writerows(curves)
    shutil.copyfile(root/'comparison.json',review/'comparison.json')
    if (root/'logs').exists():shutil.copytree(root/'logs',review/'logs',dirs_exist_ok=True)
    status='INCOMPLETE' if issues else 'PASS'
    md=['# PSMC asm20 alignment-length comparison','',f'Individual: **{spec["sample"]}**. Completion: **{status}**.','',
        'The 50 kb control is reused. New 20 kb and 10 kb runs change mask and SAM-filter lengths together; mapper preset, 10 kb ambiguity rule, MAPQ, scope, genotype filters, fitting and bootstrap settings stay fixed.','',
        '| Run | Callable Mb | Scope % | Heterozygosity | Retained Mb | Scaffolds | Fit |','| --- | ---: | ---: | ---: | ---: | ---: | --- |']
    for r in rows:md.append(f"| [{r['run']}](runs/{r['run']}/report/report.md) | {r['callable_bp']/1e6:.2f} | {r['callable_percent']:.2f} | {r['heterozygosity_percent']:.3f}% | {r['retained_callable_bp']/1e6:.2f} | {r['scaffolds']} | {'Produced' if r['status']=='ready' else 'Not fitted'} |")
    md+=['','![Observed histories](length_histories.png)','',
        'Completion PASS verifies completed runs/audits, not biological acceptance. Select using genotype agreement, sampled read support, loss categories and history sensitivity; greater retained sequence alone is insufficient.','']
    for trial in spec['trials']:md.append(f"- [50 kb versus {trial['min_alignment_bp']//1000} kb reliability](reliability/{trial['label']}/README.md)")
    if issues:md+=['','Incomplete items: '+ '; '.join(issues)]
    (review/'README.md').write_text('\n'.join(md)+'\n',encoding='utf-8')
    if curves:
        source=json.loads((Path(spec['source'])/'manifest.json').read_text())
        try:subprocess.run([str(Path(source['plot_env'])/'bin/Rscript'),str(Path(spec['repository'])/'r_scripts/plot_psmc_length_comparison.R'),str(review)],check=True)
        except (OSError,subprocess.CalledProcessError):issues.append('Summary renderer failed');status='INCOMPLETE'
    (review/'validation.json').write_text(json.dumps(dict(status=status,issues=issues),indent=2)+'\n')
    if issues:(review/'README.md').write_text((review/'README.md').read_text().replace('Completion: **PASS**','Completion: **INCOMPLETE**'))
    with tarfile.open(str(root)+'.review.tar.gz','w:gz') as t:t.add(review,arcname='review')
    if issues:raise ValueError('; '.join(issues))
    print('Combined review archive:',str(root)+'.review.tar.gz')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--comparison',required=True)
    collect(p.parse_args().comparison)
