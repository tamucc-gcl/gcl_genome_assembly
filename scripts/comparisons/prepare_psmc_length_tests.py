"""Prepare bounded asm20 alignment-length trials against a retained 50 kb control."""
import argparse
import json
from pathlib import Path
from psmc_reliability import hifi_bam


def prepare(source,out,hifi):
    source=Path(source).resolve();out=Path(out).resolve();hifi=Path(hifi).resolve()
    if out.exists():raise ValueError('Use a new length-test output directory')
    comparison=json.loads((source/'comparison.json').read_text());sample=comparison['sample']
    control=source/'asm20'
    if (control/'exit_status.txt').read_text().strip()!='0':raise ValueError('asm20 control did not complete successfully')
    params=json.loads((control/'params.json').read_text())
    if params['psmc_assembly_preset']!='asm20' or params.get('psmc_assembly_min_alignment_bp',50000)!=50000:
        raise ValueError('Expected asm20 / 50000 bp control')
    plan=json.loads((control/'results/post_assembly/psmc/psmc_plan.json').read_text())
    jobs=[j for j in plan['jobs'] if j['sample']==sample and j['mode']=='assemblies']
    if len(jobs)!=1:raise ValueError('Expected one matching assembly-route individual')
    job=jobs[0]
    for path in [params['psmc_final_eligibility'],params['psmc_manifest'],*job['assemblies']]:
        if not Path(path).is_file():raise ValueError('Required retained input unavailable: '+str(path))
    bam=hifi_bam(hifi,sample,job['reference']['sha256'])
    settings=control/'results/post_assembly/psmc'/str(job['taxid'])/sample/'genotypes/dipcall_settings.json'
    if json.loads(settings.read_text())['assembly_preset']!='asm20':raise ValueError('Control caller provenance disagrees')
    out.mkdir(parents=True);trials=[];directories={'asm20_50kb':str(control)}
    for minimum in [20000,10000]:
        label=f'asm20_{minimum//1000}kb';dest=out/label;dest.mkdir();trial=dict(params)
        trial.update(psmc_assembly_min_alignment_bp=minimum,outdir=str(dest/'results'))
        (dest/'params.json').write_text(json.dumps(trial,indent=2)+'\n')
        trials.append(dict(label=label,min_alignment_bp=minimum));directories[label]=str(dest)
    spec=dict(schema_version=1,mode='alignment_lengths',repository=str(Path(__file__).resolve().parents[2]),
        sample=sample,source=comparison['source'],control_source=str(source),trials=trials,
        run_directories=directories,hifi_run=str(hifi),hifi_bam=str(bam),baseline_label='asm20_50kb')
    (out/'comparison.json').write_text(json.dumps(spec,indent=2)+'\n')
    print(f'Prepared two full asm20 length trials for {sample}; existing 50 kb control reused.')
    return spec


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['source','out','hifi-run']:p.add_argument('--'+key,required=True)
    a=p.parse_args();prepare(a.source,a.out,a.hifi_run)


if __name__=='__main__':main()
