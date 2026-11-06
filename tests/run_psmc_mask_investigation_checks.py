"""Native removable-helper checks using an already completed PSMC fixture run."""
import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['native-root','analysis-env','plot-env','out']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();repo=Path(__file__).resolve().parents[1]
    helper=repo/'scripts/comparisons/psmc_mask_investigation.py'
    root=a.out.resolve();root.mkdir(parents=True,exist_ok=False)
    baseline=root/'baseline';baseline.mkdir()
    shutil.copytree(a.native_root/'results',baseline/'results')
    plan_path=baseline/'results/post_assembly/psmc/psmc_plan.json'
    plan=json.loads(plan_path.read_text());job=next(j for j in plan['jobs'] if j['mode']=='assemblies')
    job['reads']=str((a.native_root/'fixtures/complete.fastq').resolve())
    plan_path.write_text(json.dumps(plan))
    pafs=list((a.native_root/'work').glob('*/*/diploid.hap1.paf.gz'))
    assert len(pafs)==1,pafs
    cache=root/'retained';cache.mkdir();(cache/'diploid.hap1.paf.gz').symlink_to(pafs[0])
    records=[dict(name=f"PSMC_ANALYSIS:PSMC_ASSEMBLY_CALL ({job['sample']})",workdir=str(cache),status='COMPLETED',exit='0'),
             dict(name='PSMC_ANALYSIS:PSMC_INDEX',workdir=str(root/'plot_task'),status='COMPLETED',exit='0')]
    with (baseline/'trace.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['name','workdir','status','exit'],delimiter='\t');w.writeheader();w.writerows(records)
    (baseline/'execution_environments.json').write_text(json.dumps(dict(tasks=[
        dict(workdir=r['workdir'],conda_prefixes=[str(a.analysis_env if i==0 else a.plot_env)]) for i,r in enumerate(records)])))
    env=dict(os.environ,PATH=str(a.analysis_env/'bin')+':'+os.environ['PATH'])
    python=str(a.analysis_env/'bin/python');out=root/'investigation'
    subprocess.run([python,str(helper),'prepare','--run',str(baseline),'--out',str(out),'--hifi-sample',job['sample']],check=True,env=env)
    manifest=out/'manifest.json';data=json.loads(manifest.read_text())
    assert data['tasks'][0]['cached_paf'] and data['tasks'][1]['cached_paf'] is None
    # Running a copied job script models Slurm's spool-directory $0 behavior.
    spool=root/'spool';spool.mkdir();script=spool/'slurm_script'
    shutil.copyfile(repo/'scripts/comparisons/run_psmc_mask_investigation.sbatch',script)
    for i in [0,1]:subprocess.run(['bash',str(script),str(manifest),'task',str(i)],check=True,env=env)
    subprocess.run([python,str(helper),'report','--manifest',str(manifest)],check=True,env=env)
    status=json.loads((out/'reports/status.json').read_text());assert status['status']=='PASS',status
    assert {s['alignment_source'] for s in status['summary']}=={'cached','regenerated_PAF_only'}
    assert (out/'reports/mask_reasons.png').stat().st_size>0
    archive=Path(str(out)+'.review.tar.gz');assert archive.is_file()
    with tarfile.open(archive) as t:assert not any(m.name.endswith(('.paf.gz','.fa','.fasta','.bam')) for m in t)
    second=out/'tasks'/f"{job['sample']}.hap2/audit.json";saved=second.with_suffix('.saved')
    second.rename(saved)
    result=subprocess.run([python,str(helper),'report','--manifest',str(manifest)],env=env,capture_output=True,text=True)
    assert result.returncode!=0 and json.loads((out/'reports/status.json').read_text())['status']=='INCOMPLETE'
    assert archive.is_file() # Failure evidence is packaged, not represented as success.
    saved.rename(second)
    subprocess.run([python,str(helper),'report','--manifest',str(manifest)],check=True,env=env)
    print('PASS: cached and regenerated PAFs, native mask reproduction, Slurm spool path, R figures, missing-task failure and compact archive')


if __name__=='__main__':main()
