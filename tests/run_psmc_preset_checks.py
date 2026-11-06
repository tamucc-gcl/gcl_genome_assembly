"""Run real asm10/asm20 mask-only trials against a completed native mask fixture."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();repo=Path(__file__).resolve().parents[1]
    source=json.loads((a.source/'manifest.json').read_text());sample=source['tasks'][0]['sample']
    env=dict(os.environ,PATH=source['tasks'][0]['env']+'/bin:'+os.environ['PATH'])
    helper=repo/'scripts/comparisons/psmc_mask_investigation.py';out=a.out.resolve()
    before=hashlib.sha256((a.source/'manifest.json').read_bytes()).hexdigest()
    subprocess.run(['python',str(helper),'prepare-presets','--source',str(a.source),'--out',str(out),'--sample',sample],check=True,env=env)
    manifest=out/'manifest.json';data=json.loads(manifest.read_text())
    assert len(data['tasks'])==4 and len(data['controls'])==2
    spool=out/'spool';spool.mkdir();script=spool/'slurm_script'
    shutil.copyfile(repo/'scripts/comparisons/run_psmc_mask_investigation.sbatch',script)
    running=[subprocess.Popen(['bash',str(script),str(manifest),'task',str(i)],env=dict(env,SLURM_CPUS_PER_TASK='2')) for i in range(4)]
    assert all(p.wait()==0 for p in running)
    subprocess.run(['bash',str(script),str(manifest),'report-presets'],check=True,env=env)
    status=json.loads((out/'reports/status.json').read_text());assert status['status']=='PASS',status
    for j in data['tasks']:
        root=out/'tasks'/f"{sample}.{j['preset']}.{j['haplotype']}"
        log=(root/'raw'/f"diploid.{j['haplotype']}.paf.gz.log").read_text()
        assert '-x'+j['preset'] in log,log
        assert not any(p.suffix in ['.sam','.bam'] or p.name.endswith('.vcf.gz') for p in (root/'raw').iterdir())
    assert (out/'reports/preset_comparison.png').stat().st_size>1000
    assert hashlib.sha256((a.source/'manifest.json').read_bytes()).hexdigest()==before
    with tarfile.open(str(out)+'.review.tar.gz') as t:
        assert not any('/raw/' in m.name or m.name.endswith(('.paf.gz','.fa','.bam')) for m in t)
    missing=out/'tasks'/f'{sample}.asm20.hap2/audit.json';saved=missing.read_bytes();missing.unlink()
    result=subprocess.run(['bash',str(script),str(manifest),'report-presets'],env=env,capture_output=True,text=True)
    assert result.returncode and json.loads((out/'reports/status.json').read_text())['status']=='INCOMPLETE'
    missing.write_bytes(saved)
    subprocess.run(['bash',str(script),str(manifest),'report-presets'],check=True,env=env)
    print('PASS: four native preset trials, preserved controls, actual commands, mask validation, R plots, failure archive')


if __name__=='__main__':main()
