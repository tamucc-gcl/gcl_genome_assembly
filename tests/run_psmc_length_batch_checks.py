"""Check preparation, paired reliability and collection using native PSMC outputs."""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['control','trial20','trial10','mask-source','hifi','plot-env','out']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();repo=Path(__file__).resolve().parents[1];root=a.out.resolve();root.mkdir(parents=True,exist_ok=False)
    hifi=root/'hifi';hifi.mkdir();(hifi/'results').symlink_to(a.hifi.resolve()/'results')
    bams=[p for p in (a.hifi/'work').glob('*/*/hifi.bam') if not p.is_symlink()];assert len(bams)==1
    with (hifi/'trace.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['name','status','exit','workdir'],delimiter='\t');w.writeheader()
        w.writerow(dict(name='PSMC_ANALYSIS:PSMC_MAP_HIFI (assembly_ind)',status='COMPLETED',exit='0',workdir=str(bams[0].parent)))
    source=root/'full';control=source/'asm20';control.mkdir(parents=True)
    (control/'results').symlink_to(a.control.resolve()/'results');(control/'exit_status.txt').write_text('0')
    (source/'comparison.json').write_text(json.dumps(dict(sample='assembly_ind',source=str(a.mask_source.resolve()))))
    params=dict(psmc_assembly_preset='asm20',psmc_assembly_min_alignment_bp=50000,
        psmc_final_eligibility=str(a.control.resolve()/'fixtures/eligibility.json'),psmc_manifest=str(a.control.resolve()/'fixtures/manifest.tsv'),
        psmc_bootstraps=2,psmc_pattern='2+2+2',outdir=str(control/'results'))
    (control/'params.json').write_text(json.dumps(params))
    sys.path.insert(0,str(repo/'scripts/comparisons'))
    from prepare_psmc_length_tests import prepare
    from psmc_length_summary import collect
    out=root/'lengths';spec=prepare(source,out,hifi)
    for label,native in [('asm20_20kb',a.trial20),('asm20_10kb',a.trial10)]:
        dest=out/label;(dest/'results').symlink_to(native.resolve()/'results');(dest/'exit_status.txt').write_text('0')
        caller=json.loads((native/'results/post_assembly/psmc/1/assembly_ind/genotypes/dipcall_settings.json').read_text())
        expected=20000 if label.endswith('20kb') else 10000
        assert caller['confidence_thresholds']['mask_min_alignment_bp']==expected
        assert caller['confidence_thresholds']['sam_filter_min_alignment_bp']==expected
        audit=out/('reliability_'+label)
        subprocess.run([sys.executable,str(repo/'scripts/comparisons/psmc_reliability.py'),'--comparison',str(out),'--out',str(audit),
            '--left-run','asm20_50kb','--right-run',label,'--hifi-run',str(hifi),'--sample-sites','20'],check=True)
        subprocess.run([str(a.plot_env/'bin/Rscript'),str(repo/'r_scripts/plot_psmc_reliability.R'),str(audit)],check=True)
        (audit/'exit_status.txt').write_text('0')
        summary=json.loads((audit/'summary.json').read_text());assert summary['labels']['right']==label
        assert summary['shared_genotype_agreement']==1
    collect(out)
    archive=Path(str(out)+'.review.tar.gz');assert archive.is_file()
    with tarfile.open(archive) as t:
        assert not any(m.name.endswith(('.fa','.fa.gz','.bam','.vcf.gz','.psmcfa')) for m in t)
    assert json.loads((out/'review/validation.json').read_text())['status']=='PASS'
    failed=out/'reliability_asm20_10kb/exit_status.txt';failed.write_text('1')
    try:collect(out)
    except ValueError:pass
    else:raise AssertionError('Failed audit must not produce a completion PASS')
    assert archive.is_file() and json.loads((out/'review/validation.json').read_text())['status']=='INCOMPLETE'
    failed.write_text('0');collect(out)
    # Capture submission dependencies without contacting a scheduler. Scientific
    # outputs above are real native tool results, not scheduler stand-ins.
    fake=root/'scheduler';fake.mkdir();capture=root/'submitted.json';capture.write_text('[]')
    sbatch=fake/'sbatch';sbatch.write_text('#!'+sys.executable+'\n'+
        'import json,os,sys\nfrom pathlib import Path\np=Path(os.environ["PSMC_SCHEDULER_CAPTURE"])\n'+
        'rows=json.loads(p.read_text());rows.append(sys.argv[1:]);p.write_text(json.dumps(rows));print(len(rows))\n');sbatch.chmod(0o755)
    env=dict(os.environ,PATH=str(fake)+':'+os.environ['PATH'],PSMC_SCHEDULER_CAPTURE=str(capture))
    subprocess.run(['bash',str(repo/'scripts/comparisons/submit_psmc_length_tests.sh'),str(source),str(root/'submission'),str(hifi)],env=env,check=True)
    submitted=json.loads(capture.read_text());assert len(submitted)==5
    assert '--dependency=afterany:1' in submitted[1]
    assert '--dependency=afterany:3' in submitted[3]
    assert '--dependency=afterany:1:2:3:4' in submitted[4]
    assert not any(x.startswith('--dependency') for i in [0,2] for x in submitted[i])
    print('PASS: threshold preparation, both native workflow outputs, paired read-support audits, R rendering, compact archive and failed-audit reporting')


if __name__=='__main__':main()
