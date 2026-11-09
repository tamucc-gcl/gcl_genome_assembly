"""Native post-run cohort QC using actual retained calls and mapped reads."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['native','paired-hifi','length-run','plot-env','out']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();repo=Path(__file__).resolve().parents[1];root=a.out.resolve();root.mkdir(parents=True,exist_ok=False)
    shutil.copytree(a.native/'results',root/'results');(root/'exit_status.txt').write_text('0')
    psmc=root/'results/post_assembly/psmc';plan=json.loads((psmc/'psmc_plan.json').read_text())
    names=sorted([j['sample'] for j in plan['jobs']]+[j['sample'] for j in plan['skipped']])
    (root/'psmc_expected.json').write_text(json.dumps(dict(samples=names)))
    provenance=root/'results/provenance';provenance.mkdir(exist_ok=True)
    (provenance/'run_provenance.json').write_text(json.dumps(dict(context=dict(parameters=dict(psmc_assembly_preset='asm20',psmc_assembly_min_alignment_bp=10000)))))
    final=root/'results/bam/hifi/final';final.mkdir(parents=True)
    for job in plan['jobs']:
        native=a.paired_hifi if job['sample']=='assembly_ind' else a.native
        bams=[b for b in (native/'work').glob('*/*/hifi.bam') if not b.is_symlink()];assert len(bams)==1
        digest=hashlib.sha256((bams[0].parent/'reference.fa').read_bytes()).hexdigest();assert digest==job['reference']['sha256']
        dest=final/(job['reference']['id']+'.sorted.bam');dest.symlink_to(bams[0]);Path(str(dest)+'.bai').symlink_to(Path(str(bams[0])+'.bai'))
        (final/(job['reference']['id']+'.reference.sha256')).write_text(digest+'  verified_fixture_reference\n')
    sys.path.insert(0,str(repo/'scripts/comparisons'))
    from psmc_cohort import audit,context,collect
    from psmc_input_support import run
    checked=audit(root);assert checked['status']=='PASS' and checked['planned']==2 and len(checked['skipped'])==1
    check=root/'psmc_checks';check.mkdir()
    for job in plan['jobs']:
        dest=check/('support_'+job['sample']);support=run(root/'results',job['sample'],dest,20)
        assert any(g['supported']>0 for g in support['groups'])
        assert next(g for g in support['groups'] if g['group']=='heterozygous')['sampled']==20
        (dest/'exit_status.txt').write_text('0')
    hifi=root/'hifi_source';hifi.mkdir();(hifi/'results').symlink_to(a.paired_hifi.resolve()/'results')
    (hifi/'exit_status.txt').write_text('0')
    bam=next(b for b in (a.paired_hifi/'work').glob('*/*/hifi.bam') if not b.is_symlink())
    with (hifi/'trace.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['name','status','exit','workdir'],delimiter='\t');w.writeheader()
        w.writerow(dict(name='PSMC_ANALYSIS:PSMC_MAP_HIFI (assembly_ind)',status='COMPLETED',exit='0',workdir=str(bam.parent)))
    context(root,a.length_run,hifi)
    for left,label in [('control_20kb','20kb_vs_cohort'),('hifi','hifi_vs_cohort')]:
        dest=check/label
        subprocess.run([sys.executable,str(repo/'scripts/comparisons/psmc_reliability.py'),'--comparison',str(check),'--out',str(dest),
            '--left-run',left,'--right-run','cohort','--hifi-run',str(hifi),'--sample-sites','20'],check=True)
        subprocess.run([str(a.plot_env/'bin/Rscript'),str(repo/'r_scripts/plot_psmc_reliability.R'),str(dest)],check=True)
        (dest/'exit_status.txt').write_text('0')
    collect(root)
    archive=Path(str(root)+'.psmc-checks.review.tar.gz');assert archive.is_file()
    with tarfile.open(archive) as t:assert not any(m.name.endswith(('.bam','.fa.gz','.vcf.gz')) for m in t)
    # Missing any individual's support evidence is explicitly incomplete.
    status=check/'support_assembly_ind/exit_status.txt';status.write_text('1')
    try:collect(root)
    except ValueError:pass
    else:raise AssertionError('Failed support audit cannot pass')
    status.write_text('0');collect(root)
    print('PASS: all planned individuals, real final-BAM support, both focused comparisons, R rendering and compact incomplete/success reviews')


if __name__=='__main__':main()
