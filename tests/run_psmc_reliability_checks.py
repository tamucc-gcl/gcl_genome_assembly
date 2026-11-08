"""Exercise the reliability report using completed native calls and HiFi mapping."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['asm10','asm20','mask-source','hifi','out','plot-env']:p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();repo=Path(__file__).resolve().parents[1];out=a.out.resolve();out.mkdir(parents=True,exist_ok=False)
    comparison=out/'comparison';comparison.mkdir()
    for preset,source in [('asm10',a.asm10),('asm20',a.asm20)]:
        dest=comparison/preset;dest.mkdir();(dest/'results').symlink_to(source.resolve()/'results');(dest/'exit_status.txt').write_text('0\n')
    (comparison/'comparison.json').write_text(json.dumps(dict(sample='assembly_ind',source=str(a.mask_source.resolve()))))
    hifi=out/'hifi';hifi.mkdir();(hifi/'results').symlink_to(a.hifi.resolve()/'results')
    bams=[p for p in (a.hifi/'work').glob('*/*/hifi.bam') if not p.is_symlink()];assert len(bams)==1,bams
    with (hifi/'trace.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['name','status','exit','workdir'],delimiter='\t');w.writeheader()
        w.writerow(dict(name='PSMC_ANALYSIS:PSMC_MAP_HIFI (assembly_ind)',status='COMPLETED',exit='0',workdir=str(bams[0].resolve().parent)))
    report=out/'report'
    subprocess.run([sys.executable,str(repo/'scripts/comparisons/psmc_reliability.py'),'--comparison',str(comparison),'--out',str(report),
        '--hifi-run',str(hifi),'--sample-sites','20'],check=True)
    subprocess.run([str(a.plot_env/'bin/Rscript'),str(repo/'r_scripts/plot_psmc_reliability.R'),str(report)],check=True)
    result=json.loads((report/'summary.json').read_text())
    assert result['shared_genotype_agreement']==1
    assert result['shared_het_genotype_agreement']==1
    with (report/'sampled_read_support.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    assert len(rows)==20 and all(r['left_support']==r['right_support'] for r in rows)
    assert any(r['left_support']=='supported' for r in rows)
    assert (report/'reliability.png').stat().st_size>1000
    print('PASS: complete consensus comparison, exact callable reproduction, mapped read support and rendered R report')


if __name__=='__main__':main()
