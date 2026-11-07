"""Run the complete pinned Dipcall route with an explicit assembly preset."""
import argparse
import json
from pathlib import Path
import re
import subprocess


def override(text,preset):
    if preset not in ['asm5','asm10','asm20']:
        raise ValueError('PSMC assembly preset must be asm5, asm10 or asm20')
    rows=re.findall(r'^MM2_OPT=(.*)$',text,re.M)
    if len(rows)!=1:raise ValueError('Expected one generated Dipcall MM2_OPT assignment')
    value,count=re.subn(r'(?<!\S)-xasm5(?!\S)','-x'+preset,rows[0])
    if count!=1:raise ValueError('Cannot isolate standard Dipcall asm5 option')
    return 'MM2_OPT='+value


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--preset',choices=['asm5','asm10','asm20'],required=True)
    p.add_argument('--threads-per-hap',type=int,required=True)
    p.add_argument('--haplotypes',nargs=2,required=True)
    a=p.parse_args()
    if a.threads_per_hap<1:raise ValueError('Positive thread count required')
    with Path('diploid.mak').open('w') as f:
        subprocess.run(['run-dipcall','-a','-t',str(a.threads_per_hap),'diploid','reference.fa',*a.haplotypes],stdout=f,check=True)
    option=override(Path('diploid.mak').read_text(),a.preset)
    command=['make','-j2','-f','diploid.mak',option]
    Path('dipcall_settings.json').write_text(json.dumps(dict(assembly_preset=a.preset,
        make_command=command,threads_per_hap=a.threads_per_hap,
        confidence_thresholds='Unmodified pinned Dipcall defaults'),indent=2)+'\n')
    # The same MM2_OPT is used by both PAF/mask and SAM/VCF targets.
    with Path('dipcall.log').open('w') as f:subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,check=True)


if __name__=='__main__':main()
