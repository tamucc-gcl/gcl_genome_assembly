"""Run the complete pinned Dipcall route with an explicit assembly preset."""
import argparse
import difflib
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


def alignment_threshold(text,minimum):
    if not isinstance(minimum,int) or minimum<10000:
        raise ValueError('Minimum assembly alignment length must be an integer >=10000 bp')
    # Change only this run's generated recipes. Installed package files and the
    # 10 kb coverage/ambiguity threshold and MAPQ 5 remain untouched.
    text,sam=re.subn(r'dipcall-aux\.js samflt (?=\$<)',f'dipcall-aux.js samflt -L {minimum} ',text)
    text,paf=re.subn(r'paftools\.js call (?=- 2>)',f'paftools.js call -L {minimum} ',text)
    if sam!=2 or paf!=2:
        raise ValueError('Expected exactly two SAM filters and two PAF mask calls in generated Dipcall makefile')
    return text


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--preset',choices=['asm5','asm10','asm20'],required=True)
    p.add_argument('--threads-per-hap',type=int,required=True)
    p.add_argument('--haplotypes',nargs=2,required=True)
    p.add_argument('--min-alignment-bp',type=int,default=50000)
    a=p.parse_args()
    if a.threads_per_hap<1:raise ValueError('Positive thread count required')
    with Path('diploid.mak').open('w') as f:
        subprocess.run(['run-dipcall','-a','-t',str(a.threads_per_hap),'diploid','reference.fa',*a.haplotypes],stdout=f,check=True)
    original=Path('diploid.mak').read_text()
    updated=alignment_threshold(original,a.min_alignment_bp)
    Path('diploid.mak').write_text(updated)
    Path('dipcall.makefile.diff').write_text(''.join(difflib.unified_diff(original.splitlines(True),updated.splitlines(True),fromfile='generated upstream recipes',tofile='explicit alignment-length recipes')))
    option=override(updated,a.preset)
    command=['make','-j2','-f','diploid.mak',option]
    Path('dipcall_settings.json').write_text(json.dumps(dict(assembly_preset=a.preset,
        make_command=command,threads_per_hap=a.threads_per_hap,
        confidence_thresholds=dict(mask_min_alignment_bp=a.min_alignment_bp,sam_filter_min_alignment_bp=a.min_alignment_bp,
            coverage_min_alignment_bp=10000,min_mapq=5,overlap_rule='Unmodified pinned Dipcall rule')),
        indent=2)+'\n')
    # The same MM2_OPT is used by both PAF/mask and SAM/VCF targets.
    with Path('dipcall.log').open('w') as f:subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,check=True)


if __name__=='__main__':main()
