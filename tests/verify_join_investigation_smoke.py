#!/usr/bin/env python3
"""Assert real-tool smoke outputs, including report links and signal direction."""
import argparse
import csv
import json
from pathlib import Path
import re


def rows(path):
    with path.open() as src:return list(csv.DictReader(src,delimiter='\t'))


def main():
    p=argparse.ArgumentParser();p.add_argument('review',type=Path);a=p.parse_args();root=a.review
    sites=rows(root/'sites.tsv')
    supported=next(r for r in sites if r['assembly']=='A_hap1' and r['site']=='C07T')
    unsupported=next(r for r in sites if r['assembly']=='A_hap2' and r['site']=='C07T')
    assert int(supported['unique_arm_bridges'])==12,supported
    assert int(unsupported['original_mapq30'])==int(unsupported['competitive_credible'])==0,unsupported
    for assembly in ('A_hap1','A_hap2'):
        hic=json.loads((root/assembly/'hic/summary.json').read_text())
        assert hic['replayed_measurements']>0
        assert {r['resolution'] for r in hic['maps']}=={10000,25000,50000,100000}
        assert {r['library'] for r in hic['maps']}=={'single','POOLED'}
    profile=rows(root/'A_hap2/hic/contact_profile.tsv')
    zero=next(r for r in profile if r['point']=='600000' and r['flank_bp']=='300000' and r['library']=='single')
    assert int(zero['cross'])==0 and min(int(zero['left_within']),int(zero['right_within']))>=100,zero
    for action in rows(root/'manual-cut-worksheet.tsv'):
        assert action['selected']=='NO' and action['cut_bp']==''
    markdown=(root/'README.md').read_text()
    links=re.findall(r'\]\(([^)]+)\)',markdown)
    for link in links:
        if not link.startswith('#'):assert (root/link).exists(),link
    assert not any(p.suffix in ('.sam','.bam','.fa') for p in root.rglob('*'))
    print('PASS: supported and unsupported junctions, observable Hi-C zero, all resolutions, single library, replay, report links, and unselected cuts')


if __name__=='__main__':main()
