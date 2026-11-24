#!/usr/bin/env python3
"""Check real FCS containers using NCBI's published adaptor fixture."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import urllib.request

# Nextflow stages these helpers together, often as symlinks into node scratch.
# Python can resolve this script back to tests/, so explicitly include the task directory.
sys.path.insert(0,str(Path.cwd()))
from contamination_review import records
from contamination_cohort_report import adaptor_rows

SOURCE='https://zenodo.org/records/10932013/files/FCS_combo_test.fa'


def fasta(path):
    return {h.split()[0]:s for h,s in records(path)}


def check(source,report,clean,rescreen,out):
    original=fasta(source)
    actual=fasta(clean)
    actions=adaptor_rows(report)
    expected_ids={f'seq_{i:05d}' for i in range(1,11)}|{'seq_00018'}
    if len(actions)!=11 or {x['accession'] for x in actions}!=expected_ids:
        raise ValueError('NCBI fixture did not produce the 11 published adaptor calls')
    expected=dict(original)
    for row in actions:
        name=row['accession'];sequence=expected.pop(name)
        if name=='seq_00018':
            if row['action']!='ACTION_EXCLUDE' or len(sequence)!=522:
                raise ValueError('Expected adaptor-dominated sequence exclusion')
            continue
        interval=(1,58) if int(name.split('_')[1])<=5 else (100001,100058)
        if row['action']!='ACTION_TRIM' or row['range']!=f'{interval[0]}..{interval[1]}':
            raise ValueError('Unexpected planted adaptor location: '+name)
        lo,hi=interval
        if lo>1:expected[f'{name}~1..{lo-1}']=sequence[:lo-1]
        if hi<len(sequence):expected[f'{name}~{hi+1}..{len(sequence)}']=sequence[hi:]
    if actual!=expected:
        raise ValueError('Cleaned records differ from exact terminal trim/internal split expectations')
    if adaptor_rows(rescreen):
        raise ValueError('Adaptor findings remain after cleaning the fixture')
    Path(out).write_text(json.dumps(dict(status='PASS',source_url=SOURCE,
        source_sha256=hashlib.sha256(Path(source).read_bytes()).hexdigest(),
        input_sequences=len(original),output_sequences=len(actual),terminal_trims=5,internal_splits=5,
        excluded_sequences=1,removed_bp=sum(map(len,original.values()))-sum(map(len,actual.values())),
        repeat_adaptor_actions=0),indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['fetch','check'])
    for key in ('source','report','clean','rescreen','out'):parser.add_argument('--'+key,type=Path)
    args=parser.parse_args()
    if args.command=='fetch':
        with urllib.request.urlopen(SOURCE,timeout=120) as response,args.out.open('wb') as handle:
            while chunk:=response.read(1024*1024):handle.write(chunk)
        if 'seq_00018' not in fasta(args.out):raise ValueError('Downloaded input is not the NCBI fixture')
    else:check(args.source,args.report,args.clean,args.rescreen,args.out)
