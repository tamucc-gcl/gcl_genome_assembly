"""Validate actual corrected sequences and reassigned maps; emit compact reviewable QC."""
import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from break_chimeras import fasta


def n50(sequences):
    total=sum(map(len,sequences.values()));running=0
    for length in sorted(map(len,sequences.values()),reverse=True):
        running+=length
        if running>=total/2:return length
    return 0


def main():
    p=argparse.ArgumentParser()
    for option in ('original','corrected','names','assembly','scenario','actions','out'):p.add_argument('--'+option,required=True)
    a=p.parse_args()
    original,corrected=fasta(a.original),fasta(a.corrected)
    with open(a.names) as handle:names=list(csv.DictReader(handle,delimiter='\t'))
    if {r['old_name'] for r in names}!=set(corrected) or len(names)!=len(corrected):raise ValueError('Reassigned name map differs from FASTA records')
    if len({r['new_name'] for r in names})!=len(names):raise ValueError('Duplicate assigned names')
    reconstruction={};pieces=[];accounted=set()
    for parent,sequence in original.items():
        if parent in corrected:
            reconstruction[parent]=corrected[parent];accounted.add(parent);continue
        fragments=[]
        for child,fragment in corrected.items():
            if child.startswith(parent+'_sub_'):
                lo,hi=map(int,child[len(parent+'_sub_'):].split('_'))
                if fragment!=sequence[lo:hi]:raise ValueError('Piece sequence does not match exact source interval')
                fragments.append((lo,hi,child,fragment))
        fragments.sort();position=0
        for lo,hi,child,fragment in fragments:
            if lo!=position:raise ValueError('Missing/overlapping source interval')
            position=hi
            accounted.add(child)
            row=next(r for r in names if r['old_name']==child)
            if 'chromosome_assignment_pending' in row.get('flags',''):raise ValueError('Post-cut assignment did not run')
            pieces.append(dict(parent=parent,piece=child,start=lo,end=hi,new_name=row['new_name'],classification=row['class'],
                               sequence_sha256=hashlib.sha256(fragment.encode()).hexdigest()))
        if position!=len(sequence):raise ValueError('Parent reconstruction incomplete')
        reconstruction[parent]=''.join(r[3] for r in fragments)
    if reconstruction!=original or accounted!=set(corrected):raise ValueError('Original sequences do not reconstruct exactly or output contains extra records')
    expected=json.loads(Path(a.actions).read_text())
    if a.scenario=='manual':
        selected=[r for r in expected if r['assembly']==a.assembly]
        expected_parents={r['scaffold'] for r in selected}
        if {p['parent'] for p in pieces}!=expected_parents or len(corrected)!=len(original)+len(selected):
            raise ValueError('Manual route made edits beyond the explicit selected actions')
        for row in selected:
            left=row['scaffold']+'_sub_0_'+str(row['cut_bp'])
            if not corrected[left].endswith('N'*(int(row['gap_end'])-int(row['gap_start']))):raise ValueError('Gap Ns were not retained')
    qc=dict(assembly=a.assembly,scenario=a.scenario,status='PASS',exact_reconstruction=True,
            input_records=len(original),output_records=len(corrected),input_bp=sum(map(len,original.values())),
            output_bp=sum(map(len,corrected.values())),input_n50=n50(original),output_n50=n50(corrected),pieces=pieces)
    Path(a.out).write_text(json.dumps(qc,indent=2)+'\n')


if __name__=='__main__':main()
