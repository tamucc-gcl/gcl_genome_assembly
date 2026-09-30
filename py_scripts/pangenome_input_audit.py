#!/usr/bin/env python3
"""Validate a same-species graph cohort before expensive construction."""
import argparse
import csv
import hashlib
import json
from pathlib import Path


def inspect_fasta(path):
    digest = hashlib.sha256()
    lengths = {}
    name = None
    with open(path, 'rb') as handle:
        for line in handle:
            digest.update(line)
            if line.startswith(b'>'):
                tokens = line[1:].split()
                if not tokens:
                    raise ValueError('Empty FASTA sequence identifier')
                name = tokens[0].decode('ascii')
                if name in lengths:
                    raise ValueError('Duplicate FASTA sequence: ' + name)
                lengths[name] = 0
            elif line.strip():
                if name is None:
                    raise ValueError('Sequence before FASTA header')
                lengths[name] += len(b''.join(line.split()))
    if not lengths or any(n == 0 for n in lengths.values()):
        raise ValueError('Empty FASTA or zero-length sequence')
    return digest.hexdigest(), lengths


def audit(cohort):
    members = cohort['members']
    if len({m['id'] for m in members}) != len(members):
        raise ValueError('Duplicate assembly IDs')
    if len({m['graph_name'] for m in members}) != len(members):
        raise ValueError('Duplicate graph names')
    if len({m['sample'] for m in members}) < cohort['min_individuals']:
        raise ValueError('Insufficient biological individuals')
    if not cohort['reference_contigs'] or len(set(cohort['reference_contigs'])) != len(cohort['reference_contigs']):
        raise ValueError('Empty or duplicate reference chromosome IDs')
    refs = [m for m in members if m['id'] == cohort['reference_id']]
    if len(refs) != 1 or refs[0]['graph_name'] != cohort['reference_name']:
        raise ValueError('Reference assembly and graph identity disagree')
    audited = []
    for m in members:
        if str(m['taxid']) != str(cohort['taxid']):
            raise ValueError('Cross-species graph input')
        digest, lengths = inspect_fasta(m['staged_fasta'])
        if m['id'] == cohort['reference_id']:
            missing = set(cohort['reference_contigs']) - set(lengths)
            if missing:
                raise ValueError('Reference chromosomes absent from FASTA: ' + ','.join(sorted(missing)))
        audited.append(dict(m, sha256=digest, sequence_count=len(lengths), bases=sum(lengths.values()),
                            is_reference=m['id'] == cohort['reference_id']))
    return audited


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input', required=True)
    a = p.parse_args()
    cohort = json.loads(Path(a.input).read_text())
    rows = audit(cohort)
    columns = ['taxid', 'id', 'sample', 'haplotype', 'graph_name', 'is_reference', 'fasta', 'sha256', 'sequence_count', 'bases']
    with open('pangenome_identity.tsv', 'w') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction='ignore', delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    Path('pangenome_input_audit.json').write_text(json.dumps(dict(cohort=cohort, audited=rows), indent=2)+'\n')
    Path('input.ok').write_text('same-species identities, FASTA records and reference chromosomes verified\n')
    Path('pangenome_input_report.md').write_text(
        '<details>\n<summary>Pangenome input audit: taxid '+str(cohort['taxid'])+'</summary>\n\n'
        'Input checks passed for '+str(len(rows))+' assemblies from '+str(len({r['sample'] for r in rows}))+' individuals.\n\n'
        'Reference: '+cohort['reference_id']+'. Selection: '+cohort.get('reference_selection', 'not recorded')+'.\n\n'
        '[Biological identities and FASTA checksums](pangenome_identity.tsv).\n\n'
        'This report validates inputs only; it does not confirm a graph was built or its biological quality.\n\n</details>\n')


if __name__ == '__main__':
    main()
