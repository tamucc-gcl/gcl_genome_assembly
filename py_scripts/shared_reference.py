#!/usr/bin/env python3
"""Record the exact final linear reference, independent of graph execution."""
import argparse
import json
from pathlib import Path
from pangenome_input_audit import inspect_fasta


def verify(reference, fasta, fai):
    digest, lengths = inspect_fasta(fasta)
    expected = {}
    for line in Path(fai).read_text(encoding='utf-8').splitlines():
        fields = line.split('\t')
        if len(fields) < 2 or fields[0] in expected:
            raise ValueError('Malformed or duplicate reference FAI record')
        expected[fields[0]] = int(fields[1])
    if lengths != expected:
        raise ValueError('Reference FASTA/FAI sequence identities or lengths disagree')
    if set(reference['reference_contigs']) - set(lengths):
        raise ValueError('Selected reference chromosomes absent from final FASTA')
    return dict(reference, sha256=digest, sequence_count=len(lengths), bases=sum(lengths.values()),
                sequences=[dict(name=n, length=s) for n, s in lengths.items()],
                coordinate_stage='final',
                mapping_compatibility='Require identical FASTA SHA256 and sequence dictionary; assembly ID alone is insufficient')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--selection', required=True)
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--fai', required=True)
    args = parser.parse_args()
    result = verify(json.loads(Path(args.selection).read_text(encoding='utf-8')), args.fasta, args.fai)
    Path('reference.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    Path('reference.md').write_text(
        '# Final linear reference\n\n'
        f"Species: {result['species']} (taxid {result['taxid']}). Reference: **{result['id']}**.\n\n"
        f"Selection: `{result['selection']}`. Ranking: `{result['ranking_basis']}`.\n\n"
        f"Final FASTA SHA256: `{result['sha256']}`.\n\n"
        'This real assembly supplies the shared linear coordinate frame. Its selection does not '
        'certify biological correctness or imply that a graph was built. Graph admission is assessed separately.\n\n'
        '[Reference provenance and sequence dictionary](reference.json). Existing mappings must match '
        'this FASTA checksum and sequence dictionary before reuse.\n', encoding='utf-8')


if __name__ == '__main__':
    main()
