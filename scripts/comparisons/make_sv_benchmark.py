"""Write small deterministic assembly-rearrangement fixtures; no biological inputs."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import random


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def fixtures(seed=20261001, suite='basic'):
    rng = random.Random(seed)
    ref = {name: ''.join(rng.choices('ACGT', k=300000)) for name in ('chr1', 'chr2')}
    if suite == 'repeats':
        # Identical 2 kb repeats on both chromosomes and at event boundaries.
        # Their origins can be unresolvable: engineered coordinates are not
        # necessarily uniquely identifiable breakpoints.
        repeat = ref['chr1'][50000:52000]
        for chrom in ref:
            sequence = ref[chrom]
            for start in (98000, 118000, 138000, 148000, 218000):
                sequence = sequence[:start] + repeat + sequence[start + 2000:]
            ref[chrom] = sequence
    elif suite != 'basic':
        raise ValueError('Unknown fixture suite')
    a, b = ref.values()
    cases = {'control': dict(ref)}
    truth = []

    def add(name, sequence, event, start, end, target, target_start, target_end):
        cases[name] = {'chr1': sequence, 'chr2': b}
        truth.append((name, event, 'chr1', start, end, target, target_start, target_end))

    add('inversion', a[:100000] + reverse_complement(a[100000:140000]) + a[140000:],
        'INV', 100000, 140000, 'chr1', 100000, 140000)
    add('tandem_duplication', a[:140000] + a[120000:140000] + a[140000:],
        'DUP_TANDEM', 120000, 140000, 'chr1', 140000, 160000)
    add('inverted_duplication', a[:220000] + reverse_complement(a[100000:120000]) + a[220000:],
        'DUP_INVERTED_DISPERSED', 100000, 120000, 'chr1', 220000, 240000)
    cases['translocation'] = {'chr1': a[:100000] + a[120000:],
                              'chr2': b[:150000] + a[100000:120000] + b[150000:]}
    truth.append(('translocation', 'CUT_AND_PASTE_INTERCHROM', 'chr1', 100000, 120000,
                  'chr2', 150000, 170000))
    # No biological rearrangement: deliberately destroy chromosome correspondence.
    cases['fragmented_control'] = {f'{name}_part{i + 1}': seq[start:start + 100000]
                                   for name, seq in ref.items()
                                   for i, start in enumerate(range(0, len(seq), 100000))}
    return ref, cases, truth


def write_fasta(path, sequences):
    with path.open('w') as handle:
        for name, sequence in sequences.items():
            handle.write('>' + name + '\n')
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--suite', choices=['basic', 'repeats'], default='basic')
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    ref, cases, truth = fixtures(suite=args.suite)
    write_fasta(args.out / 'reference.fa', ref)
    for name, sequences in cases.items():
        write_fasta(args.out / (name + '.fa'), sequences)
    with (args.out / 'truth.tsv').open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['case', 'event', 'reference_chromosome', 'reference_start', 'reference_end',
                         'query_chromosome', 'query_start', 'query_end'])
        writer.writerows(truth)
    manifest = {'seed': 20261001, 'suite': args.suite, 'coordinates': 'zero-based half-open',
                'cases': list(cases), 'negative_controls': ['control', 'fragmented_control'],
                'scope': 'Bounded synthetic benchmark; not genomewide accuracy or graph-allele validation',
                'repeat_design': ('Identical 2kb blocks at 98000,118000,138000,148000,218000 on both chromosomes; breakpoint origin may be ambiguous'
                                  if args.suite == 'repeats' else 'No deliberately planted repeats'),
                'sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                           for p in sorted(args.out.glob('*.fa'))}}
    (args.out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')


if __name__ == '__main__':
    main()
