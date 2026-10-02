"""Prepare a bounded paired FASTQ prefix for routing tests, not biological inference."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re


def record(handle):
    header = handle.readline()
    if not header:
        return None
    seq, plus, qual = (handle.readline() for _ in range(3))
    if not header.startswith('@') or not plus.startswith('+') or not qual:
        raise ValueError('Malformed/truncated four-line FASTQ')
    if len(seq.rstrip()) != len(qual.rstrip()):
        raise ValueError('Sequence and quality lengths differ')
    key = re.sub(r'/[12]$', '', header.split()[0][1:])
    return key, header + seq + plus + qual


def prepare(r1, r2, out, pairs, sample, taxid):
    if pairs < 1 or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', sample):
        raise ValueError('Invalid pair count or sample ID')
    if r1.resolve() == r2.resolve():
        raise ValueError('Mate files must differ')
    out.mkdir(parents=True, exist_ok=False)
    dest = [out / 'reads_1.fastq.gz', out / 'reads_2.fastq.gz']
    digests = [hashlib.sha256(), hashlib.sha256()]
    count = 0
    opener = lambda p: gzip.open(p, 'rt') if str(p).endswith('.gz') else open(p)
    with opener(r1) as a, opener(r2) as b, gzip.open(dest[0], 'wt') as x, gzip.open(dest[1], 'wt') as y:
        for _ in range(pairs):
            left, right = record(a), record(b)
            if left is None and right is None:
                break
            if left is None or right is None or left[0] != right[0]:
                raise ValueError('Mates are missing or out of order at pair ' + str(count + 1))
            for rec, handle, digest in zip((left, right), (x, y), digests):
                handle.write(rec[1])
                digest.update(rec[1].encode())
            count += 1
    if not count:
        raise ValueError('No read pairs')
    with open(out / 'samples.csv', 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['sample_id', 'taxid', 'ploidy', 'hifi_bam', 'hic_r1', 'hic_r2', 'sr_r1', 'sr_r2'])
        writer.writerow([sample, taxid, 2, '', '', '', str(dest[0].resolve()), str(dest[1].resolve())])
    (out / 'subset.json').write_text(json.dumps({
        'method': 'first N paired records; not random or representative',
        'requested_pairs': pairs, 'written_pairs': count,
        'sources': [str(r1.resolve()), str(r2.resolve())],
        'uncompressed_subset_sha256': [d.hexdigest() for d in digests],
        'scope': 'FASTQ validation covers selected prefix only'}, indent=2) + '\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('r1', type=Path)
    p.add_argument('r2', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--pairs', type=int, default=2000000)
    p.add_argument('--full-input', action='store_true',
                   help='Reference original FASTQs without copying or subsetting')
    p.add_argument('--sample', required=True)
    p.add_argument('--taxid', type=int, required=True)
    a = p.parse_args()
    if a.full_input:
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', a.sample):
            p.error('Invalid sample ID')
        paths = [a.r1.resolve(strict=True), a.r2.resolve(strict=True)]
        if paths[0] == paths[1] or any(not x.is_file() or x.stat().st_size == 0 for x in paths):
            p.error('Need distinct, nonempty mate files')
        a.output.mkdir(parents=True, exist_ok=False)
        with open(a.output / 'samples.csv', 'w', newline='') as handle:
            writer = csv.writer(handle)
            writer.writerow(['sample_id', 'taxid', 'ploidy', 'hifi_bam', 'hic_r1', 'hic_r2', 'sr_r1', 'sr_r2'])
            writer.writerow([a.sample, a.taxid, 2, '', '', '', *map(str, paths)])
        (a.output / 'subset.json').write_text(json.dumps({
            'method': 'full original FASTQs; no subsetting or copying',
            'sources': list(map(str, paths)),
            'scope': 'Existence checks only; FASTQ processing occurs in the pipeline'}, indent=2) + '\n')
    else:
        prepare(a.r1, a.r2, a.output, a.pairs, a.sample, a.taxid)


if __name__ == '__main__':
    main()
