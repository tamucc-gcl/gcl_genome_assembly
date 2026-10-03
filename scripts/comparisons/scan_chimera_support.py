#!/usr/bin/env python3
"""Read-only diagnostic of retained context SAMs/BAMs; never authorizes cuts."""
import argparse
import bisect
import csv
import gzip
import json
import re
import shutil
import subprocess
from pathlib import Path


def alignment(fields):
    pos = int(fields[3])-1
    start = pos
    blocks, chain, columns = [], None, 0
    for size, op in re.findall(r'(\d+)([MIDNSHP=X])', fields[5]):
        size = int(size)
        if op in 'M=X':
            if chain is None:
                chain = pos
            pos += size
        elif op in 'DN' or (op == 'I' and size > 50):
            if op == 'N' or size > 50:
                if chain is not None:
                    blocks.append((chain, pos))
                chain = None
            if op in 'DN':
                pos += size
        if op in 'M=XID':
            columns += size
    if chain is not None:
        blocks.append((chain, pos))
    tags = dict(f.split(':', 2)[::2] for f in fields[11:] if f.count(':') >= 2)
    rate = int(tags['NM'])/columns if 'NM' in tags and columns else None
    return start, pos, blocks, rate


def scan(sam, interval, sequence, step=1000, anchor=1000):
    positions = list(range(interval['start'], interval['end'], step))
    tracks = {key: [0]*(len(positions)+1) for key in
              ('primary_span_any_mapq', 'primary_span_mapq20', 'chain_mapq20_nm1pct',
               'chain_mapq20_nm2pct')}
    names = set()
    def add(key, lo, hi):
        left = bisect.bisect_left(positions, lo+anchor)
        right = bisect.bisect_right(positions, hi-anchor)
        if left < right:
            tracks[key][left] += 1
            tracks[key][right] -= 1
    with sam.open() as handle:
        for line in handle:
            if line.startswith('@'):
                continue
            f = line.rstrip().split('\t')
            if len(f) < 11:
                raise ValueError('Malformed SAM record')
            names.add(f[0])
            if int(f[1]) & 0xF04 or f[2] != interval['scaffold']:
                continue
            lo, hi, blocks, rate = alignment(f)
            add('primary_span_any_mapq', lo, hi)
            if int(f[4]) < 20 or int(f[4]) == 255:
                continue
            add('primary_span_mapq20', lo, hi)
            for limit, key in ((0.01, 'chain_mapq20_nm1pct'), (0.02, 'chain_mapq20_nm2pct')):
                if rate is not None and rate <= limit:
                    for lo, hi in blocks:
                        add(key, lo, hi)
    totals = dict.fromkeys(tracks, 0)
    output = []
    for i, pos in enumerate(positions):
        for key in tracks:
            totals[key] += tracks[key][i]
        offset = pos-interval['start']
        bases = sequence[max(0, offset-anchor):offset+anchor]
        output.append(dict(position_0based=pos,
                           inside_transition=interval['lo'] <= pos < interval['hi'],
                           complete_sequence_window=offset >= anchor and offset+anchor <= len(sequence),
                           n_bases_in_2kb=bases.upper().count('N'), **totals))
    return names, output


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('context', type=Path)
    p.add_argument('out', type=Path)
    p.add_argument('--samtools', default='samtools')
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=False)
    manifests = sorted(a.context.glob('*.sequence_context/provenance.json'))
    if not manifests:
        p.error('No retained sequence-context provenance found')
    completed = []
    for manifest in manifests:
        meta = json.loads(manifest.read_text())
        if not meta.get('hifi'):
            continue
        if meta.get('coordinate_stage') != 'pre_finishing':
            raise ValueError('Expected pre-finishing context')
        root = manifest.parent
        dest = a.out/meta['assembly']
        dest.mkdir()
        shutil.copy2(manifest, dest/'provenance.json')
        names = set()
        for key, iv in meta['intervals'].items():
            fasta = root/(key+'.fa')
            sequence = ''.join(line.strip() for line in fasta.read_text().splitlines()
                               if not line.startswith('>'))
            if len(sequence) != iv['end']-iv['start']:
                raise ValueError('Candidate sequence length disagrees with provenance')
            shutil.copy2(fasta, dest/fasta.name)
            ids, report = scan(root/(key+'.sam'), iv, sequence)
            names.update(ids)
            with (dest/(key+'.support.tsv')).open('w') as handle:
                writer = csv.DictWriter(handle, fieldnames=list(report[0]), delimiter='\t')
                writer.writeheader()
                writer.writerows(report)
        name_file = dest/'regional_read_names.txt'
        name_file.write_text(''.join(n+'\n' for n in sorted(names)))
        bam = root/'hifi.bam'
        subprocess.run([a.samtools, 'quickcheck', str(bam)], check=True)
        # Scan each retained BAM once, retaining all reported placements for the
        # selected read names. SEQ/QUAL omitted; CIGAR and all tags retained.
        with (dest/'samtools.stderr.txt').open('w') as err, gzip.open(dest/'all_placements.tsv.gz', 'wt') as out:
            out.write('qname\tflag\trname\tpos_1based\tmapq\tcigar\trnext\tpnext\ttlen\ttags\n')
            proc = subprocess.Popen([a.samtools, 'view', '-N', str(name_file), str(bam)],
                                    stdout=subprocess.PIPE, stderr=err, text=True)
            for line in proc.stdout:
                f = line.rstrip().split('\t')
                out.write('\t'.join(f[:9]+[';'.join(f[11:])])+'\n')
            if proc.wait():
                raise RuntimeError('samtools placement export failed; see stderr')
        completed.append(dict(assembly=meta['assembly'], candidates=len(meta['intervals']),
                              regional_read_names=len(names), assessment_sha256=meta['sha256']))
    if not completed:
        raise ValueError('No HiFi-enabled contexts found')
    (a.out/'status.json').write_text(json.dumps(dict(status='SUCCESS', assemblies=completed,
        step_bp=1000, anchor_bp=1000, max_internal_indel_bp=50,
        interpretation='Descriptive sensitivity tracks, not cut calls. NM rate is alignment-wide. '
        'Small indels are allowed; mapping quality does not establish unique anchors. '
        'Overlapping windows are not independent observations. Assembly SHA is recorded provenance, '
        'not an independent rehash of the unavailable assessment FASTA.'), indent=2)+'\n')


if __name__ == '__main__':
    main()
