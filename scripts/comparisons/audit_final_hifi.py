"""Read-only final BAM provenance audit; does not map reads or scan BAM alignments."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('results', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--samtools', required=True)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    data = json.loads((a.results / 'assembly_eligibility.json').read_text())
    root = a.results / 'bam/hifi/final'
    rows = []
    for record in data['assemblies']:
        aid = record['id']
        bam = root / (aid + '.sorted.bam')
        row = dict(assembly=aid, status='NOT_AVAILABLE', reason='', expected_sha256='', actual_sha256='')
        if not bam.exists():
            row['reason'] = 'no_published_final_BAM'
            rows.append(row)
            continue
        try:
            fai = root / (aid + '.reference.fai')
            checksum = root / (aid + '.reference.sha256')
            final = a.results / 'assembly/final' / (aid + '.fasta')
            for path in (bam, Path(str(bam) + '.bai'), fai, checksum, final):
                if not path.is_file() or path.stat().st_size == 0:
                    raise ValueError('Missing/empty file: ' + str(path))
            subprocess.run([a.samtools, 'quickcheck', '-v', str(bam)], check=True, capture_output=True, text=True)
            row['expected_sha256'] = checksum.read_text().split()[0]
            row['actual_sha256'] = sha256(final)
            if row['expected_sha256'] != row['actual_sha256']:
                raise ValueError('BAM mapping reference checksum differs from published final FASTA')
            expected = [tuple(line.split('\t')[:2]) for line in fai.read_text().splitlines()]
            header = subprocess.check_output([a.samtools, 'view', '-H', str(bam)], text=True)
            observed = []
            for line in header.splitlines():
                if line.startswith('@SQ\t'):
                    tags = dict(field.split(':', 1) for field in line.split('\t')[1:])
                    observed.append((tags['SN'], tags['LN']))
            if observed != expected:
                raise ValueError('BAM dictionary differs from retained reference FAI')
            (a.output / (aid + '.header.sam')).write_text(header)
            stats = subprocess.check_output([a.samtools, 'idxstats', str(bam)], text=True)
            (a.output / (aid + '.idxstats.tsv')).write_text(stats)
            row.update(status='PASS_PROVENANCE_CHECKS', reason='checksum_dictionary_index_quickcheck')
        except (ValueError, OSError, subprocess.CalledProcessError) as exc:
            row.update(status='REVIEW', reason=str(exc))
        rows.append(row)
    with open(a.output / 'final_hifi_audit.tsv', 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['assembly', 'status', 'reason', 'expected_sha256', 'actual_sha256'], delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    (a.output / 'README.md').write_text('# Final HiFi BAM provenance audit\n\n'
        'PASS verifies the recorded mapping reference checksum, BAM dictionary, index access and quickcheck. '
        'It does not scan every alignment, verify read identities against raw inputs, or establish biological SV support. '
        'NOT_AVAILABLE is explicit and must be reconciled against the expected HiFi cohort. '
        'Final FASTAs are read for hashing; BAM alignments are not remapped.\n')
    if any(r['status'] == 'REVIEW' for r in rows):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
