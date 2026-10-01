"""Audit exported standard GREF VCF identity and published-tool record counts."""
import argparse
import csv
import json
from pathlib import Path


def sample_mapping(rows, samples):
    expected = {}
    for row in rows:
        name = row['graph_name']
        base, sep, hap = name.rpartition('.')
        column = base if sep and hap.isdigit() else name
        expected.setdefault(column, []).append(row)
    if len(samples) != len(set(samples)) or set(samples) != set(expected):
        raise ValueError('VCF sample columns do not match the graph identity ledger')
    output = []
    for column in samples:
        members = expected[column]
        individuals = {r['sample'] for r in members}
        if len(individuals) != 1:
            raise ValueError('VCF column spans multiple biological individuals')
        output.append((column, next(iter(individuals)), ','.join(r['id'] for r in members), len(members)))
    return output


def audit(directory, identities):
    directory = Path(directory)
    with open(identities) as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    samples = (directory / 'vcf_samples.txt').read_text().splitlines()
    mapping = sample_mapping(rows, samples)
    stats = {}
    for line in (directory / 'bcftools_stats.txt').read_text().splitlines():
        fields = line.split('\t')
        if fields[0] == 'SN':
            if len(fields) != 4 or fields[1] != '0':
                raise ValueError('Unexpected bcftools summary format')
            stats[fields[2].rstrip(':')] = int(fields[3])
    required = ('number of samples', 'number of records')
    if any(k not in stats for k in required) or stats['number of samples'] != len(samples):
        raise ValueError('Missing or inconsistent bcftools summary')
    indexed = 0
    coordinates = set()
    for line in (directory / 'coordinate_records.tsv').read_text().splitlines():
        fields = line.split('\t')
        if len(fields) != 3 or fields[0] in coordinates or int(fields[2]) < 0:
            raise ValueError('Unexpected index statistics format')
        coordinates.add(fields[0])
        indexed += int(fields[2])
    if indexed != stats['number of records']:
        raise ValueError('Index counts disagree with full VCF statistics')
    with (directory / 'vcf_sample_identity.tsv').open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(('vcf_column', 'individual', 'assembly_ids', 'included_haplotypes'))
        writer.writerows(mapping)
    with (directory / 'record_summary.tsv').open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(('metric', 'count'))
        writer.writerows(stats.items())
    result = {'status': 'PASS', 'vcf_columns': len(samples), 'individuals': len({r['sample'] for r in rows}),
              'haplotypes': len(rows), 'records': indexed, 'indexed_coordinate_sequences': len(coordinates),
              'representation': 'Cactus standard GREF(CLIP) export, unchanged'}
    (directory / 'variant_audit.json').write_text(json.dumps(result, indent=2) + '\n')
    report = ['# GREF variant catalog checkpoint', '',
              f"{indexed:,} VCF records; {len(samples)} VCF columns representing {result['individuals']} individuals and {len(rows)} input haplotypes.", '',
              '<details><summary>Catalog statistics and identity</summary>', '',
              '| Metric | Count |', '|---|---:|']
    report.extend(f'| {key} | {value:,} |' for key, value in stats.items())
    report.extend(['', 'Record categories can overlap at multiallelic sites; they are not an exclusive biological event classification. '
                   'These are statistics of the exported representation, not normalized SNP/SV counts. '
                   'The separate reference VCF column remains mapped to its biological individual, not an extra individual. '
                   'This audit does not validate genotype ploidy, allele traversal correspondence or REF sequence against a FASTA.', '',
                   'GREF coordinate sequences include graph-derived sequence. Coordinate counts are not chromosome-only counts; '
                   'do not normalize this VCF against the ordinary assembly FASTA.', '',
                   '[Column identities](vcf_sample_identity.tsv) · [Record summary](record_summary.tsv) · '
                   '[Coordinate index counts](coordinate_records.tsv) · [Full bcftools statistics](bcftools_stats.txt)', '',
                   'Frequency/genotype interpretation and inversion/duplication/translocation calling remain outside this checkpoint.', '', '</details>'])
    (directory / 'variant_report.md').write_text('\n'.join(report) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', required=True)
    parser.add_argument('--identities', required=True)
    args = parser.parse_args()
    audit(args.directory, args.identities)
