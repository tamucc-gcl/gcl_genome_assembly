"""Collect deterministic small GREF windows and genotype/annotation diagnostics."""
import argparse
from collections import Counter
import csv
import gzip
import json
from pathlib import Path
import re
import shutil
import subprocess


def choose_windows(rows, count=12, width=10000):
    """Sample coordinate space, not variant mechanisms or population frequencies."""
    selected = []
    categories = {'reference_named': [], 'alt_named': [], 'other': []}
    for chrom, length, records in rows:
        if length <= 0 or records <= 0:
            continue
        category = 'reference_named' if re.fullmatch(r'chr\d+_\d+', chrom) else 'alt_named' if chrom.endswith('_alt') else 'other'
        categories[category].append((chrom, length, records))
    for category, candidates in categories.items():
        # Span length ranks to include short and long graph-derived coordinates.
        candidates.sort(key=lambda row: (row[1], row[0]))
        ranks = sorted({round(i * (len(candidates) - 1) / max(1, min(count, len(candidates)) - 1))
                        for i in range(min(count, len(candidates)))})
        for rank in ranks:
            chrom, length, records = candidates[rank]
            for fraction in (0, 0.5, 1):
                start = int(max(0, length - width) * fraction)
                selected.append((chrom, start, min(length, start + width), category))
    merged = []
    for chrom, start, end, category in sorted(set(selected)):
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            previous = merged[-1]
            merged[-1] = (chrom, previous[1], max(previous[2], end), category)
        else:
            merged.append((chrom, start, end, category))
    return merged


def ledger_slots(rows):
    """Expected positions under the numeric graph-haplotype export convention.

    This is a hypothesis checked against GT, not inferred biological ploidy.
    """
    slots = {}
    for row in rows:
        name = row['graph_name']
        match = re.fullmatch(r'(.+)\.(\d+)', name)
        column, index = (match[1], int(match[2]) - 1) if match else (name, 0)
        if index < 0 or index in slots.setdefault(column, set()):
            raise ValueError('Duplicate or invalid graph haplotype: ' + name)
        slots[column].add(index)
    return slots


def slot_accounting(gt, expected):
    alleles = re.split(r'[/|]', gt)
    if expected is None:
        return Counter(unresolved_records=1, unknown_column=1)
    if max(expected) >= len(alleles) or ('/' in gt and expected != set(range(len(alleles)))):
        return Counter(unresolved_records=1, unresolved_biological_slots=len(expected))
    counts = Counter(resolved_records=1)
    for index, allele in enumerate(alleles):
        kind = 'biological' if index in expected else 'placeholder'
        counts[kind + ('_missing' if allele == '.' else '_called')] += 1
    return counts


def inspect_vcf(path, out, label, expected_slots=None):
    samples = []
    shapes = Counter()
    flags = Counter()
    examples = []
    slot_counts = {}
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            if line.startswith('#CHROM'):
                samples = line.rstrip().split('\t')[9:]
                if expected_slots is not None and set(samples) != set(expected_slots):
                    raise ValueError('VCF columns and identity ledger disagree')
            if line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            chrom, pos, identifier, ref, alt = fields[:5]
            alts = alt.split(',')
            info = dict(x.split('=', 1) if '=' in x else (x, '') for x in fields[7].split(';'))
            formats = fields[8].split(':') if len(fields) > 8 else []
            issues = []
            ac = [0] * len(alts)
            an = 0
            if 'GT' not in formats:
                issues.append('no_GT')
            else:
                gt_index = formats.index('GT')
                for sample, value in zip(samples, fields[9:]):
                    values = value.split(':')
                    gt = values[gt_index] if gt_index < len(values) else '.'
                    if expected_slots is not None:
                        accounting = slot_accounting(gt, expected_slots.get(sample))
                        slot_counts.setdefault(sample, Counter()).update(accounting)
                        if accounting['placeholder_called']:
                            issues.append('unexpected_called_placeholder')
                        if accounting['unresolved_records']:
                            issues.append('unresolved_slot_assignment')
                    alleles = re.split(r'[/|]', gt)
                    shape = (len(alleles), 'phased' if '|' in gt else 'unphased_or_haploid',
                             sum(x == '.' for x in alleles))
                    shapes[(sample, *shape)] += 1
                    for allele in alleles:
                        if allele != '.':
                            value = int(allele)
                            if not 0 <= value <= len(alts):
                                issues.append('GT_out_of_range')
                                continue
                            an += 1
                            if value:
                                ac[value - 1] += 1
                if 'AN' in info and info['AN'] != '.' and int(info['AN']) != an:
                    issues.append('AN_differs_from_exported_GT')
                if 'AC' in info and info['AC'] != '.' and info['AC'].split(',') != list(map(str, ac)):
                    issues.append('AC_differs_from_exported_GT')
            if 'AT' in info and len(info['AT'].split(',')) != len(alts) + 1:
                issues.append('AT_cardinality_mismatch')
            flags['records'] += 1
            flags.update(set(issues))
            if len(examples) < 100 or (issues and len(examples) < 300):
                examples.append((chrom, pos, identifier, len(ref), ','.join(map(lambda x: str(len(x)), alts)),
                                 info.get('LV', ''), info.get('PS', ''), info.get('RC', ''), ';'.join(issues)))
    for suffix, header, rows in [
        ('genotype_shapes.tsv', ('vcf_column', 'GT_slots', 'separator', 'missing_slots', 'records'),
         [(*key, value) for key, value in sorted(shapes.items())]),
        ('examples.tsv', ('coordinate', 'position', 'id', 'REF_length', 'ALT_lengths', 'LV', 'PS', 'RC', 'flags'), examples)]:
        with (out / (label + '.' + suffix)).open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(header)
            writer.writerows(rows)
    if expected_slots is not None:
        columns = ['resolved_records', 'unresolved_records', 'biological_called',
                   'biological_missing', 'placeholder_missing', 'placeholder_called',
                   'unresolved_biological_slots']
        with (out / (label + '.biological_slots.tsv')).open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(['vcf_column', 'expected_slots_1based', *columns])
            for sample in samples:
                counts = slot_counts.get(sample, Counter())
                writer.writerow([sample, ','.join(str(x + 1) for x in sorted(expected_slots[sample])),
                                 *(counts[x] for x in columns)])
    return dict(flags)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pangenome-dir', type=Path, required=True)
    parser.add_argument('--taxid', required=True)
    parser.add_argument('--bcftools', default='bcftools')
    parser.add_argument('--conda-cache', type=Path,
                        help='If bcftools is absent from PATH, find an existing 1.21 executable here; no install')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if not shutil.which(args.bcftools) and args.conda_cache:
        for candidate in sorted(args.conda_cache.glob('env-*/bin/bcftools')):
            version = subprocess.run([str(candidate), '--version-only'], capture_output=True, text=True)
            if version.returncode == 0 and re.match(r'^1\.21(?:\+|\s|$)', version.stdout):
                args.bcftools = str(candidate.resolve())
                break
    if not shutil.which(args.bcftools):
        raise ValueError('Supply --bcftools with an existing executable or --conda-cache; nothing was installed')
    args.out.mkdir(parents=True, exist_ok=False)
    source = args.pangenome_dir.resolve()
    with (source / 'pangenome_identity.tsv').open() as handle:
        expected_slots = ledger_slots(csv.DictReader(handle, delimiter='\t'))
    vcf = source / (args.taxid + '.gref.vcf.gz')
    commands = []

    def run(arguments, destination=None):
        command = [args.bcftools, *map(str, arguments)]
        commands.append(command)
        if destination:
            with destination.open('w') as handle:
                subprocess.run(command, stdout=handle, check=True)
        else:
            subprocess.run(command, check=True)

    run(['--version'], args.out / 'bcftools_version.txt')
    run(['index', '--stats', vcf], args.out / 'coordinate_records.tsv')
    with (args.out / 'coordinate_records.tsv').open() as handle:
        records = [(f[0], int(f[1]), int(f[2])) for f in csv.reader(handle, delimiter='\t') if f[1] != '.']
    windows = choose_windows(records)
    if not windows:
        raise ValueError('No populated coordinate windows found')
    bed = args.out / 'windows.bed'
    with bed.open('w') as handle:
        for chrom, start, end, category in windows:
            handle.write(f'{chrom}\t{start}\t{end}\t{category}\n')
    results = {}
    for label, candidate in [('standard', vcf), ('raw', source / (args.taxid + '.gref.raw.vcf.gz'))]:
        if not candidate.exists() or not any(Path(str(candidate) + x).exists() for x in ('.tbi', '.csi')):
            if label == 'standard':
                raise ValueError('Standard GREF VCF/index missing')
            results[label] = {'status': 'SKIPPED', 'reason': 'No existing indexed raw GREF VCF; not creating a new index'}
            continue
        run(['view', '-h', candidate], args.out / (label + '.header.txt'))
        output = args.out / (label + '.windows.vcf.gz')
        run(['view', '--regions-overlap', '0', '-R', bed, '-Oz', '-o', output, candidate])
        results[label] = inspect_vcf(output, args.out, label, expected_slots)
    for name in ('pangenome_identity.tsv', 'pangenome_manifest.tsv'):
        shutil.copy2(source / name, args.out / name)
    # Small prefix establishes the actual segment-map schema; full map remains on cluster.
    segments = source / (args.taxid + '.gref.gref-segs.tsv.gz')
    if segments.exists():
        with gzip.open(segments, 'rt') as handle, (args.out / 'gref_segments_prefix.tsv').open('w') as dest:
            for _, line in zip(range(101), handle):
                dest.write(line)
    (args.out / 'inspection.json').write_text(json.dumps({'source': str(source), 'windows': len(windows),
        'selection': 'Length-rank coordinate strata with start/middle/end 10kb windows; not random or representative',
        'scope': 'Exported GT slots and annotation consistency only, not independent biological ploidy or SV accuracy',
        'results': results, 'commands': commands}, indent=2) + '\n')


if __name__ == '__main__':
    main()
