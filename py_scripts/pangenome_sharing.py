"""Identity adapter and summary of Panacus outputs; never counts graph nodes."""
import argparse
import csv
import gzip
import json
import math
from pathlib import Path


def prepare(gfa, ledger, out):
    with open(ledger) as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    aliases = {}
    for row in rows:
        name = row['graph_name']
        base, sep, hap = name.rpartition('.')
        key = (base, hap) if sep and hap.isdigit() else (name, '0')
        if key in aliases or not row['sample'] or not row['id']:
            raise ValueError('Invalid or duplicate identity: ' + name)
        aliases[key] = row
    if not rows or len({r['id'] for r in rows}) != len(rows):
        raise ValueError('Empty ledger or duplicate assembly IDs')
    paths = {}
    seen = set()
    opener = gzip.open if str(gfa).endswith('.gz') else open
    with opener(gfa, 'rt') as handle:
        for line in handle:
            if line.startswith('W\t'):
                fields = line.split('\t', 6)
                key = (fields[1], fields[2])
                path = '#'.join(fields[1:4])
            elif line.startswith('P\t'):
                path = line.split('\t', 2)[1]
                parts = path.split('#')
                if len(parts) != 3:
                    raise ValueError('Non-PanSN path: ' + path)
                key = tuple(parts[:2])
            else:
                continue
            if key not in aliases:
                raise ValueError('Graph path absent from identity ledger: ' + path)
            row = aliases[key]
            paths[path] = row
            seen.add(row['id'])
    missing = {r['id'] for r in rows} - seen
    if missing:
        raise ValueError('Ledger assemblies missing from graph: ' + ', '.join(sorted(missing)))
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    for unit, column in [('haplotype', 'id'), ('individual', 'sample')]:
        with (out / (unit + '.groups.tsv')).open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerows((p, r[column]) for p, r in sorted(paths.items()))
    metadata = {'haplotype': len(rows), 'individual': len({r['sample'] for r in rows}),
                'path_names': len(paths)}
    (out / 'denominators.json').write_text(json.dumps(metadata, indent=2) + '\n')
    return metadata


def read_hist(path, n):
    bins = {}
    headers = []
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith('#'):
            continue
        fields = line.split('\t')
        if fields[0].isdigit():
            if len(fields) != 2:
                raise ValueError('Expected one bp histogram column')
            k, value = map(int, fields)
            if k in bins or not 0 <= k <= n or value < 0:
                raise ValueError('Invalid histogram bin: ' + line)
            bins[k] = value
        else:
            headers.append(fields)
    if ['count', 'bp'] not in headers or set(bins) != set(range(n + 1)):
        raise ValueError('Panacus histogram format/denominator mismatch: ' + str(path))
    return bins


def cutoff(value, n):
    if not math.isfinite(value) or value == 0:
        raise ValueError('Tier thresholds must be finite and nonzero')
    return max(1, min(n, math.ceil(n * value if 0 < value <= 1 else n + value if value < 0 else value)))


def summarize(directory, core, softcore, shell):
    directory = Path(directory)
    meta = json.loads((directory / 'denominators.json').read_text())
    records = []
    totals = []
    for unit in ('haplotype', 'individual'):
        n = meta[unit]
        bins = read_hist(directory / (unit + '.hist.tsv'), n)
        cuts = [cutoff(v, n) for v in (core, softcore, shell)]
        if cuts[0] < cuts[1]:
            raise ValueError('Core cutoff must be at least the softcore cutoff')
        total = sum(v for k, v in bins.items() if k > 0)
        totals.append(total + bins[0])
        classes = dict.fromkeys(('core', 'softcore', 'shell', 'cloud'), 0)
        for k, bp in bins.items():
            if k:
                category = next((label for label, cut in zip(('core', 'softcore', 'shell'), cuts) if k >= cut), 'cloud')
                classes[category] += bp
        # Private is a separate, potentially overlapping measure, not a fifth tier.
        classes.update(private=bins.get(1, 0), unrepresented=bins[0])
        for category, bp in classes.items():
            records.append((unit, n, category, bp, '' if category == 'unrepresented' or not total else f'{100 * bp / total:.6f}', *cuts))
    if totals[0] != totals[1]:
        raise ValueError('Haplotype and individual histograms disagree on total graph bp')
    with (directory / 'sharing_summary.tsv').open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(('unit', 'denominator', 'class', 'bp', 'percent_represented_graph_bp', 'core_min', 'softcore_min', 'shell_min'))
        writer.writerows(records)
    lines = ['# CLIP sequence sharing', '',
             f"Cohort: {meta['haplotype']} haplotypes from {meta['individual']} individuals.", '',
             'Counts describe distinct graph sequence, not assembly lengths, copy number or independent population samples. '
             'Individual presence means present in at least one included haplotype. Missing/clipped sequence can affect apparent absence.', '',
             '<details><summary>Sharing categories and denominators</summary>', '',
             '| Unit | N | Class | bp | % represented graph bp |', '|---|---:|---|---:|---:|']
    lines.extend(f'| {u} | {n} | {c} | {bp} | {pct} |' for u, n, c, bp, pct, *_ in records)
    lines.extend(['', 'Core, softcore, shell and cloud are disjoint, assigned in that order. '
                  'Private (coverage = 1) is a separate overlapping statistic. Unrepresented (coverage = 0) is excluded from percentages. '
                  'Realized minimum counts are recorded in sharing_summary.tsv; small cohorts can have empty tiers.', '',
                  'Growth curves are descriptive sequence accumulation, not evidence for an open/closed genome or a population-genetic inference.', '', '</details>', ''])
    (directory / 'sharing_report.md').write_text('\n'.join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    prep = commands.add_parser('prepare')
    prep.add_argument('--gfa', required=True)
    prep.add_argument('--ledger', required=True)
    prep.add_argument('--out', required=True)
    summary = commands.add_parser('summarize')
    summary.add_argument('--directory', required=True)
    for name, default in [('core', 1), ('softcore', -1), ('shell', 2)]:
        summary.add_argument('--' + name, type=float, default=default)
    args = parser.parse_args()
    if args.command == 'prepare':
        prepare(args.gfa, args.ledger, args.out)
    else:
        summarize(args.directory, args.core, args.softcore, args.shell)


if __name__ == '__main__':
    main()
