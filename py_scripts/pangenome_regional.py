"""Aggregate Panacus global node coverage onto biological paths; no coverage inference."""
import argparse
from array import array
from collections import defaultdict
import csv
import gzip
import json
import math
from pathlib import Path
import re
import tempfile

import numpy as np


def text_open(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def chromosome_label(sequence):
    match = re.fullmatch(r'chr(\d+)_\d+', sequence)
    if match:
        return 'chr' + match.group(1)
    if '+' in sequence and all(re.fullmatch(r'chr\d+_\d+', part) for part in sequence.split('+')):
        return 'composite:' + sequence
    return 'unplaced'


def aggregate(gfa, totals, groups_file, histogram, out, core=1.0, softcore=-1.0, shell=2.0):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    with open(groups_file) as handle:
        pairs = list(csv.reader(handle, delimiter='\t'))
    groups = dict(pairs)
    if len(groups) != len(pairs) or any(len(row) != 2 for row in pairs):
        raise ValueError('Invalid path group map')
    haplotypes = sorted(set(groups.values()))
    n = len(haplotypes)
    if not n:
        raise ValueError('Empty haplotype cohort')
    def cutoff(value):
        if not math.isfinite(value) or value == 0:
            raise ValueError('Tier thresholds must be finite and nonzero')
        return max(1, min(n, math.ceil(n * value if 0 < value <= 1 else n + value if value < 0 else value)))
    cuts = list(map(cutoff, (core, softcore, shell)))
    if cuts[0] < cuts[1]:
        raise ValueError('Core cutoff below softcore')
    ids, lengths = array('Q'), array('Q')
    offsets = defaultdict(list)
    with tempfile.TemporaryDirectory(prefix='regional-', dir=out.parent) as tmp:
        # Retain only path walks on local scratch; no node-by-haplotype matrix.
        with open(Path(tmp) / 'walks.txt', 'w+') as walks:
            with text_open(gfa) as handle:
                for line in handle:
                    if line.startswith('S\t'):
                        fields = line.rstrip('\n').split('\t')
                        ids.append(int(fields[1]))
                        if fields[2] == '*':
                            tag = next((x[5:] for x in fields[3:] if x.startswith('LN:i:')), None)
                            if tag is None:
                                raise ValueError('Segment lacks sequence and LN tag')
                            lengths.append(int(tag))
                        else:
                            lengths.append(len(fields[2]))
                    elif line.startswith(('W\t', 'P\t')):
                        fields = line.rstrip('\n').split('\t')
                        if fields[0] == 'W':
                            path = '#'.join(fields[1:4])
                            walk = fields[6]
                            kind = 'W'
                        else:
                            path, walk, kind = fields[1], fields[2], 'P'
                        if path not in groups:
                            raise ValueError('Unmapped path: ' + path)
                        offsets[path].append(walks.tell())
                        walks.write(kind + '\t' + walk + '\n')
            if set(offsets) != set(groups):
                raise ValueError('Graph paths and audited group map disagree')
            node_ids = np.asarray(ids, dtype=np.uint64)
            node_lengths = np.asarray(lengths, dtype=np.uint64)
            order = np.argsort(node_ids)
            node_ids, node_lengths = node_ids[order], node_lengths[order]
            del ids, lengths, order
            if not len(node_ids) or np.any(node_ids[1:] == node_ids[:-1]):
                raise ValueError('Empty graph or duplicate segment IDs')

            def indices(values):
                values = np.asarray(values, dtype=np.uint64)
                index = np.searchsorted(node_ids, values)
                if np.any(index >= len(node_ids)) or np.any(node_ids[index] != values):
                    raise ValueError('Unknown segment ID')
                return index

            coverage = np.full(len(node_ids), np.iinfo(np.uint32).max, dtype=np.uint32)
            batch_ids, batch_cov = [], []

            def put_batch():
                index = indices(batch_ids)
                if len(np.unique(index)) != len(index) or np.any(coverage[index] != np.iinfo(np.uint32).max):
                    raise ValueError('Duplicate node in Panacus total table')
                coverage[index] = batch_cov

            header_seen = False
            with text_open(totals) as handle:
                for line in handle:
                    if not line.strip() or line.startswith('#'):
                        continue
                    fields = line.rstrip().split('\t')
                    if fields == ['panacus', 'table']:
                        continue
                    if fields == ['bp', 'total']:
                        header_seen = True
                        continue
                    if not header_seen or len(fields) != 2:
                        raise ValueError('Unexpected Panacus total table format')
                    node, count = map(int, fields)
                    if not 0 <= count <= n:
                        raise ValueError('Global coverage outside cohort bounds')
                    batch_ids.append(node)
                    batch_cov.append(count)
                    if len(batch_ids) == 100000:
                        put_batch()
                        batch_ids, batch_cov = [], []
            if batch_ids:
                put_batch()
            if not header_seen or np.any(coverage > n):
                raise ValueError('Missing Panacus coverage for graph segments')
            expected = {}
            with open(histogram) as handle:
                for line in handle:
                    fields = line.rstrip().split('\t')
                    if fields[0].isdigit():
                        expected[int(fields[0])] = int(fields[1])

            def bins(index):
                return np.bincount(coverage[index].astype(np.int64),
                                   weights=node_lengths[index], minlength=n + 1).astype(np.uint64)

            global_hist = np.zeros(n + 1, dtype=np.uint64)
            for start in range(0, len(node_ids), 1000000):
                global_hist += bins(slice(start, start + 1000000))
            if expected != dict(enumerate(map(int, global_hist))):
                raise ValueError('Global node coverage does not reconcile with accepted Panacus histogram')
            seen = np.zeros(len(node_ids), dtype=np.uint32)
            chromosome_seen = np.zeros(len(node_ids), dtype=np.uint32)
            stamp = 0

            def walk_indices(paths):
                for path in paths:
                    for offset in offsets[path]:
                        walks.seek(offset)
                        kind, walk = walks.readline().rstrip().split('\t', 1)
                        if kind == 'W':
                            values = np.fromiter((int(x.group(1)) for x in re.finditer(r'[<>](\d+)', walk)), dtype=np.uint64)
                        else:
                            values = np.fromiter((int(x[:-1]) for x in walk.split(',')), dtype=np.uint64)
                        if not len(values):
                            raise ValueError('Empty or unsupported graph walk')
                        yield np.unique(indices(values))

            by_chromosome = defaultdict(lambda: defaultdict(list))
            by_haplotype = defaultdict(list)
            for path, hap in groups.items():
                sequence = path.split('#', 2)[2]
                by_chromosome[chromosome_label(sequence)][hap].append(path)
                by_haplotype[hap].append(path)
            rows = []
            summaries = []

            def emit(scope, hap, chromosome, counts):
                total = int(counts.sum())
                rows.extend((scope, hap, chromosome, n, k, int(bp), total,
                             100 * int(bp) / total if total else '') for k, bp in enumerate(counts) if k > 0)
                classified = dict.fromkeys(('core', 'softcore', 'shell', 'cloud'), 0)
                for k in range(1, n + 1):
                    category = next((label for label, cut in zip(('core', 'softcore', 'shell'), cuts) if k >= cut), 'cloud')
                    classified[category] += int(counts[k])
                classified['private'] = int(counts[1])
                summaries.extend((scope, hap, chromosome, n, label, bp, total,
                                  100 * bp / total if total else '', *cuts) for label, bp in classified.items())

            for cstamp, (chromosome, members) in enumerate(sorted(by_chromosome.items()), 1):
                chrom_counts = np.zeros(n + 1, dtype=np.uint64)
                for hap, paths in sorted(members.items()):
                    stamp += 1
                    counts = np.zeros(n + 1, dtype=np.uint64)
                    for index in walk_indices(paths):
                        fresh = index[seen[index] != stamp]
                        seen[fresh] = stamp
                        counts += bins(fresh)
                        fresh_chrom = index[chromosome_seen[index] != cstamp]
                        chromosome_seen[fresh_chrom] = cstamp
                        chrom_counts += bins(fresh_chrom)
                    emit('chromosome_haplotype', hap, chromosome, counts)
                emit('chromosome', '', chromosome, chrom_counts)
            reconstructed = np.zeros(n + 1, dtype=np.uint64)
            for hap, paths in sorted(by_haplotype.items()):
                stamp += 1
                counts = np.zeros(n + 1, dtype=np.uint64)
                for index in walk_indices(paths):
                    fresh = index[seen[index] != stamp]
                    seen[fresh] = stamp
                    counts += bins(fresh)
                reconstructed += counts
                emit('haplotype', hap, '', counts)
            if any(int(reconstructed[k]) != k * int(global_hist[k]) for k in range(1, n + 1)):
                raise ValueError('Haplotype attribution fails global coverage reconciliation')
    with (out / 'regional_coverage.tsv').open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(('scope', 'haplotype', 'chromosome', 'denominator', 'coverage', 'bp', 'unit_distinct_bp', 'percent'))
        writer.writerows(rows)
    with (out / 'regional_sharing.tsv').open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(('scope', 'haplotype', 'chromosome', 'denominator', 'class', 'bp', 'unit_distinct_bp', 'percent', 'core_min', 'softcore_min', 'shell_min'))
        writer.writerows(summaries)
    (out / 'regional_audit.json').write_text(json.dumps({
        'status': 'PASS', 'haplotypes': n, 'nodes': len(node_ids), 'global_histogram_bp': list(map(int, global_hist)),
        'chromosome_policy': 'chrN_part -> chrN; composite labels kept separate; all other paths -> unplaced',
        'additivity': 'Chromosome unions can overlap; haplotype unions deduplicate across all chromosomes.',
    }, indent=2) + '\n')
    (out / 'regional_report.md').write_text(
        '# CLIP regional sharing\n\n'
        '<details><summary>Haplotype and chromosome attribution</summary>\n\n'
        f'Coverage uses all {n} haplotypes. Core/softcore/shell minimum counts: {cuts}.\n\n'
        '[Sharing categories](regional_sharing.tsv) · [Exact coverage bins](regional_coverage.tsv) · [Audit](regional_audit.json)\n\n'
        'Each reporting unit counts distinct graph sequence once. Repeated traversals are not copy number. '
        'Chromosome unions can overlap and must not be summed as an additive genome partition. '
        'Composite scaffolds remain separate; unplaced includes paths without a harmonized chromosome name. '
        'Absent chromosome/haplotype combinations have no row, not an inferred biological zero. '
        'Private is a separate overlapping statistic. Both global histogram and haplotype attribution reconciliations passed.\n\n'
        '</details>\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('gfa', 'totals', 'groups', 'histogram', 'out'):
        parser.add_argument('--' + name, required=True)
    for name, default in [('core', 1), ('softcore', -1), ('shell', 2)]:
        parser.add_argument('--' + name, type=float, default=default)
    args = parser.parse_args()
    aggregate(args.gfa, args.totals, args.groups, args.histogram, args.out, args.core, args.softcore, args.shell)


if __name__ == '__main__':
    main()
