#!/usr/bin/env python3
"""Auditable sampling and nucleotide routing for optional post-FCS review.

Sampling nominates evidence, never removal. Coordinates remain on the exact
post-FCS FASTA. Unselected windows remain explicitly unsearched in reports.
"""
import argparse
import csv
import json
from pathlib import Path
import re

from contamination_review import digest, records, write_record, write_table
from contamination_evidence import table, search_hits, grouped_hits, relationship


def coding_score(sequence):
    """Longest stop-free run in six frames; a sampling heuristic, not a gene call."""
    sequence = sequence.upper()
    best = 0
    for strand in (sequence, sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]):
        # Ambiguity breaks a run, rather than making N-rich regions look coding-rich.
        for fragment in re.findall('[ACGT]+', strand):
            stops = [[], [], []]
            for match in re.finditer(r'(?=TAA|TAG|TGA)', fragment):
                stops[match.start() % 3].append(match.start())
            for frame, positions in enumerate(stops):
                previous = frame - 3
                for stop in positions + [len(fragment) - (len(fragment) - frame) % 3]:
                    best = max(best, (stop - previous - 3) // 3)
                    previous = stop
    return best


def depth_outliers(path, fold):
    if not path or not path.stat().st_size:
        return set(), None
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    values = sorted((float(row['meandepth']), int(row['endpos'])) for row in rows)
    total = sum(length for _, length in values)
    cumulative, baseline = 0, None
    for depth, length in values:
        cumulative += length
        if cumulative >= total / 2:
            baseline = depth
            break
    if baseline is None or baseline <= 0:
        return set(), baseline
    return {row['#rname'] for row in rows
            if float(row['meandepth']) <= baseline / fold or
            float(row['meandepth']) >= baseline * fold}, baseline


def write_queries(fasta, rows, output, chunk_bp):
    """Write only nominated windows, preserving source coordinates and IDs."""
    by_sequence = {}
    for row in rows:
        if row['status'] == 'search':
            by_sequence.setdefault(row['sequence'], []).append(row)
    chunks, handle, size, count, name = [], None, 0, 0, ''
    try:
        for header, sequence in records(fasta):
            for row in by_sequence.get(header.split()[0], []):
                fragment = sequence[int(row['start0']):int(row['end0'])]
                if handle is None or size + len(fragment) > chunk_bp:
                    if handle is not None:
                        handle.close()
                        chunks.append(dict(chunk=name, queries=count, query_bp=size))
                    name = f'protein_{len(chunks) + 1:05d}.fasta'
                    handle = (output / name).open('w', encoding='ascii')
                    size, count = 0, 0
                write_record(handle, row['query'], fragment)
                row['chunk'] = name
                size += len(fragment)
                count += 1
        if handle is not None:
            chunks.append(dict(chunk=name, queries=count, query_bp=size))
    finally:
        if handle is not None:
            handle.close()
    write_table(output / 'chunks.tsv', chunks, ['chunk', 'queries', 'query_bp'])
    write_table(output / 'windows.tsv', rows,
                ['query', 'sequence', 'start0', 'end0', 'acgt_bp', 'status', 'chunk', 'selection_reason'])
    return chunks


def sample(args):
    if not (100 <= args.window_bp <= args.spacing_bp and args.chunk_bp >= args.window_bp and
            args.short_contig_bp >= 0 and 0 < args.gc_delta <= 1 and args.depth_fold > 1):
        raise ValueError('Invalid staged sampling sizes or outlier thresholds')
    source = json.loads((args.metrics / 'source.json').read_text())
    if digest(args.fasta) != source['source_sha256']:
        raise ValueError('Sampling FASTA differs from metrics source')
    metrics = table(args.metrics / 'sequence_metrics.tsv')
    acgt = sum(int(row['acgt_bp']) for row in metrics)
    gc = sum(float(row['gc_fraction']) * int(row['acgt_bp'])
             for row in metrics if row['gc_fraction'] != 'NA') / acgt if acgt else None
    depth, baseline = depth_outliers(args.coverage, args.depth_fold)
    rows = []
    for header, sequence in records(args.fasta):
        name = header.split()[0]
        local, blocks = [], {}
        for lo in range(0, len(sequence), args.window_bp):
            fragment = sequence[lo:lo + args.window_bp]
            upper = fragment.upper()
            count = sum(upper.count(base) for base in 'ACGT')
            reasons = []
            if count >= 100:
                if len(sequence) <= args.short_contig_bp:
                    reasons.append('short-contig')
                if gc is not None and abs((upper.count('G') + upper.count('C')) / count - gc) >= args.gc_delta:
                    reasons.append('gc-outlier')
                if name in depth:
                    reasons.append('depth-outlier')
                score = (coding_score(fragment), count, -lo)
                block = lo // args.spacing_bp
                if block not in blocks or score > blocks[block][0]:
                    blocks[block] = (score, len(local))
            local.append(dict(query=f'Q{len(rows) + len(local) + 1:09d}', sequence=name,
                              start0=lo, end0=lo + len(fragment), acgt_bp=count,
                              status='insufficient_acgt' if count < 100 else 'unsearched',
                              chunk='', selection_reason=';'.join(reasons)))
        for _, index in blocks.values():
            row = local[index]
            row['selection_reason'] = ';'.join(filter(None, [row['selection_reason'], 'distributed-coding-sample']))
        for row in local:
            if row['selection_reason']:
                row['status'] = 'search'
        rows.extend(local)
    args.out.mkdir()
    chunks = write_queries(args.fasta, rows, args.out, args.chunk_bp)
    source.update(window_bp=args.window_bp, chunk_bp=args.chunk_bp, minimum_acgt_bp=100,
                  windows=len(rows), search_queries=sum(row['queries'] for row in chunks),
                  search_query_bp=sum(row['query_bp'] for row in chunks), chunks=len(chunks),
                  boundary_policy='non-overlapping fixed windows; no inference across a window boundary',
                  search_mode='staged', sampling=dict(spacing_bp=args.spacing_bp,
                      short_contig_bp=args.short_contig_bp, gc_delta=args.gc_delta,
                      genome_gc=gc, depth_fold=args.depth_fold, depth_baseline=baseline,
                      coverage_available=bool(args.coverage and args.coverage.stat().st_size),
                      representative='longest six-frame stop-free run per spacing block; earliest tie'))
    if digest(args.fasta) != source['source_sha256']:
        raise ValueError('FASTA changed during sampling')
    (args.out / 'source.json').write_text(json.dumps(source, indent=2) + '\n')


def route(args):
    rows = table(args.prepared / 'windows.tsv')
    windows = {row['query']: row for row in rows}
    hits = grouped_hits(search_hits(args.hits, windows, args))
    chosen = set()
    by_sequence = {}
    for index, row in enumerate(rows):
        by_sequence.setdefault(row['sequence'], []).append(index)
    for indices in by_sequence.values():
        for offset, index in enumerate(indices):
            row = rows[index]
            if row['status'] != 'search':
                continue
            relation = relationship(hits.get(row['query'], []))
            # A host hit does not veto composition/depth evidence.
            outlier = any(reason in row['selection_reason'] for reason in ('gc-outlier', 'depth-outlier'))
            if relation != 'host-compatible' or outlier:
                chosen.update(indices[max(0, offset - 1):offset + 2])
    nominated = []
    for index, row in enumerate(rows):
        if index in chosen and int(row['acgt_bp']) >= 100:
            nominated.append(dict(row, status='search', chunk='', selection_reason='unresolved-nonhost-or-outlier-with-flanks'))
    args.out.mkdir()
    source = json.loads((args.prepared / 'source.json').read_text())
    if digest(args.fasta) != source['source_sha256']:
        raise ValueError('Nucleotide routing FASTA differs from prepared source')
    chunks = write_queries(args.fasta, nominated, args.out, source['chunk_bp'])
    source.update(search_queries=len(nominated), search_query_bp=sum(row['query_bp'] for row in chunks), chunks=len(chunks))
    (args.out / 'source.json').write_text(json.dumps(source, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    for command in ('sample', 'route'):
        sub = commands.add_parser(command)
        sub.add_argument('--fasta', type=Path, required=True)
        sub.add_argument('--out', type=Path, required=True)
        if command == 'sample':
            sub.add_argument('--metrics', type=Path, required=True)
            sub.add_argument('--coverage', type=Path)
            sub.add_argument('--window-bp', type=int, default=10000)
            sub.add_argument('--spacing-bp', type=int, default=1000000)
            sub.add_argument('--chunk-bp', type=int, default=5000000)
            sub.add_argument('--short-contig-bp', type=int, default=100000)
            sub.add_argument('--gc-delta', type=float, default=.15)
            sub.add_argument('--depth-fold', type=float, default=4)
        else:
            sub.add_argument('--prepared', type=Path, required=True)
            sub.add_argument('--hits', nargs='*', type=Path, default=[])
            sub.add_argument('--taxdump', type=Path, required=True)
            sub.add_argument('--taxid', type=int, required=True)
            sub.add_argument('--host-rank', default='phylum')
    args = parser.parse_args()
    {'sample': sample, 'route': route}[args.command](args)


if __name__ == '__main__':
    main()
