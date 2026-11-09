#!/usr/bin/env python3
"""Prepare source-bound contamination review sheets and apply reviewed removals.

Standalone tool, deliberately independent of taxonomic classification and plotting.
Coordinates are 0-based, half-open on the exact input FASTA. Internal removals
split retained flanks; they never manufacture an adjacency.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import tempfile


SCHEMA = 'contamination-review-v1'
FIELDS = ['action_id', 'assembly', 'stage', 'source_sha256', 'sequence',
          'action', 'start0', 'end0', 'selected', 'reviewer', 'reason']
LIFT_FIELDS = ['sequence', 'source_start0', 'source_end0', 'disposition',
               'output_sequence', 'output_start0', 'output_end0', 'sha256', 'action_id']


def digest(path):
    with open(path, 'rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def records(path):
    """Read one sequence at a time; preserve base case and complete headers."""
    seen, header, parts = set(), None, []
    with open(path, encoding='ascii') as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    if not parts:
                        raise ValueError('Empty FASTA sequence: ' + header)
                    yield header, ''.join(parts)
                header, parts = line[1:], []
                if not header or not header.split():
                    raise ValueError('Empty FASTA header')
                name = header.split()[0]
                if name in seen:
                    raise ValueError('Duplicate FASTA identifier: ' + name)
                seen.add(name)
            else:
                if header is None or not re.fullmatch('[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+', line):
                    raise ValueError('Invalid DNA FASTA sequence or missing header')
                parts.append(line)
    if header is None or not parts:
        raise ValueError('Empty FASTA or final sequence')
    yield header, ''.join(parts)


def write_table(path, rows, fields):
    with open(path, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def read_actions(path):
    with open(path, newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader((line for line in handle if line.strip() and not line.lstrip().startswith('#')), delimiter='\t')
        columns = reader.fieldnames or []
        if len(columns) != len(set(columns)) or not set(FIELDS) <= set(columns):
            raise ValueError('Review sheet requires unique columns: ' + ', '.join(FIELDS))
        result = list(reader)
        if any(None in row or None in row.values() for row in result):
            raise ValueError('Malformed review row')
        return result


def number(value):
    if not re.fullmatch('[0-9]+', str(value)):
        raise ValueError('Coordinates must be non-negative integers')
    return int(value)


def validate_actions(rows, assembly, stage, checksum, lengths):
    selected, ids = {}, set()
    for row in rows:
        if row['selected'] not in ('YES', 'NO') or not row['assembly']:
            raise ValueError('Every row needs assembly and selected=YES or NO')
        if row['assembly'] != assembly:
            continue
        if row['stage'] != stage or row['source_sha256'] != checksum:
            raise ValueError('Review sheet does not match source FASTA/stage')
        if not row['action_id'].strip() or row['action_id'] in ids:
            raise ValueError('Missing or duplicate action_id within assembly')
        ids.add(row['action_id'])
        if row['sequence'] not in lengths:
            raise ValueError('Unknown source sequence: ' + row['sequence'])
        if row['selected'] == 'NO':
            continue
        if not row['reviewer'].strip() or not row['reason'].strip():
            raise ValueError('Selected removal requires reviewer and reason')
        if row['action'] not in ('REMOVE_SEQUENCE', 'REMOVE_INTERVAL'):
            raise ValueError('Unknown selected action')
        lo, hi = number(row['start0']), number(row['end0'])
        if not 0 <= lo < hi <= lengths[row['sequence']]:
            raise ValueError('Removal outside source sequence')
        if row['action'] == 'REMOVE_SEQUENCE' and (lo, hi) != (0, lengths[row['sequence']]):
            raise ValueError('REMOVE_SEQUENCE must span the complete sequence')
        selected.setdefault(row['sequence'], []).append((lo, hi, row['action_id']))
    for name, intervals in selected.items():
        intervals.sort()
        if any(right[0] < left[1] for left, right in zip(intervals, intervals[1:])):
            raise ValueError('Overlapping selected removals: ' + name)
    return selected


def partitions(length, intervals):
    cursor = 0
    for lo, hi, action_id in intervals:
        if cursor < lo:
            yield cursor, lo, 'retained', ''
        yield lo, hi, 'removed', action_id
        cursor = hi
    if cursor < length:
        yield cursor, length, 'retained', ''


def write_record(handle, header, sequence):
    handle.write('>' + header + '\n')
    for offset in range(0, len(sequence), 80):
        handle.write(sequence[offset:offset+80] + '\n')


def prepare(args, output):
    checksum, metrics, template = digest(args.fasta), [], []
    for index, (header, sequence) in enumerate(records(args.fasta), 1):
        name, upper = header.split()[0], sequence.upper()
        acgt = sum(upper.count(base) for base in 'ACGT')
        metrics.append(dict(sequence=name, length_bp=len(sequence), acgt_bp=acgt,
                            n_bp=upper.count('N'), other_ambiguous_bp=len(sequence)-acgt-upper.count('N'),
                            gc_fraction=(upper.count('G')+upper.count('C'))/acgt if acgt else 'NA'))
        template.append(dict(zip(FIELDS, [f'A{index:05d}', args.assembly, args.stage, checksum,
                            name, 'REMOVE_SEQUENCE', 0, len(sequence), 'NO', '', ''])))
    if digest(args.fasta) != checksum:
        raise ValueError('Input FASTA changed during preparation')
    write_table(output/'sequence_metrics.tsv', metrics, list(metrics[0]))
    write_table(output/'review.tsv', template, FIELDS)
    (output/'source.json').write_text(json.dumps(dict(schema=SCHEMA, assembly=args.assembly,
        stage=args.stage, source_sha256=checksum, coordinates='0-based half-open',
        sequences=len(metrics), source_bp=sum(row['length_bp'] for row in metrics)), indent=2)+'\n')


def search_input(args, output):
    """Make bounded, non-overlapping coordinate windows and parallel search chunks.

    This only prepares primary protein searches. It makes no taxonomy or removal
    decisions; the evidence stage will select nucleotide follow-up separately.
    """
    if args.window_bp < 1 or args.chunk_bp < args.window_bp or args.min_acgt_bp < 1:
        raise ValueError('Require positive window/minimum ACGT and chunk >= window')
    checksum = digest(args.fasta)
    window_fields = ['query', 'sequence', 'start0', 'end0', 'acgt_bp', 'status', 'chunk']
    chunk_fields = ['chunk', 'queries', 'query_bp']
    chunks, handle, chunk_name, chunk_bp, chunk_queries, total_bp = [], None, '', 0, 0, 0
    index = 0
    try:
        with open(output/'windows.tsv', 'w', newline='') as table:
            writer = csv.DictWriter(table, fieldnames=window_fields, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            for header, sequence in records(args.fasta):
                name = header.split()[0]
                total_bp += len(sequence)
                for lo in range(0, len(sequence), args.window_bp):
                    index += 1
                    fragment = sequence[lo:lo+args.window_bp]
                    acgt = sum(fragment.upper().count(base) for base in 'ACGT')
                    query = f'Q{index:09d}'
                    eligible = acgt >= args.min_acgt_bp
                    if eligible:
                        if handle is None or chunk_bp+len(fragment) > args.chunk_bp:
                            if handle is not None:
                                handle.close()
                                chunks.append(dict(chunk=chunk_name, queries=chunk_queries, query_bp=chunk_bp))
                            chunk_name = f'protein_{len(chunks)+1:05d}.fasta'
                            handle = open(output/chunk_name, 'w', encoding='ascii')
                            chunk_bp, chunk_queries = 0, 0
                        write_record(handle, query, fragment)
                        chunk_bp += len(fragment)
                        chunk_queries += 1
                    writer.writerow(dict(zip(window_fields, [query, name, lo, lo+len(fragment), acgt,
                        'search' if eligible else 'insufficient_acgt', chunk_name if eligible else ''])))
        if handle is not None:
            chunks.append(dict(chunk=chunk_name, queries=chunk_queries, query_bp=chunk_bp))
    finally:
        if handle is not None:
            handle.close()
    if digest(args.fasta) != checksum:
        raise ValueError('Input FASTA changed during window preparation')
    write_table(output/'chunks.tsv', chunks, chunk_fields)
    (output/'source.json').write_text(json.dumps(dict(schema=SCHEMA, assembly=args.assembly,
        stage=args.stage, source_sha256=checksum, coordinates='0-based half-open',
        window_bp=args.window_bp, chunk_bp=args.chunk_bp, minimum_acgt_bp=args.min_acgt_bp,
        source_bp=total_bp, windows=index, search_queries=sum(row['queries'] for row in chunks),
        search_query_bp=sum(row['query_bp'] for row in chunks), chunks=len(chunks),
        boundary_policy='non-overlapping fixed windows; no inference across a window boundary'), indent=2)+'\n')


def apply(args, output):
    checksum = digest(args.fasta)
    lengths = {header.split()[0]:len(sequence) for header, sequence in records(args.fasta)}
    actions = read_actions(args.review)
    chosen = validate_actions(actions, args.assembly, args.stage, checksum, lengths)
    if set(chosen) == set(lengths) and all(sum(hi-lo for lo, hi, _ in chosen[name]) == length
                                         for name, length in lengths.items()):
        raise ValueError('Refusing to remove the entire assembly')
    lifts, totals, generated = [], dict(retained=0, removed=0), set()
    with open(output/'retained.fasta', 'w', encoding='ascii') as retained, open(output/'removed.fasta', 'w', encoding='ascii') as removed:
        for header, sequence in records(args.fasta):
            name = header.split()[0]
            for lo, hi, disposition, action_id in partitions(len(sequence), chosen.get(name, [])):
                unchanged = disposition == 'retained' and name not in chosen
                piece_name = name if unchanged else f'{name}__{disposition}_{lo}_{hi}'
                if piece_name in generated or (not unchanged and piece_name in lengths):
                    raise ValueError('Output identifier collision: ' + piece_name)
                generated.add(piece_name)
                fragment = sequence[lo:hi]
                write_record(retained if disposition == 'retained' else removed,
                             header if unchanged else piece_name, fragment)
                totals[disposition] += len(fragment)
                lifts.append(dict(zip(LIFT_FIELDS, [name, lo, hi, disposition, piece_name,
                                    0, hi-lo, hashlib.sha256(fragment.encode('ascii')).hexdigest(), action_id])))
    if not chosen:
        shutil.copyfile(args.fasta, output/'retained.fasta')
    if sum(totals.values()) != sum(lengths.values()) or digest(args.fasta) != checksum:
        raise ValueError('Source changed or sequence accounting failed')
    write_table(output/'coordinate_map.tsv', lifts, LIFT_FIELDS)
    write_table(output/'applied_actions.tsv', [row for row in actions if row['assembly'] == args.assembly and row['selected'] == 'YES'], list(actions[0]) if actions else FIELDS)
    (output/'verification.json').write_text(json.dumps(dict(schema=SCHEMA, assembly=args.assembly,
        stage=args.stage, source_sha256=checksum, review_sha256=digest(args.review),
        coordinates='0-based half-open', source_bp=sum(lengths.values()),
        retained_bp=totals['retained'], removed_bp=totals['removed'],
        selected_actions=sum(len(value) for value in chosen.values()),
        retained_sha256=digest(output/'retained.fasta'), removed_sha256=digest(output/'removed.fasta')),
        indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    for name in ('prepare', 'apply', 'search-input'):
        sub = commands.add_parser(name)
        sub.add_argument('--fasta', type=Path, required=True)
        sub.add_argument('--assembly', required=True)
        sub.add_argument('--stage', required=True)
        sub.add_argument('--out', type=Path, required=True)
        if name == 'apply':
            sub.add_argument('--review', type=Path, required=True)
        if name == 'search-input':
            sub.add_argument('--window-bp', type=int, default=100000)
            sub.add_argument('--chunk-bp', type=int, default=5000000)
            sub.add_argument('--min-acgt-bp', type=int, default=100)
    args = parser.parse_args()
    if args.out.exists() or args.out.is_symlink():
        parser.error('Output already exists; choose a new directory')
    if not args.assembly.strip() or not args.stage.strip():
        parser.error('Assembly and stage must be non-empty')
    args.out.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix='.contamination-', dir=args.out.parent))
    try:
        {'prepare':prepare, 'apply':apply, 'search-input':search_input}[args.command](args, temporary)
        os.rename(temporary, args.out)
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)


if __name__ == '__main__':
    main()
