"""Advisory QC for native SyRI events; no calls removed or labelled validated.

Coordinates are converted from 1-based inclusive to 0-based half-open internally.
Direct AL children only: local CPG/CPL/TDM lack equivalent alignment support rows.
"""
import argparse
import bisect
import collections
import csv
import hashlib
import json
from pathlib import Path
import re

TYPES = {'INV', 'TRANS', 'INVTR', 'DUP', 'INVDP', 'CPG', 'CPL', 'TDM'}
FIELDS = 'ref rstart rend refseq altseq query qstart qend id parent type copy_status'.split()


def interval(start, end):
    a, b = int(start), int(end)
    if min(a, b) < 1:
        raise ValueError('Invalid nonpositive SyRI coordinate')
    return min(a, b) - 1, max(a, b)


def union_size(intervals):
    total, end = 0, -1
    for a, b in sorted(intervals):
        total += max(0, b - max(a, end))
        end = max(end, b)
    return total


def fasta_index(path):
    """Scan sequences without retaining bases; coalesce N runs across FASTA lines."""
    lengths, gaps = {}, collections.defaultdict(list)
    digest = hashlib.sha256()
    name = None
    with open(path, 'rb') as handle:
        for raw in handle:
            digest.update(raw)
            line = raw.decode('ascii').strip()
            if not line:
                continue
            if line.startswith('>'):
                name = line[1:].split()[0]
                if name in lengths:
                    raise ValueError('Duplicate FASTA name: ' + name)
                lengths[name] = 0
            else:
                if name is None:
                    raise ValueError('Sequence before FASTA header')
                seq = ''.join(line.split())
                offset = lengths[name]
                for match in re.finditer('[Nn]+', seq):
                    a, b = offset + match.start(), offset + match.end()
                    if gaps[name] and gaps[name][-1][1] == a:
                        gaps[name][-1] = (gaps[name][-1][0], b)
                    else:
                        gaps[name].append((a, b))
                lengths[name] += len(seq)
    if not lengths or any(v == 0 for v in lengths.values()):
        raise ValueError('Empty FASTA sequence')
    return lengths, gaps, digest.hexdigest()


def boundary_gap_distance(bounds, gaps):
    if not gaps:
        return None
    starts = [a for a, b in gaps]
    distances = []
    for p in bounds:
        k = bisect.bisect_right(starts, p)
        for a, b in gaps[max(0, k - 1):k + 1]:
            distances.append(max(a - p, p - b, 0))
    return min(distances)


def assess(event, children, refs, queries, near_bp, min_fraction):
    result = {k: event[k] for k in FIELDS if k not in ('refseq', 'altseq')}
    flags = []
    for role, name_key, start_key, end_key, index in (
            ('reference', 'ref', 'rstart', 'rend', refs),
            ('query', 'query', 'qstart', 'qend', queries)):
        lengths, gaps, _ = index
        name = event[name_key]
        if name not in lengths:
            raise ValueError('Call references unknown sequence: ' + name)
        a, b = interval(event[start_key], event[end_key])
        if b > lengths[name]:
            raise ValueError('Call outside FASTA: ' + event['id'])
        result[role + '_span_bp'] = b - a
        result[role + '_end_distance_bp'] = min(a, lengths[name] - b)
        gap_distance = boundary_gap_distance((a, b), gaps[name])
        result[role + '_boundary_N_distance_bp'] = gap_distance
        if result[role + '_end_distance_bp'] <= near_bp:
            flags.append(role + '_near_sequence_end')
        if gap_distance is not None and gap_distance <= near_bp:
            flags.append(role + '_boundary_near_N')
        clips = []
        for child in children:
            if child[name_key] != name:
                raise ValueError('Child chromosome mismatch: ' + child['id'])
            c, d = interval(child[start_key], child[end_key])
            if max(a, c) < min(b, d):
                clips.append((max(a, c), min(b, d)))
        fraction = union_size(clips) / (b - a) if children else None
        result[role + '_assigned_span_fraction'] = fraction
        if fraction is not None and fraction < min_fraction:
            flags.append(role + '_low_assigned_span_fraction')
    result['direct_alignment_children'] = len(children)
    if not children:
        flags.append('direct_alignment_metrics_unavailable')
    result['qc_status'] = 'REVIEW' if flags else 'NO_FLAGS_IN_IMPLEMENTED_CHECKS'
    result['qc_flags'] = ';'.join(flags) or 'none'
    result['read_support'] = 'NOT_ASSESSED'
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('source', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--near-bp', type=int, default=1000)
    p.add_argument('--min-assigned-fraction', type=float, default=0.5)
    a = p.parse_args()
    if a.near_bp < 0 or not 0 <= a.min_assigned_fraction <= 1:
        p.error('Invalid advisory threshold')
    a.output.mkdir(parents=True, exist_ok=True)
    events, children, seen = [], collections.defaultdict(list), set()
    with open(a.source / 'syri.out') as handle:
        for row in csv.reader(handle, delimiter='\t'):
            if not row or row[0].startswith('#'):
                continue
            if len(row) != 12:
                raise ValueError('Expected 12 SyRI columns')
            item = dict(zip(FIELDS, row))
            if item['type'] not in TYPES and item['type'] not in {'SYNAL', 'INVAL', 'TRANSAL', 'INVTRAL', 'DUPAL', 'INVDPAL'}:
                continue
            if item['id'] in seen:
                raise ValueError('Duplicate SyRI ID: ' + item['id'])
            seen.add(item['id'])
            if item['type'] in TYPES:
                events.append(item)
            elif item['type'].endswith('AL'):
                children[item['parent']].append(item)
    refs = fasta_index(a.source / 'reference.fa')
    queries = fasta_index(a.source / 'query.fa')
    results = [assess(e, children[e['id']], refs, queries, a.near_bp,
                      a.min_assigned_fraction) for e in events]
    # All events preserved; no allele strings duplicated into the small report.
    with open(a.output / 'sv_qc.tsv', 'w', newline='') as handle:
        fields = list(results[0]) if results else ['id', 'type', 'qc_status', 'qc_flags', 'read_support']
        w = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(results)
    counts = collections.Counter(r['qc_status'] for r in results)
    metadata = {'near_bp': a.near_bp, 'min_assigned_fraction': a.min_assigned_fraction,
                'reference_sha256': refs[2], 'query_sha256': queries[2],
                'events': len(events), 'status_counts': dict(counts),
                'source': str(a.source.resolve())}
    (a.output / 'sv_qc_settings.json').write_text(json.dumps(metadata, indent=2) + '\n')
    (a.output / 'sv_qc.md').write_text(
        '# Advisory SV QC\n\n' + f'{len(events)} native events retained. Status counts: {dict(counts)}.\n\n' +
        'These are configurable screening flags, not calibrated biological filters. '
        'NO_FLAGS_IN_IMPLEMENTED_CHECKS is not a high-confidence or validated classification.\n\n'
        'Checks: event-boundary proximity to Ns/sequence ends and the union of direct '
        'SyRI alignment-child spans clipped to the event, separately in both assemblies. '
        'Span coverage includes internal alignment gaps; it is not identity. '
        'CPG/CPL/TDM often lack direct alignment children and remain unassessed by this metric.\n\n'
        'Read support, repeat ambiguity, competing placements, alignment identity and '
        'breakpoint adjacency validation are NOT assessed in this first checkpoint. '
        'Original syri.out remains the authoritative unfiltered callset.\n')


if __name__ == '__main__':
    main()
