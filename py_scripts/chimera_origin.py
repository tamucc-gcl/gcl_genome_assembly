#!/usr/bin/env python3
"""Pure helpers for read-only join-origin diagnostics. Coordinates: 0-based half-open."""
import csv
import hashlib
import re
from collections import defaultdict


def table(path):
    with open(path) as handle:
        reader = csv.DictReader((s for s in handle if s.strip() and not s.startswith('#')),
                                delimiter='\t')
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError(f'Missing or duplicate table header: {path}')
        result = list(reader)
        if any(None in r or None in r.values() for r in result):
            raise ValueError(f'Malformed table row: {path}')
        return result


def write_table(path, fields, rows):
    with open(path, 'w') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t',
                                lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def sha(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def revcomp(sequence):
    return sequence.upper().translate(str.maketrans('ACGTRYMKSWBDHVN', 'TGCAYRKMSWVHDBN'))[::-1]


def agp_rows(path):
    """Reject malformed/unknown-orientation AGPs instead of guessing a projection."""
    rows, ends, parts = [], {}, {}
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            f = line.rstrip().split('\t')
            if len(f) != 9:
                raise ValueError(f'AGP needs nine fields: {path}')
            obj, lo, hi, part = f[0], int(f[1]) - 1, int(f[2]), int(f[3])
            if lo != ends.get(obj, 0) or hi <= lo or part != parts.get(obj, 0) + 1:
                raise ValueError(f'Noncontiguous AGP: {path}: {obj}')
            ends[obj], parts[obj] = hi, part
            row = dict(object=obj, start=lo, end=hi, part=part, kind=f[4],
                       component='.', component_start='.', component_end='.', orientation='.')
            if f[4] in ('N', 'U'):
                if int(f[5]) != hi-lo:
                    raise ValueError('AGP gap length mismatch')
            else:
                start, end = int(f[6])-1, int(f[7])
                if start < 0 or end-start != hi-lo or f[8] not in ('+', '-'):
                    raise ValueError('AGP component has invalid length or unresolved orientation')
                row.update(component=f[5], component_start=start,
                           component_end=end, orientation=f[8])
            rows.append(row)
    return rows


def component_interval(row, lo, hi):
    if row['kind'] in ('N', 'U') or not row['start'] <= lo < hi <= row['end']:
        raise ValueError('Interval is not contained in one AGP component')
    if row['orientation'] == '+':
        return row['component_start'] + lo-row['start'], row['component_start'] + hi-row['start']
    return row['component_end'] - (hi-row['start']), row['component_end'] - (lo-row['start'])


def paf_rows(path):
    result = []
    with open(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            f = line.rstrip().split('\t')
            if len(f) < 12 or f[4] not in ('+', '-'):
                raise ValueError('Malformed PAF')
            qlen, qs, qe = map(int, f[1:4])
            tlen, ts, te, matches, block, mq = map(int, f[6:12])
            if not (0 <= qs < qe <= qlen and 0 <= ts < te <= tlen and block > 0):
                raise ValueError('Invalid PAF bounds')
            tags = {s.split(':', 2)[0]: s.split(':', 2)[2] for s in f[12:] if s.count(':') >= 2}
            result.append(dict(query=f[0], query_length=qlen, query_start=qs, query_end=qe,
                               strand=f[4], target=f[5], target_length=tlen,
                               target_start=ts, target_end=te, matches=matches,
                               block=block, mapq=mq, cigar=tags.get('cg', '.'),
                               alignment_type=tags.get('tp', '.')))
    return result


def exact_projection(hit, lo, hi):
    """Project only through a contiguous = CIGAR run; never bridge indels/edits.

    PAF CIGAR advances along target and along the oriented query (reverse query
    coordinates decrease). This function is intentionally more conservative than
    an approximate lift. Endpoints alone cannot establish continuity.
    """
    cigar = hit['cigar']
    ops = re.findall(r'(\d+)([=XMID])', cigar)
    if ''.join(n+op for n, op in ops) != cigar:
        return None
    q = hit['query_start'] if hit['strand'] == '+' else hit['query_end']
    t = hit['target_start']
    answer = None
    for count, op in ops:
        n = int(count)
        if n <= 0:
            raise ValueError('Zero-length CIGAR operation')
        if op == '=':
            if hit['strand'] == '+' and q <= lo < hi <= q+n:
                answer = (t+lo-q, t+hi-q)
            elif hit['strand'] == '-' and q-n <= lo < hi <= q:
                answer = (t+q-hi, t+q-lo)
        if op in '=XMI':
            q += n if hit['strand'] == '+' else -n
        if op in '=XMD':
            t += n
    expected_q = hit['query_end'] if hit['strand'] == '+' else hit['query_start']
    if q != expected_q or t != hit['target_end']:
        raise ValueError('CIGAR disagrees with PAF spans')
    return answer


def exact_span(hit, sequence, fetch):
    """Verify actual letters, including Ns, for a full ungapped-size placement.

    Exact sequence presence is historical continuity of this representation,
    not biological continuity across Ns, and not proof of exhaustive uniqueness.
    """
    if (hit['query_start'] != 0 or hit['query_end'] != len(sequence)
            or hit['query_length'] != len(sequence)
            or hit['target_end']-hit['target_start'] != len(sequence)):
        return False
    target = fetch(hit['target'], hit['target_start'], hit['target_end']).upper()
    return sequence.upper() == (target if hit['strand'] == '+' else revcomp(target))


def source_to_query(hit, lo, hi):
    """Use only AFTER full-span literal sequence verification."""
    if not hit['target_start'] <= lo <= hi <= hit['target_end']:
        raise ValueError('Source interval outside verified window')
    if hit['strand'] == '+':
        return lo-hit['target_start'], hi-hit['target_start']
    return hit['target_end']-hi, hit['target_end']-lo


def read_support(lines, scaffold, lo, hi, max_indel=50, min_mapq=20, max_nm=0.02):
    """Distinct molecules spanning two specified anchors, not midpoint coverage.

    Input may be a region-limited SAM: read-length observations are explicitly
    local and cannot estimate expected bridge counts or genome-wide alternatives.
    """
    spans, chains, eligible_lengths, observed = set(), set(), {}, set()
    for line in lines:
        if line.startswith('@') or not line.strip():
            continue
        f = line.rstrip().split('\t')
        if len(f) < 11:
            raise ValueError('Malformed SAM')
        flag = int(f[1])
        if flag & 0xF04 or f[2] != scaffold:
            continue
        name, pos = f[0], int(f[3])-1
        ops = re.findall(r'(\d+)([MIDNSHP=X])', f[5])
        if not ops or ''.join(n+op for n, op in ops) != f[5]:
            raise ValueError('Malformed SAM CIGAR')
        observed.add(name)
        eligible_lengths[name] = max(eligible_lengths.get(name, 0),
                                    sum(int(n) for n, op in ops if op in 'MIS=X'))
        start, chain_start, blocks, columns = pos, None, [], 0
        for n, op in ops:
            n = int(n)
            if op in 'M=X':
                if chain_start is None:
                    chain_start = pos
                pos += n
            elif op in 'IDN':
                if op == 'N' or n > max_indel:
                    if chain_start is not None:
                        blocks.append((chain_start, pos))
                    chain_start = None
                if op in 'DN':
                    pos += n
            if op in 'M=XID':
                columns += n
        if chain_start is not None:
            blocks.append((chain_start, pos))
        if start <= lo and pos >= hi:
            spans.add(name)
        tags = {s.split(':', 2)[0]: s.split(':', 2)[2] for s in f[11:] if s.count(':') >= 2}
        rate = int(tags['NM']) / columns if 'NM' in tags and columns else None
        if (min_mapq <= int(f[4]) < 255 and rate is not None and rate <= max_nm
                and any(a <= lo and b >= hi for a, b in blocks)):
            chains.add(name)
    return dict(primary_molecules_observed=len(observed),
                molecules_bracketing_anchors=len(spans),
                molecules_passing_chain_screen=len(chains),
                local_read_length_ge_anchor_span=sum(v >= hi-lo for v in eligible_lengths.values()),
                longest_local_primary_read=max(eligible_lengths.values(), default=0),
                expected_bridge_count='not_estimated',
                support_interpretation='descriptive_not_unique_anchor_or_fusion_validation')
