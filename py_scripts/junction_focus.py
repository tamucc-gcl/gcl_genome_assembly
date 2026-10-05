"""Coordinate helpers for focused, read-only junction diagnostics."""
import re


def union(intervals):
    result = []
    for lo, hi in sorted(intervals):
        if hi <= lo:
            continue
        if result and lo <= result[-1][1]:
            result[-1] = (result[-1][0], max(hi, result[-1][1]))
        else:
            result.append((lo, hi))
    return result


def intersection_length(a, b):
    return sum(max(0, min(y, v)-max(x, u)) for x, y in union(a) for u, v in union(b))


def read_blocks(pos, cigar, reverse=False):
    """Aligned reference/read blocks in original molecule coordinates, including H clips."""
    ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    if not ops or ''.join(n+op for n, op in ops) != cigar or any(int(n) <= 0 for n, _ in ops):
        raise ValueError('Invalid CIGAR')
    length = sum(int(n) for n, op in ops if op in 'MISH=X')
    q, blocks = 0, []
    for n, op in ops:
        n = int(n)
        if op in 'M=X':
            blocks.append((pos, pos+n, q, q+n))
        if op in 'MDN=X':
            pos += n
        if op in 'MISH=X':
            q += n
    if reverse:
        blocks = [(a,b,length-d,length-c) for a,b,c,d in blocks]
    return length, blocks


def focal_query(blocks, lo, hi, reverse=False):
    result = []
    for a,b,c,d in blocks:
        x,y = max(lo,a),min(hi,b)
        if x < y:
            result.append((d-(y-a),d-(x-a)) if reverse else (c+x-a,c+y-a))
    return union(result)


def anchor_intervals(start, end, length, size, offset):
    """Only complete windows: clipped windows are not comparable anchors."""
    left, right = (start-offset-size,start-offset), (end+offset,end+offset+size)
    return {side: (lo,hi) if 0 <= lo < hi <= length else None
            for side,(lo,hi) in [('left',left),('right',right)]}


def anchor_pair(left, right):
    if left['target'] != right['target']:
        return 'different_targets', None
    if left['strand'] != right['strand']:
        return 'orientation_discordant', None
    distance = (right['target_start']-left['target_end'] if left['strand']=='+'
                else left['target_start']-right['target_end'])
    return ('same_target_ordered' if distance >= 0 else 'overlap_or_reordered'), distance
