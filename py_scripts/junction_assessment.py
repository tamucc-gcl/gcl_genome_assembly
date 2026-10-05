"""Read-only junction diagnostics. All internal coordinates are 0-based half-open."""
from collections import defaultdict, Counter
import hashlib
import re


def overlap(a, b, c, d):
    return a < d and c < b


def stable_id(assembly, scaffold, lo, hi):
    return 'gap_' + hashlib.sha256(f'{assembly}\t{scaffold}\t{lo}\t{hi}'.encode()).hexdigest()[:16]


def validate_registry(rows):
    ids = set()
    for r in rows:
        if not r.get('sample') or not re.fullmatch(r'[A-Za-z0-9_.-]+',r['assembly']) or r['assembly'] in ('.','..'):
            raise ValueError('Sample and safe assembly ID required')
        if not re.fullmatch(r'[A-Za-z0-9_.-]+', r['id']) or r['id'] in ids:
            raise ValueError('Invalid or duplicate tracking ID')
        ids.add(r['id'])
        if not re.fullmatch(r'[a-f0-9]{64}', r['assessment_sha256']):
            raise ValueError('Registry needs the exact assessment SHA256')
        if not 0 <= int(r['start']) < int(r['end']):
            raise ValueError('Invalid registry interval')
    return rows


def project_gap(row, length, reverse=False):
    return (length-row['end'], length-row['start']) if reverse else (row['start'], row['end'])


def choose_gaps(seeds, gaps, controls=3, nearby=1000000):
    """Select intersecting gaps and nearby controls; never snap a seed to a gap."""
    selected, links = {}, []
    suspect = set()
    for s in seeds:
        for g in gaps:
            if s['scaffold'] == g['scaffold'] and overlap(int(s['start']), int(s['end']), g['start'], g['end']):
                suspect.add(g['id'])
                selected[g['id']] = dict(g, role='candidate_gap')
                links.append(dict(seed=s['id'], gap=g['id'], relationship='intersects_seed_interval'))
    for s in seeds:
        center = (int(s['start'])+int(s['end']))/2
        candidates = [g for g in gaps if g['scaffold'] == s['scaffold'] and g['id'] not in suspect
                      and min(abs(g['start']-int(s['end'])), abs(g['end']-int(s['start']))) <= nearby]
        # Controls are comparison joins, not validated true joins. Prefer equal gap lengths.
        lengths = {g['end']-g['start'] for g in gaps if g['id'] in suspect and g['scaffold'] == s['scaffold']}
        candidates.sort(key=lambda g: (0 if g['end']-g['start'] in lengths else 1,
                                       abs((g['start']+g['end'])/2-center), g['id']))
        for g in candidates[:controls]:
            selected.setdefault(g['id'], dict(g, role='nearby_control'))
            links.append(dict(seed=s['id'], gap=g['id'], relationship='nearby_control_not_a_proposed_cut'))
        if not any(l['seed'] == s['id'] and l['relationship'] == 'intersects_seed_interval' for l in links):
            links.append(dict(seed=s['id'], gap='.', relationship='no_verified_gap_in_interval'))
    return list(selected.values()), links


class Lift:
    """Project a pairs end through one AGP; duplicated source placements are unavailable."""
    def __init__(self, rows):
        self.by_source = defaultdict(list)
        for r in rows:
            if r['kind'] not in ('N', 'U'):
                self.by_source[r['component']].append(r)

    def locate(self, name, pos):
        rows = [r for r in self.by_source[name] if r['component_start'] <= pos < r['component_end']]
        if len(rows) != 1:
            return None
        r = rows[0]
        p = r['start']+pos-r['component_start'] if r['orientation'] == '+' else r['end']-1-(pos-r['component_start'])
        return r['object'], p


def library_for(qname, manifests):
    values = [r['library_id'] for r in manifests if qname.startswith(r['qname_prefix'])]
    if len(values) > 1:
        raise ValueError('Read-set prefixes overlap')
    return values[0] if values else 'UNASSIGNED'


class Contacts:
    """Same-width flank contacts and alternative bins; descriptive, not a likelihood."""
    def __init__(self, gaps, lengths, flank=100000):
        self.flank, self.windows = flank, {}
        self.index = defaultdict(list)
        self.total, self.cross, self.margin, self.partners = Counter(), Counter(), Counter(), Counter()
        for g in gaps:
            for side, lo, hi in [('left', max(0, g['start']-flank), g['start']),
                                  ('right', g['end'], min(lengths[g['scaffold']], g['end']+flank))]:
                self.windows[g['id'], side] = (g['scaffold'], lo, hi)
                for b in range(lo//flank, (hi-1)//flank+1) if hi > lo else []:
                    self.index[g['scaffold'], b].append((g['id'], side, lo, hi))

    def ends(self, chrom, pos):
        return {(g, side) for g, side, lo, hi in self.index[chrom, pos//self.flank] if lo <= pos < hi}

    def add(self, library, x, y):
        self.total[library] += 1
        a, b = self.ends(*x), self.ends(*y)
        connected = {g for g, side in a if (g, 'right' if side == 'left' else 'left') in b}
        for g in connected:
            self.cross[library, g] += 1
        for tags, partner in ((a, y), (b, x)):
            for g, side in tags:
                # Exclude pairs entirely inside the originating flank from alternative counts.
                c, lo, hi = self.windows[g, side]
                if partner[0] == c and lo <= partner[1] < hi:
                    continue
                self.margin[library, g, side] += 1
                self.partners[library, g, side, partner[0], partner[1]//self.flank] += 1


def decision_rows(seeds, links):
    return [dict(id=s['id'], assembly=s['assembly'], scaffold=s['scaffold'],
                 interval_start=s['start'], interval_end=s['end'], assessment_sha256=s['assessment_sha256'],
                 verified_gaps=','.join(sorted({l['gap'] for l in links if l['seed'] == s['id'] and l['relationship'] == 'intersects_seed_interval'})) or '.',
                 interpretation='UNRESOLVED', action='REVIEW', chosen_gap='.',
                 evidence_for='', evidence_against='', missing_evidence='', rationale='', reviewer='', cut_authorized='no') for s in seeds]


def validate_decisions(decisions, seeds, gaps):
    seedmap, gapmap = {r['id']: r for r in seeds}, {r['id']: r for r in gaps}
    if len(decisions) != len(seedmap) or {r['id'] for r in decisions} != set(seedmap):
        raise ValueError('Decision table must contain each seed exactly once')
    for d in decisions:
        s = seedmap[d['id']]
        if any(str(d[k]) != str(s[k]) for k in ('assembly', 'scaffold', 'assessment_sha256')):
            raise ValueError('Decision identity/checksum mismatch')
        if int(d['interval_start']) != int(s['start']) or int(d['interval_end']) != int(s['end']):
            raise ValueError('Decision interval changed')
        if d['action'] not in ('REVIEW', 'KEEP', 'BREAK', 'DEFER'):
            raise ValueError('Unknown decision action')
        if d['action'] in ('KEEP', 'BREAK') and not all(d[k].strip() for k in ('rationale', 'reviewer', 'evidence_for', 'evidence_against')):
            raise ValueError('Reviewed decisions require rationale, reviewer and both evidence fields (use none observed if applicable)')
        if d['action'] == 'BREAK':
            g = gapmap.get(d['chosen_gap'])
            if not g or g['assembly'] != s['assembly'] or g['scaffold'] != s['scaffold'] or not overlap(int(s['start']), int(s['end']), int(g['start']), int(g['end'])):
                raise ValueError('BREAK requires a verified gap inside the seed interval')
    return True
