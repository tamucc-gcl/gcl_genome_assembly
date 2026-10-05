"""Experimental all-junction review helpers; no automatic assembly edits."""
import re
from junction_focus import anchor_pair


def n_runs(sequence):
    return [(m.start(), m.end()) for m in re.finditer('N+', sequence.upper())]


def unique_anchor(hits, coverage=0.8, mapq=20):
    """Uniqueness among emitted near-full alignments, not proof of genomic uniqueness."""
    qualified = []
    for h in hits:
        cigar = h.get('cigar', '')
        aligned = sum(int(n) for n, op in re.findall(r'(\d+)([=XMID])', cigar) if op in '=XM')
        if aligned / h['query_length'] >= coverage:
            qualified.append(h)
    if len(qualified) != 1 or not mapq <= qualified[0]['mapq'] < 255:
        return None
    return qualified[0]


def compare_anchors(left, right):
    l, r = unique_anchor(left), unique_anchor(right)
    if l is None or r is None:
        return dict(status='ambiguous_or_missing', separation_bp='.')
    status, distance = anchor_pair(l, r)
    return dict(status=status, separation_bp=distance if distance is not None else '.',
                left_target=l['target'], right_target=r['target'],
                left_start=l['target_start'], left_end=l['target_end'],
                right_start=r['target_start'], right_end=r['target_end'],
                left_strand=l['strand'], right_strand=r['strand'],
                left_mapq=l['mapq'], right_mapq=r['mapq'])


def recommend(row, sister):
    """Triage, not calibrated fusion probabilities or executable break instructions."""
    usable = [r for r in sister if r.get('anchor_context') == 'clean']
    states = {r['status'] for r in usable}
    conflict = bool(states & {'different_targets', 'orientation_discordant', 'overlap_or_reordered'})
    ordered = 'same_target_ordered' in states
    chains = row.get('hifi_chain_molecules', '.')
    molecular = chains != '.' and int(chains) > 0
    if conflict:
        return ('REVIEW_CONFLICT' if molecular or ordered else 'PRIORITIZE_MISJOIN_REVIEW',
                'Sister has an alternative arrangement; check sister fragmentation, phase and both boundary contexts.',
                'Local uniquely anchored molecule review and opportunity-adjusted Hi-C boundary review')
    if molecular:
        return ('SUPPORT_PRESENT', 'Screened HiFi chains present; repeat uniqueness and phase are not certified.',
                'Review molecule anchors before assigning RETAIN')
    if ordered:
        return ('COMPARATIVE_SUPPORT_ONLY', 'Sister flanks are ordered on one scaffold; intervening structure may differ.',
                'Assess direct boundary support; do not infer continuity from flanking anchors alone')
    return ('NOT_RESOLVED', 'Missing or ambiguous evidence is not evidence against the join.',
            'Assess uniquely mappable flanks and informative read/contact opportunity')


def apply_reviews(rows, reviews):
    lookup = {r['id']: r for r in rows}
    seen = set()
    allowed = {'RETAIN', 'BREAK_PROBABLE_MISJOIN', 'UNJOIN_UNSUPPORTED', 'UNRESOLVED'}
    for review in reviews:
        ident = review['id']
        if ident not in lookup or ident in seen:
            raise ValueError('Unknown or duplicate reviewed junction')
        seen.add(ident)
        r = lookup[ident]
        for key in ('assembly', 'scaffold', 'start', 'end', 'assessment_sha256'):
            if str(review[key]) != str(r[key]):
                raise ValueError('Reviewed coordinates or assessment checksum changed')
        action = review['decision']
        if action not in allowed or not all(review.get(k, '').strip() for k in
                ('reviewer', 'reason', 'evidence_for', 'evidence_against')):
            raise ValueError('Decision needs reviewer, reason and evidence for/against')
        if action.startswith('BREAK') or action == 'UNJOIN_UNSUPPORTED':
            if r['kind'] != 'verified_agp_gap':
                raise ValueError('This diagnostic only accepts unjoining verified AGP gaps; internal cuts need separate localization')
        r.update({k: review[k] for k in ('decision', 'reviewer', 'reason', 'evidence_for', 'evidence_against')})
