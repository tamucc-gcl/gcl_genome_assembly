"""Presentation adapter for retained junction assays; no alignments or decisions."""
from pathlib import Path
from misassembly_core import table, write_table
from misassembly_junction import credible_bridge


def export(packet, event):
    packet = Path(packet)
    evidence = event.get('junction_evidence', {})
    if evidence.get('hifi_status') != 'measured':
        return False
    root = packet / 'junction'
    original = {r['site']: r for r in table(root/'hifi/original_bridges.tsv')
                if r['event'] == event['id'] and r['anchor_bp'] == '1000'}
    competitive = table(root/'hifi/competitive_bridges.tsv')
    rows = []
    for site in evidence['sites']:
        old = original.get(site['id'], {})
        reads = [r for r in competitive if r['event'] == event['id'] and r['site'] == site['id']]
        applicable = site['credible_bridges'] is not None
        rows.append(dict(site=site['id'], point=site['point'],
            primary_all=old.get('mapq0', ''), primary_high=old.get('mapq30', ''),
            credible=site['credible_bridges'] if applicable else '',
            unique_arms=site['unique_arm_bridges'] if applicable else '',
            competing=sum(int(r['same_haplotype_alternatives']) > 0 for r in reads) if applicable else '',
            other_haplotype=sum(int(r['other_haplotype_alternatives']) > 0 for r in reads) if applicable else '',
            minimum_anchor_span_bp=site['minimum_anchor_span_bp'] or '',
            local_reads_long_enough=site['local_reads_long_enough'] if site['local_reads_long_enough'] is not None else ''))
        # Catch stale/incompatible evidence rather than silently plotting new counts.
        if applicable and sum(credible_bridge(r) for r in reads) != site['credible_bridges']:
            raise ValueError('Retained credible bridge count disagrees at '+site['id'])
    if not rows:
        return False
    write_table(root/(event['id']+'.placement_counts.tsv'), rows, list(rows[0]))
    return True
