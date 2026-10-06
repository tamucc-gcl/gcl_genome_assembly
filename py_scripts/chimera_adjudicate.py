#!/usr/bin/env python3
"""Generate auditable decisions and source-bound actions, never using votes as permission."""
import argparse
import csv
import hashlib
import html
import json
from pathlib import Path

ACTION_FIELDS = ['id', 'assembly', 'coordinate_stage', 'assessment_sha256', 'scaffold', 'action',
                 'cut_bp', 'gap_start', 'gap_end', 'decision_source', 'evidence_packet_id', 'reason',
                 'auto_eligible', 'policy_version']


def decide(measurement):
    """Missing assays fail closed. Independent individuals, not haplotypes, are counted."""
    required = ['verified_gap', 'independent_discordant_individuals', 'hifi_informative',
                'hifi_spanning_molecules', 'hic_informative', 'hic_support_loss',
                'matched_controls_pass', 'alternative_placements_checked', 'graph_contradiction']
    if measurement.get('hifi_spanning_molecules', 0) >= 2 or measurement.get('graph_contradiction') is True:
        return 'RETAIN', False, 'local continuity contradicts a cut; chromosome-scale interpretation remains unresolved'
    missing = [key for key in required if key not in measurement]
    if missing:
        return 'UNRESOLVED', False, 'missing evidence: '+', '.join(missing)
    if measurement['hifi_spanning_molecules'] >= 2 or measurement['graph_contradiction'] is True:
        return 'RETAIN', False, 'local continuity contradicts a cut; chromosome-scale interpretation remains unresolved'
    checks = [measurement['verified_gap'] is True,
              measurement['independent_discordant_individuals'] >= 3,
              measurement['hifi_informative'] is True, measurement['hifi_spanning_molecules'] == 0,
              measurement['hic_informative'] is True, measurement['hic_support_loss'] is True,
              measurement['matched_controls_pass'] is True,
              measurement['alternative_placements_checked'] is True,
              measurement['graph_contradiction'] is False]
    if all(checks):
        return 'UNJOIN_UNSUPPORTED', True, 'verified gap; independent chromosome discordance and informative support loss without continuity contradiction'
    return 'UNRESOLVED', False, 'evidence does not satisfy conservative gap-v1 eligibility'


def write(path, rows, fields):
    with open(path, 'w') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--calls', required=True)
    p.add_argument('--assembly', required=True)
    p.add_argument('--context', required=True)
    p.add_argument('--out', required=True)
    a = p.parse_args()
    out = Path(a.out); out.mkdir()
    with open(a.calls) as handle:
        calls = list(csv.DictReader((x for x in handle if not x.startswith('#')), delimiter='\t'))
    (out/'assessed_candidates.tsv').write_bytes(Path(a.calls).read_bytes())
    evidence_path = Path(a.context)/'decision_measurements.json'
    measurements = json.loads(evidence_path.read_text()) if evidence_path.exists() else {}
    decisions, actions = [], []
    for index, row in enumerate(calls):
        key = hashlib.sha256(json.dumps(row, sort_keys=True).encode()).hexdigest()[:20]
        measurement = measurements.get(key, {})
        action, eligible, reason = decide(measurement)
        decision = dict(id=key, assembly=a.assembly, scaffold=row['scaffold'], action=action,
                        chromosome_scale_status='UNRESOLVED', auto_eligible='yes' if eligible else 'no',
                        reason=reason, policy_version='gap-v1', evidence_packet_id=a.assembly+'.sequence_context')
        decisions.append(decision)
        if eligible:
            for field in ('gap_start', 'gap_end', 'cut_bp', 'assessment_sha256'):
                if field not in measurement:
                    raise ValueError('Eligible action lacks source-bound coordinate '+field)
            if measurement['assessment_sha256'] != row.get('assembly_sha256'):
                raise ValueError('Measurement checksum differs from candidate')
            actions.append(dict(id=key, assembly=a.assembly, scaffold=row['scaffold'], action=action,
                                coordinate_stage='pre_finishing', assessment_sha256=row['assembly_sha256'],
                                cut_bp=measurement['cut_bp'], gap_start=measurement['gap_start'], gap_end=measurement['gap_end'],
                                decision_source='automatic', evidence_packet_id=decision['evidence_packet_id'],
                                reason=reason, auto_eligible='yes', policy_version='gap-v1'))
    fields = ['id', 'assembly', 'scaffold', 'action', 'chromosome_scale_status', 'auto_eligible', 'reason', 'policy_version', 'evidence_packet_id']
    write(out/'decisions.tsv', decisions, fields)
    write(out/'actions.tsv', actions, ACTION_FIELDS)
    (out/'measurement_audit.json').write_text(json.dumps(dict(policy='gap-v1', measurements=measurements,
        measurement_status='available' if evidence_path.exists() else 'not_available',
        candidate_file_sha256=hashlib.sha256(Path(a.calls).read_bytes()).hexdigest()), indent=2)+'\n')
    table = ''.join('<tr>'+''.join('<td>'+html.escape(str(row[k]))+'</td>' for k in fields)+'</tr>' for row in decisions)
    (out/'report.html').write_text('<!doctype html><meta charset="utf-8"><title>Chimera decisions</title>'
        '<h1>'+html.escape(a.assembly)+'</h1><p>Coordinates and evidence refer to the pre-finishing assembly. '
        'Unresolved candidates remain uncut. Local continuity does not establish biological fusion.</p>'
        '<table border="1"><tr>'+''.join('<th>'+k+'</th>' for k in fields)+'</tr>'+table+'</table>')


if __name__ == '__main__':
    main()
