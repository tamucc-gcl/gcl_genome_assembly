"""Synthetic integration fixtures; positive eligibility is a software control, not biological evidence."""
import argparse
import csv
import hashlib
import json
from pathlib import Path


def table(path, rows):
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--out', required=True)
    p.add_argument('--scenario', choices=['unresolved', 'eligible', 'manual'], required=True)
    a = p.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    sequence = '>composite\nAAAANNCCCC\n>other\nGGGG\n'
    (out/'assessment.fa').write_text(sequence)
    digest = hashlib.sha256((out/'assessment.fa').read_bytes()).hexdigest()
    table(out/'names.tsv', [dict(old_name='composite', new_name='chr1+chr2+chr3', orient='-', order=1,
        length=10, **{'class':'chromosome'}, ref_span='.', flags='.', chromosome_member='yes'),
        dict(old_name='other', new_name='chr4', orient='+', order=2, length=4,
             **{'class':'chromosome'}, ref_span='.', flags='.', chromosome_member='yes')])
    row = dict(assembly='synthetic', assembly_sha256=digest, coordinate_stage='pre_finishing',
               scaffold='composite', transition_lo='4', transition_hi='6', candidate_verdict='REVIEW')
    table(out/'calls.tsv', [row])
    key = hashlib.sha256(json.dumps({k:str(v) for k,v in row.items()}, sort_keys=True).encode()).hexdigest()[:20]
    context = out/'context'; context.mkdir(exist_ok=True)
    measurement = dict(assessment_sha256=digest, gap_start=4, gap_end=6, cut_bp=6,
        assessment_scaffold_length=10,verified_gap=True, independent_discordant_individuals=3, hifi_informative=True,
        hifi_spanning_molecules=0, hic_informative=True, hic_support_loss=True,
        matched_controls_pass=True, alternative_placements_checked=True, graph_contradiction=False)
    (context/'decision_measurements.json').write_text(json.dumps({key:measurement} if a.scenario=='eligible' else {}))
    action = dict(selected='NO' if a.scenario=='unresolved' else 'YES',reviewer='test-reviewer',localization_status='proposed_gap',id=key, assembly='synthetic', coordinate_stage='pre_finishing', assessment_sha256=digest,
        scaffold='composite', action='UNJOIN_UNSUPPORTED', cut_bp=6, gap_start=4, gap_end=6,
        decision_source='review', evidence_packet_id='synthetic-fixture', reason='synthetic gap-edge control')
    table(out/'manual-actions.tsv', [action])


if __name__ == '__main__':
    main()
