"""Plan conservative reference-versus-assembly SyRI comparisons from eligibility."""
import argparse
import csv
import json
from pathlib import Path
import re


def chromosomes(row):
    names = row.get('reference_contigs', [])
    # Only already harmonized, single-component chromosome assignments.
    return sorted(names) if names and all(re.fullmatch(r'chr[0-9]+_1', n) for n in names) and len(set(names)) == len(names) else []


def chromosome_issue(row):
    names = row.get('reference_contigs', [])
    if any('+' in n for n in names):
        return 'composite_chromosome_assignment_not_supported'
    if names and all(re.fullmatch(r'chr[0-9]+_[0-9]+', n) for n in names):
        if any(not n.endswith('_1') for n in names):
            return 'split_chromosome_representation_not_supported'
    return 'missing_or_unrecognized_chromosome_names'


def plan(data, selected=None):
    assemblies = {a['id']: a for a in data['assemblies']}
    if len(assemblies) != len(data['assemblies']):
        raise ValueError('Duplicate assembly IDs')
    rows, pairs = [], []
    for cohort in data['cohorts']:
        taxid = str(cohort['taxid'])
        members = sorted((a for a in assemblies.values() if str(a['taxid']) == taxid), key=lambda a: a['id'])
        ref = assemblies.get(cohort.get('reference_id'))
        ref_chroms = chromosomes(ref) if ref else []
        for query in members:
            reason = ''
            if cohort.get('missing'):
                reason = 'incomplete_species_assembly_set'
            elif ref is None:
                reason = 'no_shared_species_reference'
            elif str(ref['taxid']) != taxid:
                raise ValueError('Reference belongs to another species')
            elif query['id'] == ref['id']:
                reason = 'reference_self_comparison'
            elif query.get('completion') != 'finalized':
                reason = 'assembly_not_finalized'
            elif not query.get('representation_eligible'):
                reason = 'unsupported_assembly_representation'
            elif not query.get('chromosome_scale') or not ref.get('chromosome_scale'):
                reason = 'not_chromosome_scale'
            elif not ref_chroms:
                reason = 'reference_' + chromosome_issue(ref)
            elif not chromosomes(query):
                reason = chromosome_issue(query)
            elif chromosomes(query) != ref_chroms:
                reason = 'chromosome_sets_differ'
            row = dict(taxid=taxid, reference=ref['id'] if ref else '', query=query['id'],
                       status='skipped' if reason else 'planned', reason=reason or 'matched_harmonized_chromosome_sets')
            rows.append(row)
            if not reason:
                pairs.append(dict(row, key='pair_' + query['id'].encode().hex(),
                    reference_fasta=ref['fasta'], query_fasta=query['fasta'], chromosomes=ref_chroms))
    if selected:
        selected = set(selected)
        unavailable = selected - {p['query'] for p in pairs}
        if unavailable:
            raise ValueError('Requested queries are unknown or ineligible: ' + ', '.join(sorted(unavailable)))
        pairs = [p for p in pairs if p['query'] in selected]
        for row in rows:
            if row['status'] == 'planned' and row['query'] not in selected:
                row.update(status='not_selected', reason='outside_requested_query_set')
    return dict(schema_version=1, pairs=pairs, outcomes=rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('manifest', type=Path)
    p.add_argument('--selection', type=Path)
    a = p.parse_args()
    selected = json.loads(a.selection.read_text()) if a.selection else []
    result = plan(json.loads(a.manifest.read_text()), selected)
    Path('sv_plan.json').write_text(json.dumps(result, indent=2) + '\n')
    with open('sv_plan.tsv', 'w', newline='') as handle:
        w = csv.DictWriter(handle, fieldnames=['taxid', 'reference', 'query', 'status', 'reason'], delimiter='\t')
        w.writeheader()
        w.writerows(result['outcomes'])
    lines = ['# Assembly structural-variation plan', '',
        'Optional SyRI branch; planning is not execution or biological validation.', '',
        'Reference selection comes from the shared assembly eligibility manifest, independent of graph execution. '
        'This initial contract requires matching, unambiguous harmonized chromosome sets. '
        'Unplaced sequence is excluded. Missing shared references and unsuitable pairs are explicit skips.', '',
        '| Taxid | Reference | Query | Status | Reason |', '|---|---|---|---|---|']
    lines += ['| ' + ' | '.join(str(r[k]) for k in ('taxid', 'reference', 'query', 'status', 'reason')) + ' |' for r in result['outcomes']]
    Path('sv_plan.md').write_text('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
