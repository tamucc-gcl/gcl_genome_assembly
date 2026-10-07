#!/usr/bin/env python3
"""Apply normalized action TSVs; retain source coordinates and reconstruction audit."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from chimera_actions import validate_actions, split_sequences


def rows(path):
    with open(path) as handle:
        reader = csv.DictReader((x for x in handle if x.strip() and not x.lstrip().startswith('#')), delimiter='\t')
        if reader.fieldnames and len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError('Duplicate TSV columns')
        result = list(reader)
        if any(None in r or None in r.values() for r in result):
            raise ValueError('Malformed TSV row')
        return result


def fasta(path):
    result, name = {}, None
    with open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                name = line[1:].split()[0]
                if name in result:
                    raise ValueError('Duplicate FASTA name')
                result[name] = []
            elif line.strip():
                if name is None:
                    raise ValueError('Sequence before FASTA header')
                result[name].append(line.strip())
    return {k: ''.join(v) for k, v in result.items()}


def write(path, records, fields):
    with open(path, 'w') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(records)


def digest(path):
    with open(path, 'rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('fasta', 'name-map', 'actions', 'assembly', 'out-fasta', 'out-name-map'):
        p.add_argument('--'+key, required=True)
    p.add_argument('--audit', required=True)
    p.add_argument('--mode', choices=('auto', 'file'), default='file')
    p.add_argument('--min-piece-bp', type=int, default=1000000)
    a = p.parse_args()
    if a.mode == 'auto':
        raise ValueError('Automated cutting is deferred; supply a reviewed file with selected=YES, reviewer and reason')
    if a.min_piece_bp < 1:
        raise ValueError('Minimum piece length must be positive')
    sequences, maps = fasta(a.fasta), rows(a.name_map)
    by_name = {r['old_name']: r for r in maps}
    if not maps or len(by_name) != len(maps) or set(by_name) != set(sequences):
        raise ValueError('Name map and FASTA must contain the same unique records')
    required = {'id', 'assembly', 'coordinate_stage', 'assessment_sha256', 'scaffold', 'action',
                'cut_bp', 'gap_start', 'gap_end', 'decision_source', 'evidence_packet_id', 'reason'}
    with open(a.actions) as handle:
        header = next((line.rstrip().split('\t') for line in handle if line.strip() and not line.lstrip().startswith('#')), [])
    if not required.issubset(header):
        raise ValueError('Action table missing required columns: '+', '.join(sorted(required-set(header))))
    action_rows = rows(a.actions)
    if any(not r.get('assembly') for r in action_rows):
        raise ValueError('Every action must identify its assembly')
    if a.mode == 'file':
        if not {'selected','reviewer','localization_status'}.issubset(header):
            raise ValueError('Review file requires selected, reviewer and localization_status columns')
        if any(r['selected'] not in ('YES','NO') for r in action_rows):
            raise ValueError('Every review row requires selected=YES or NO')
        action_rows=[r for r in action_rows if r['selected']=='YES']
        if any(not r.get('reviewer','').strip() or not r.get('reason','').strip() or r.get('decision_source')!='review' for r in action_rows):
            raise ValueError('Every selected cut requires reviewer, reason and decision_source=review')
    selected_rows = [r for r in action_rows if r['assembly'] == a.assembly]
    selected = validate_actions(selected_rows, digest(a.fasta), sequences, a.min_piece_bp, a.mode)
    output, lift = split_sequences(sequences, selected)
    new_maps = []
    for entry in lift:
        row = dict(by_name[entry['parent']])
        if entry['parent'] in selected:
            row.update(old_name=entry['piece'], new_name=entry['piece'], length=str(entry['piece_end']), orient='+')
            if 'class' in row:
                row['class'] = 'unplaced'
            if 'chromosome_member' in row:
                row['chromosome_member'] = 'no'
            if 'ref_span' in row:
                row['ref_span'] = '.'
            if 'flags' in row:
                row['flags'] = (row.get('flags', '')+';chimera_broken;chromosome_assignment_pending').strip(';')
        new_maps.append(row)
    with open(a.out_fasta, 'w') as handle:
        for name, sequence in output.items():
            handle.write('>'+name+'\n')
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start:start+60]+'\n')
    # Check the actual written artifact, not only the in-memory slices.
    actual = fasta(a.out_fasta)
    if actual != output:
        raise ValueError('Written FASTA verification failed')
    write(a.out_name_map, new_maps, list(maps[0]))
    write(a.out_fasta+'.coordinate_lift.tsv', lift, list(lift[0]))
    verification = dict(assembly=a.assembly, input_sha256=digest(a.fasta), output_sha256=digest(a.out_fasta),
                        input_records=len(sequences), output_records=len(output), total_bp=sum(map(len, output.values())),
                        exact_parent_reconstruction=True, coordinate_stage='pre_finishing', selected_actions=selected_rows)
    Path(a.out_fasta+'.verification.json').write_text(json.dumps(verification, indent=2)+'\n')
    write(a.audit, [dict(metric=k, value=v) for k, v in [('assembly', a.assembly), ('mode', a.mode),
          ('scaffolds_broken', len(selected)), ('cuts_total', sum(map(len, selected.values()))),
          ('exact_parent_reconstruction', 'yes')]], ['metric', 'value'])


if __name__ == '__main__':
    main()
