"""Validate source-bound actions and generate sequence-preserving neutral pieces."""
import hashlib
import re


def integer(value):
    if not re.fullmatch(r'\d+', str(value)):
        raise ValueError('Expected an integer coordinate: ' + str(value))
    return int(value)


def validate_actions(rows, checksum, sequences, minimum, mode):
    selected, identities = {}, set()
    for row in rows:
        if row.get('assessment_sha256') != checksum or row.get('coordinate_stage') != 'pre_finishing':
            raise ValueError('Action does not match this pre-finishing FASTA')
        scaffold = row.get('scaffold')
        if scaffold not in sequences:
            raise ValueError('Unknown action scaffold: ' + str(scaffold))
        cut = integer(row.get('cut_bp'))
        sequence = sequences[scaffold]
        if not 0 < cut < len(sequence):
            raise ValueError('Cut outside scaffold: ' + scaffold)
        if not all(row.get(k) for k in ('id', 'decision_source', 'evidence_packet_id', 'reason')):
            raise ValueError('Action requires identity, provenance and rationale')
        if row['id'] in identities:
            raise ValueError('Duplicate action identity')
        identities.add(row['id'])
        if row['action'] == 'UNJOIN_UNSUPPORTED':
            lo, hi = integer(row.get('gap_start')), integer(row.get('gap_end'))
            if not (0 < lo < hi < len(sequence) and lo <= cut <= hi and set(sequence[lo:hi].upper()) == {'N'}):
                raise ValueError('Unjoin requires a literal all-N gap and a cut within its edges')
        elif row['action'] == 'BREAK_PROBABLE_MISJOIN':
            if mode == 'auto' or not row.get('reviewer') or row.get('localization_status') != 'localized':
                raise ValueError('Internal cutting requires explicit review and localization')
        else:
            raise ValueError('Unsupported selected action')
        if mode == 'auto' and (row.get('auto_eligible') != 'yes' or row.get('policy_version') != 'gap-v1'):
            raise ValueError('Automatic action lacks decision-stage eligibility')
        if cut in selected.setdefault(scaffold, []):
            raise ValueError('Duplicate/conflicting cut')
        selected[scaffold].append(cut)
    for scaffold, cuts in selected.items():
        cuts.sort()
        bounds = [0] + cuts + [len(sequences[scaffold])]
        if any(hi-lo < minimum for lo, hi in zip(bounds, bounds[1:])):
            raise ValueError('Cut set creates a piece below minimum length')
    return selected


def split_sequences(sequences, selected):
    output, lift = {}, []
    for parent, sequence in sequences.items():
        bounds = [0] + selected.get(parent, []) + [len(sequence)]
        reconstruction = []
        for lo, hi in zip(bounds, bounds[1:]):
            name = '%s_sub_%d_%d' % (parent, lo, hi) if parent in selected else parent
            if name in output or (name != parent and name in sequences):
                raise ValueError('Piece name collision')
            fragment = sequence[lo:hi]
            output[name] = fragment
            reconstruction.append(fragment)
            lift.append(dict(parent=parent, piece=name, parent_start=lo, parent_end=hi,
                             piece_start=0, piece_end=hi-lo,
                             sequence_sha256=hashlib.sha256(fragment.encode()).hexdigest()))
        if ''.join(reconstruction) != sequence:
            raise ValueError('Parent reconstruction failed')
    return output, lift
