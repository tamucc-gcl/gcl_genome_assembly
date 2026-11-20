#!/usr/bin/env python3
"""Run a bounded, synthetic acceptance test through the real post-FCS cleanup workflow."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tarfile
import time


FIELDS = ['action_id', 'assembly', 'stage', 'source_sha256', 'sequence', 'action',
          'start0', 'end0', 'selected', 'reviewer', 'reason']


def checksum(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fasta_records(path):
    records, name = {}, None
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            name = line[1:].split()[0]
            if name in records:
                raise AssertionError('Duplicate output identifier: ' + name)
            records[name] = ''
        elif line.strip():
            if name is None:
                raise AssertionError('Sequence without a FASTA header')
            records[name] += line.strip()
    return records


def table(path):
    with Path(path).open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2) + '\n')


def write_sheet(path, rows):
    with Path(path).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def package_review(out):
    archive = Path(str(out) + '.review.tar.gz')
    with tarfile.open(archive, 'w:gz') as bundle:
        for path in sorted(out.iterdir()):
            if path.name != 'work' and not path.name.startswith('.nextflow'):
                bundle.add(path, arcname=out.name + '/' + path.name)
    return archive


def create_fixture(directory):
    directory.mkdir()
    originals = {
        'whole_hap1': {'host': 'ACGTacgtNN', 'foreign': 'TTTTCCCC'},
        'internal_hap1': {'joined': 'AACCGGTTNN', 'spare': 'GATTACA'},
        'trim_hap1': {'ends': 'AACCGGTTNN', 'spare': 'TGCATGCA'},
        'unselected_hap1': {'same': 'aCGTNNtt'},
        'untouched_hap2': {'same': 'ccGGNNaA'},
    }
    cases = []
    for identity, records in originals.items():
        source = directory / (identity + '.fasta')
        # Deliberate descriptions, mixed case, and wrapping check byte-preserving bypass.
        source.write_text(''.join(f'>{name} fixture description\n{seq[:3]}\n{seq[3:]}\n'
                                  for name, seq in records.items()))
        cases.append(dict(id=identity, source=str(source), original=records,
                          source_sha256=checksum(source),
                          expected={mode: dict(records) for mode in ('baseline', 'changed', 'noop')}))
    by_id = {row['id']: row for row in cases}
    by_id['whole_hap1']['expected']['baseline'] = {'host': 'ACGTacgtNN'}
    by_id['internal_hap1']['expected']['baseline'] = {
        'joined__retained_0_2': 'AA', 'joined__retained_6_10': 'TTNN', 'spare': 'GATTACA'}
    by_id['trim_hap1']['expected']['baseline'] = {'ends__retained_2_8': 'CCGGTT', 'spare': 'TGCATGCA'}
    for row in cases:
        row['expected']['changed'] = dict(row['expected']['baseline'])
    by_id['internal_hap1']['expected']['changed'] = {
        'joined__retained_0_3': 'AAC', 'joined__retained_6_10': 'TTNN', 'spare': 'GATTACA'}
    decisions = [
        ('whole_hap1', 'foreign', 'REMOVE_SEQUENCE', 0, 8, 'YES'),
        ('internal_hap1', 'joined', 'REMOVE_INTERVAL', 2, 6, 'YES'),
        ('trim_hap1', 'ends', 'REMOVE_INTERVAL', 0, 2, 'YES'),
        ('trim_hap1', 'ends', 'REMOVE_INTERVAL', 8, 10, 'YES'),
        ('unselected_hap1', 'same', 'REMOVE_INTERVAL', 2, 4, 'NO'),
    ]
    rows = [dict(zip(FIELDS, [f'fixture_{i}', identity, 'post_fcs',
                  by_id[identity]['source_sha256'], name, action, lo, hi, selected,
                  'acceptance fixture', 'Synthetic test only']))
            for i, (identity, name, action, lo, hi, selected) in enumerate(decisions)]
    write_sheet(directory / 'baseline.tsv', rows)
    changed = [dict(row) for row in rows]
    changed[1]['start0'] = 3
    write_sheet(directory / 'changed.tsv', changed)
    write_sheet(directory / 'noop.tsv', [dict(row, selected='NO') for row in rows])
    for label, field, value in [('wrong_checksum', 'source_sha256', '0' * 64),
                                ('wrong_stage', 'stage', 'final'),
                                ('unknown_assembly', 'assembly', 'not_an_input_hap1')]:
        invalid = [dict(rows[1], **{field: value})]
        write_sheet(directory / (label + '.tsv'), invalid)
    write_sheet(directory / 'overlap.tsv', [rows[1], dict(rows[1], action_id='overlap', start0=4, end0=8)])
    write_json(directory / 'manifest.json', dict(cases=cases))
    return cases


def verify_one(args):
    expected = json.loads(args.expected.read_text())
    actual = fasta_records(args.fasta)
    if actual != expected['retained']:
        raise AssertionError(f'{expected["assembly"]}: downstream received incorrect sequences: {actual}')
    source = Path(expected['source'])
    if expected['changed']:
        if args.fasta.name != 'retained.fasta' or args.fasta.resolve() == source.resolve():
            raise AssertionError('Reviewed assembly bypassed the retained FASTA')
    elif args.fasta.resolve() != source.resolve() or args.fasta.read_bytes() != source.read_bytes():
        raise AssertionError('Unselected assembly did not preserve its original input path and bytes')
    write_json(args.out, dict(assembly=expected['assembly'], input=str(args.fasta.resolve()),
                             sha256=checksum(args.fasta), records=len(actual),
                             retained_bp=sum(map(len, actual.values())), changed=expected['changed']))


def verify_phase(output, fixture, cases, scenario):
    findings = []
    rows = table(fixture / ('changed.tsv' if scenario == 'changed' else scenario + '.tsv'))
    for case in cases:
        identity = case['id']
        if checksum(case['source']) != case['source_sha256']:
            raise AssertionError('The source FASTA was modified')
        report = json.loads((output / 'downstream' / (identity + '.json')).read_text())
        expected = case['expected'][scenario]
        if report['retained_bp'] != sum(map(len, expected.values())):
            raise AssertionError(identity + ': downstream base count')
        audit = output / 'assembly/contamination' / identity / 'cleanup/audit'
        selected = [row for row in rows if row['assembly'] == identity and row['selected'] == 'YES']
        if not selected:
            if audit.exists():
                raise AssertionError(identity + ': unselected assembly entered cleanup')
            findings.append(dict(assembly=identity, changed=False, original_path_and_bytes_preserved=True))
            continue
        retained, removed = fasta_records(audit / 'retained.fasta'), fasta_records(audit / 'removed.fasta')
        if retained != expected:
            raise AssertionError(identity + ': published retained FASTA differs from downstream input')
        pieces = {**retained, **removed}
        coordinates = table(audit / 'coordinate_map.tsv')
        if {row['sequence'] for row in coordinates} != set(case['original']):
            raise AssertionError(identity + ': missing or unknown source in coordinate audit')
        if {row['output_sequence'] for row in coordinates} != set(pieces):
            raise AssertionError(identity + ': missing or unknown output in coordinate audit')
        for name, source_seq in case['original'].items():
            lifted = [row for row in coordinates if row['sequence'] == name]
            cursor, reconstructed = 0, ''
            for row in lifted:
                lo, hi = int(row['source_start0']), int(row['source_end0'])
                piece = pieces[row['output_sequence']]
                outputs = retained if row['disposition'] == 'retained' else removed
                if row['disposition'] not in ('retained', 'removed') or row['output_sequence'] not in outputs:
                    raise AssertionError(identity + ': incorrect coordinate disposition')
                if lo != cursor or hi <= lo or piece != source_seq[lo:hi]:
                    raise AssertionError(identity + ': invalid coordinate partition')
                if (int(row['output_start0']), int(row['output_end0'])) != (0, len(piece)):
                    raise AssertionError(identity + ': invalid output coordinates')
                if hashlib.sha256(piece.encode('ascii')).hexdigest() != row['sha256']:
                    raise AssertionError(identity + ': invalid coordinate checksum')
                reconstructed += piece
                cursor = hi
            if reconstructed != source_seq:
                raise AssertionError(identity + ': audit cannot reconstruct original sequence')
        verification = json.loads((audit / 'verification.json').read_text())
        if verification['source_bp'] != verification['retained_bp'] + verification['removed_bp']:
            raise AssertionError(identity + ': base accounting')
        if (verification['source_bp'], verification['retained_bp'], verification['removed_bp']) != (
                sum(map(len, case['original'].values())), sum(map(len, retained.values())),
                sum(map(len, removed.values()))):
            raise AssertionError(identity + ': audit totals differ from actual sequences')
        if verification['source_sha256'] != case['source_sha256'] or verification['selected_actions'] != len(selected):
            raise AssertionError(identity + ': incorrect action/source audit')
        if verification['retained_sha256'] != checksum(audit / 'retained.fasta') or verification['removed_sha256'] != checksum(audit / 'removed.fasta'):
            raise AssertionError(identity + ': output checksums')
        applied = table(audit / 'applied_actions.tsv')
        if {tuple(row[field] for field in FIELDS) for row in applied} != {
                tuple(row[field] for field in FIELDS) for row in selected}:
            raise AssertionError(identity + ': applied actions')
        findings.append(dict(assembly=identity, changed=True, **{
            key: verification[key] for key in ('retained_bp', 'removed_bp', 'selected_actions')}))
    return findings


def run_acceptance(args):
    repo = Path(__file__).resolve().parents[2]
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    fixture = out / 'fixture'
    cases = create_fixture(fixture)
    config = out / 'acceptance.config'
    config.write_text("""process {
    withName: '.*CONTAMINATION_DECISIONS|.*APPLY_CONTAMINATION_REVIEW|CHECK_CLEANUP_INPUT' {
        cpus = 1
        memory = '1 GB'
        time = '10m'
        errorStrategy = 'terminate'
        scratch = false
    }
}
dag.enabled = false
report.enabled = false
timeline.enabled = false
""" + ('conda.enabled = false\nsingularity.enabled = false\n' if args.no_conda else ''))
    summaries, session = [], None
    phases = [('baseline', 'baseline'), ('repeat', 'baseline'), ('changed', 'changed'),
              ('noop', 'noop'), ('disabled', 'noop')]
    if not args.skip_negative:
        phases += [(name, 'baseline') for name in ('wrong_checksum', 'wrong_stage', 'unknown_assembly', 'overlap')]
    for label, scenario in phases:
        phase = out / label
        phase.mkdir()
        sheet = {'repeat': 'baseline', 'disabled': 'noop'}.get(label, label)
        command = [args.nextflow, '-log', str(phase / 'nextflow.log'), 'run', str(repo / 'tests/contamination_cleanup.nf'),
                   '-profile', args.profile, '-c', str(repo / 'nextflow.config'),
                   '-c', str(config), '-w', str(out / 'work'), '-ansi-log', 'false',
                   '-with-trace', str(phase / 'trace.tsv'),
                   '--cleanup_fixture', str(fixture / 'manifest.json'), '--cleanup_scenario', scenario,
                   '--contamination_review', 'false' if label == 'disabled' else str(fixture / (sheet + '.tsv')),
                   '--run_contamination_evidence', 'false', '--outdir', str(phase)]
        if session:
            command += ['-resume', session]
        write_json(phase / 'command.json', command)
        print(f'Running cleanup acceptance phase: {label}', flush=True)
        started = time.monotonic()
        with (phase / 'console.log').open('w') as log:
            result = subprocess.run(command, cwd=out, stdout=log, stderr=subprocess.STDOUT)
        negative = label not in ('baseline', 'repeat', 'changed', 'noop', 'disabled')
        if negative:
            expected_error = {'wrong_checksum': 'does not match source FASTA/stage',
                              'wrong_stage': 'does not match source FASTA/stage',
                              'unknown_assembly': 'reference unknown assemblies',
                              'overlap': 'Overlapping selected removals'}[label]
            console = (phase / 'console.log').read_text()
            if result.returncode == 0 or expected_error not in console:
                raise RuntimeError(f'{label}: expected specific rejection; see {phase / "console.log"}')
            if (phase / 'assembly/contamination/internal_hap1/cleanup/audit').exists():
                raise AssertionError(label + ': invalid decisions published a cleanup audit')
            summaries.append(dict(phase=label, expected_rejection=True))
            continue
        if result.returncode:
            raise RuntimeError(f'{label} failed; see {phase / "console.log"}')
        findings = verify_phase(phase, fixture, cases, scenario)
        trace = table(phase / 'trace.tsv')
        statuses = Counter(row['status'] for row in trace)
        if any(row['name'].split(' (')[0].split(':')[-1] not in
               {'CONTAMINATION_DECISIONS', 'APPLY_CONTAMINATION_REVIEW', 'CHECK_CLEANUP_INPUT'} for row in trace):
            raise AssertionError('The cleanup fixture unexpectedly launched evidence/assembly work')
        if label == 'repeat' and (len(trace) != 9 or statuses != {'CACHED': 9}):
            raise AssertionError('Unchanged repeat did not reuse every task: ' + str(statuses))
        if label == 'changed':
            cached_bypass = {row['name'].split(' (')[-1].rstrip(')') for row in trace
                             if row['name'].startswith('CHECK_CLEANUP_INPUT') and row['status'] == 'CACHED'}
            if cached_bypass != {'unselected_hap1', 'untouched_hap2'}:
                raise AssertionError('Unchanged downstream assemblies did not retain their caches')
        if label == 'disabled' and (len(trace) != 5 or statuses != {'CACHED': 5}):
            raise AssertionError('Disabling cleanup did not preserve the original downstream inputs')
        if session is None:
            session = re.search(r'Session UUID: ([a-f0-9-]+)', (phase / 'nextflow.log').read_text()).group(1)
        summaries.append(dict(phase=label, statuses=dict(statuses), assemblies=findings,
                              seconds=round(time.monotonic() - started, 2)))
        write_json(phase / 'acceptance.json', summaries[-1])
    write_json(out / 'acceptance_summary.json', dict(status='PASS', phases=summaries))
    archive = package_review(out)
    print('PASS: cleanup routing, exact sequences, coordinate audits, cache reuse, and rejection checks', flush=True)
    print(out / 'acceptance_summary.json', flush=True)
    print('Review archive: ' + str(archive), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    run = commands.add_parser('run')
    run.add_argument('--out', type=Path, required=True, help='New output directory, outside production outputs')
    run.add_argument('--nextflow', default='nextflow', help='Nextflow executable')
    run.add_argument('--profile', choices=['local', 'slurm'], default='local')
    run.add_argument('--no-conda', action='store_true', help='Use an existing Python environment for local tests')
    run.add_argument('--skip-negative', action='store_true', help='Skip expected-failure workflow runs')
    verify = commands.add_parser('verify-one')
    verify.add_argument('--fasta', type=Path, required=True)
    verify.add_argument('--expected', type=Path, required=True)
    verify.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    (run_acceptance if args.command == 'run' else verify_one)(args)


if __name__ == '__main__':
    main()
