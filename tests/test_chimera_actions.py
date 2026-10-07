import sys
import unittest
import tempfile
import subprocess
import csv
import hashlib
import json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'py_scripts'))
from chimera_actions import validate_actions, split_sequences


class ActionTests(unittest.TestCase):
    def row(self, **changes):
        row = dict(id='discovered-boundary', assessment_sha256='hash', coordinate_stage='pre_finishing',
                   scaffold='composite', action='UNJOIN_UNSUPPORTED', cut_bp='6', gap_start='4', gap_end='6',
                   decision_source='review', evidence_packet_id='packet', reason='reviewed gap')
        row.update(changes)
        return row

    def test_gap_edge_preserves_every_base_and_neutral_identity(self):
        seq = {'composite': 'AAAANNCCCC', 'other': 'GGGG'}
        cuts = validate_actions([self.row()], 'hash', seq, 4, 'file')
        output, lift = split_sequences(seq, cuts)
        self.assertEqual(output['composite_sub_0_6'], 'AAAANN')
        self.assertEqual(output['composite_sub_6_10'], 'CCCC')
        self.assertEqual(''.join(output[e['piece']] for e in lift if e['parent']=='composite'), seq['composite'])

    def test_reject_wrong_hash_non_gap_fraction_duplicate_and_small_piece(self):
        for change in [dict(assessment_sha256='wrong'), dict(gap_start='3'), dict(cut_bp='6.0')]:
            with self.assertRaises(ValueError):
                validate_actions([self.row(**change)], 'hash', {'composite':'AAAANNCCCC'}, 4, 'file')
        with self.assertRaises(ValueError):
            validate_actions([self.row(), self.row(id='second')], 'hash', {'composite':'AAAANNCCCC'}, 4, 'file')
        with self.assertRaises(ValueError):
            validate_actions([self.row()], 'hash', {'composite':'AAAANNCCCC'}, 5, 'file')

    def test_auto_requires_policy_and_internal_review_cannot_authorize_auto(self):
        with self.assertRaises(ValueError):
            validate_actions([self.row()], 'hash', {'composite':'AAAANNCCCC'}, 4, 'auto')
        with self.assertRaises(ValueError):
            validate_actions([self.row(action='BREAK_PROBABLE_MISJOIN', reviewer='human', localization_status='localized',
                                       auto_eligible='yes', policy_version='gap-v2')], 'hash', {'composite':'AAAANNCCCC'}, 4, 'auto')

    def test_collision(self):
        with self.assertRaises(ValueError):
            split_sequences({'composite':'AAAANNCCCC', 'composite_sub_0_6':'GG'}, {'composite':[6]})

    def test_pipeline_entrypoint_round_trip(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root/'source.fa'
            source.write_text('>composite\nAAAANNCCCC\n>other\nGGGG\n')
            name_map = root/'names.tsv'
            name_map.write_text('old_name\tnew_name\torient\torder\tlength\tclass\tref_span\tflags\tchromosome_member\n'
                                'composite\tchr1+chr2+chr3\t-\t1\t10\tchromosome\told\t-\tyes\n'
                                'other\tchr4\t+\t2\t4\tchromosome\told\t-\tyes\n')
            row = self.row(assembly='assembly', assessment_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
            actions = root/'actions.tsv'
            with actions.open('w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=list(row), delimiter='\t')
                writer.writeheader(); writer.writerow(row)
            script = Path(__file__).resolve().parents[1]/'py_scripts'/'break_chimeras.py'
            output = root/'output.fa'
            subprocess.run([sys.executable, str(script), '--fasta', str(source), '--name-map', str(name_map),
                            '--actions', str(actions), '--assembly', 'assembly', '--out-fasta', str(output),
                            '--out-name-map', str(root/'output.tsv'), '--audit', str(root/'audit.tsv'),
                            '--mode', 'file', '--min-piece-bp', '4'], check=True)
            verification = json.loads(Path(str(output)+'.verification.json').read_text())
            self.assertTrue(verification['exact_parent_reconstruction'])
            self.assertEqual(verification['total_bp'], 14)
            with (root/'output.tsv').open() as handle:
                rewritten = list(csv.DictReader(handle, delimiter='\t'))
            self.assertEqual(rewritten[0]['new_name'], 'composite_sub_0_6')
            self.assertEqual(rewritten[1]['class'], 'unplaced')
            self.assertEqual(rewritten[1]['orient'], '+')
            self.assertEqual(rewritten[2]['new_name'], 'chr4')
            # Empty auto selections must successfully pass through every record.
            actions.write_text('\t'.join(row)+'\n')
            subprocess.run([sys.executable, str(script), '--fasta', str(source), '--name-map', str(name_map),
                            '--actions', str(actions), '--assembly', 'assembly', '--out-fasta', str(output),
                            '--out-name-map', str(root/'output.tsv'), '--audit', str(root/'audit.tsv'),
                            '--mode', 'auto', '--min-piece-bp', '4'], check=True)
            verification = json.loads(Path(str(output)+'.verification.json').read_text())
            self.assertEqual(verification['output_records'], 2)


if __name__ == '__main__':
    unittest.main()

