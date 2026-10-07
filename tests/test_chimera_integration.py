"""Exercise fixture -> decision -> action -> verification with real command-line tools."""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class IntegrationTests(unittest.TestCase):
    def test_all_three_routes(self):
        for scenario in ('unresolved', 'eligible', 'manual'):
            with self.subTest(scenario=scenario), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                def run(script, *args):
                    result = subprocess.run([sys.executable, str(ROOT/script), *map(str,args)], capture_output=True, text=True)
                    self.assertEqual(result.returncode, 0, result.stderr)
                run('tests/prepare_chimera_integration.py', '--out', root, '--scenario', scenario)
                run('py_scripts/chimera_adjudicate.py', '--calls', root/'calls.tsv', '--assembly', 'synthetic',
                    '--context', root/'context', '--out', root/'decisions', '--min-piece-bp',4)
                actions = root/'manual-actions.tsv'
                if scenario=='unresolved':
                    actions.write_text(actions.read_text().splitlines()[0]+'\n')
                run('py_scripts/break_chimeras.py', '--fasta', root/'assessment.fa', '--name-map', root/'names.tsv',
                    '--actions', actions, '--assembly', 'synthetic', '--out-fasta', root/'output.fa',
                    '--out-name-map', root/'output.tsv', '--audit', root/'audit.tsv',
                    '--mode', 'file', '--min-piece-bp', 4)
                verify = json.loads((root/'output.fa.verification.json').read_text())
                self.assertEqual(verify['output_records'], 2 if scenario=='unresolved' else 3)
                self.assertEqual(verify['total_bp'], 14)
                self.assertTrue(verify['exact_parent_reconstruction'])
                self.assertTrue((root/'decisions/report.html').exists())
                self.assertTrue((root/'decisions/measurement_audit.json').exists())


if __name__ == '__main__':
    unittest.main()
