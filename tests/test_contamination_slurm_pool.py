from argparse import Namespace
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'py_scripts'))
import contamination_slurm_pool as pool


class PoolTests(unittest.TestCase):
    def args(self, root):
        plan = root / 'plan.json'
        plan.write_text(json.dumps(dict(tasks=[{'assembly': 'A', 'name': 'q'}], threads=8)))
        return Namespace(plan=plan, out=root / 'pool', helper=root / 'contamination_workers.py', nodes=2,
                         slots=8, scratch='/scratch', partition='normal', nodelist='node1,node2',
                         walltime='48:00:00', timeout=100, poll_seconds=0.001)

    def test_complete_resume_requests_no_nodes(self):
        with tempfile.TemporaryDirectory() as temp, patch.object(pool, 'Results') as results, patch.object(pool, 'command') as command:
            args = self.args(Path(temp))
            results.return_value.completed.return_value = Path('verified.hits')
            pool.run(args)
            command.assert_not_called()
            self.assertEqual(json.loads((args.out / 'pool_summary.json').read_text())['remaining'], 0)

    def test_completed_chunks_retire_pending_but_let_running_finish(self):
        with tempfile.TemporaryDirectory() as temp, patch.object(pool, 'Results') as results, \
             patch.object(pool, 'command', side_effect=['101', '102']) as command, \
             patch.object(pool, 'active_jobs', side_effect=[{'101': 'RUNNING', '102': 'PENDING'}, {'101': 'RUNNING'}, {}]), \
             patch.object(pool.subprocess, 'run') as cancel:
            args = self.args(Path(temp))
            cancel.return_value.returncode = 0
            results.return_value.completed.side_effect = [None, 'verified.hits']
            pool.run(args)
            self.assertEqual(cancel.call_args_list[0].args[0], ['scancel', '102'])
            self.assertEqual(cancel.call_count, 1)
            submitted = command.call_args_list[0].args[0]
            self.assertIn('--exclusive', submitted)
            self.assertIn('--mem=0', submitted)
            self.assertIn('--cpus-per-task=64', submitted)
            self.assertIn('--nodelist=node1,node2', submitted)

    def test_incomplete_pool_fails_and_preserves_results(self):
        with tempfile.TemporaryDirectory() as temp, patch.object(pool, 'Results') as results, \
             patch.object(pool, 'command', side_effect=['101', '102']), patch.object(pool, 'active_jobs', return_value={}):
            args = self.args(Path(temp))
            results.return_value.completed.return_value = None
            with self.assertRaisesRegex(RuntimeError, '1 chunks remain'):
                pool.run(args)
            self.assertEqual(json.loads((args.out / 'pool_summary.json').read_text())['remaining'], 1)

    def test_other_jobs_are_never_owned(self):
        with patch.object(pool, 'command', return_value='101|PENDING|different-name\n999|RUNNING|other\n102|PENDING|ours'):
            self.assertEqual(pool.active_jobs({'101': {'name': 'ours'}, '102': {'name': 'ours'}}), {'102': 'PENDING'})

    def test_interruption_cancels_only_owned_requests(self):
        with tempfile.TemporaryDirectory() as temp, patch.object(pool, 'Results') as results, \
             patch.object(pool, 'command', side_effect=['101', '102']), \
             patch.object(pool, 'active_jobs', return_value={'101': 'RUNNING', '102': 'PENDING'}), \
             patch.object(pool.time, 'sleep', side_effect=KeyboardInterrupt), patch.object(pool.subprocess, 'run') as cancel:
            results.return_value.completed.return_value = None
            cancel.return_value.returncode = 0
            with self.assertRaises(KeyboardInterrupt):
                pool.run(self.args(Path(temp)))
            self.assertEqual([call.args[0] for call in cancel.call_args_list], [['scancel', '101'], ['scancel', '102']])


if __name__ == '__main__':
    unittest.main()
