"""Check recovery ordering with isolated Slurm command doubles; no real jobs."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


@unittest.skipUnless(os.name == 'posix', 'Recovery uses Bash and Slurm commands')
class RecoveryTests(unittest.TestCase):
    def run_recovery(self, controller_stops):
        script = Path(__file__).resolve().parents[1] / 'scripts/comparisons/recover_contamination_staged_cohort.sh'
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tools = root / 'bin'
            tools.mkdir()
            state = root / 'state.json'
            state.write_text(json.dumps(dict(alive=True, events=[])))
            source = '''import json,os,pathlib,sys
p=pathlib.Path(os.environ['RECOVERY_STATE']); s=json.loads(p.read_text())
name=pathlib.Path(sys.argv[0]).name; a=sys.argv[1:]
if name=='squeue':
    job=a[a.index('-j')+1]; field=a[a.index('-o')+1]
    if field=='%u': print(os.environ['USER'])
    elif field=='%j': print('contam_cohort' if job=='123' else 'nf-CONTAMINATION_MEGABLAST_NODE')
    elif field=='%T': print('PENDING')
    elif field=='%i' and s['alive']: print('123')
elif name=='python3': s['events'].append('verified')
elif name=='scancel':
    s['events'].append('cancel:'+','.join(a))
    if a==['123'] and os.environ['CONTROLLER_STOPS']=='yes': s['alive']=False
elif name=='sbatch':
    if s['alive']: raise SystemExit('Submitted while controller alive')
    s['events'].append('submit'); print('789')
elif name=='sleep': pass
p.write_text(json.dumps(s))
'''
            for name in ('squeue', 'scancel', 'sbatch', 'python3', 'sleep'):
                path = tools / name
                path.write_text('#!' + sys.executable + '\n' + source)
                path.chmod(0o755)
            startup = root / 'startup.sh'
            startup.write_text('module() { :; }\n')
            result = subprocess.run(['bash', str(script), '123', '456', 'unused-plan.json'], cwd=root,
                                    env=dict(os.environ, PATH=str(tools) + ':' + os.environ['PATH'],
                                             BASH_ENV=str(startup), RECOVERY_STATE=str(state),
                                             CONTROLLER_STOPS='yes' if controller_stops else 'no'),
                                    capture_output=True, text=True)
            return result, json.loads(state.read_text())['events']

    def test_job_cancellation_precedes_child_cancellation_and_submission(self):
        result, events = self.run_recovery(True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(events, ['verified', 'cancel:123', 'cancel:456', 'submit'])

    def test_live_controller_blocks_child_cancellation_and_resubmission(self):
        result, events = self.run_recovery(False)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(events, ['verified', 'cancel:123'])


if __name__ == '__main__':
    unittest.main()
