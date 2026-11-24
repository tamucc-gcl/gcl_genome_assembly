"""Exercise the acceptance helper through Nextflow-style scratch symlinks."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class ScratchImports(unittest.TestCase):
    def test_staged_cli_imports_without_pythonpath_or_source_checkout_modules(self):
        root=Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as temp:
            scratch=Path(temp)/'node scratch'
            scratch.mkdir()
            for relative in ('tests/fcs_adaptor_acceptance.py','py_scripts/contamination_review.py',
                             'py_scripts/contamination_cohort_report.py'):
                source=root/relative
                (scratch/source.name).symlink_to(source)
            # Isolated Python ignores PYTHONPATH and automatic cwd imports. Only the
            # helper's explicit task-directory import can find both staged modules.
            result=subprocess.run([sys.executable,'-I','fcs_adaptor_acceptance.py','--help'],
                                  cwd=scratch,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            self.assertIn('{fetch,check}',result.stdout)


if __name__=='__main__':unittest.main()
