from pathlib import Path
import sys
import unittest
import tempfile
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from organelle_command import validate_log


class OrganelleCommandTests(unittest.TestCase):
    def test_completed_no_recovery_is_not_an_execution_error(self):
        validate_log('Extracting animal_mt from the assemblies failed.\nThank you!\n')

    def test_logged_errors_fail_even_with_completion_marker(self):
        for text in ['2026-10-25 - ERROR: missing database\nThank you!',
                     'ERROR: dependency missing\nThank you!',
                     'Traceback (most recent call last):\nThank you!']:
            with self.assertRaisesRegex(ValueError,'reported an error'):validate_log(text)

    def test_missing_completion_marker_fails(self):
        with self.assertRaisesRegex(ValueError,'completion marker missing'):validate_log('Starting assembly\n')

    def test_warning_does_not_fail_completed_analysis(self):
        validate_log('2026-10-25 - WARNING: no circular path\nThank you!')

    def test_empty_retry_requires_verified_artifacts_and_no_errors(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)
            command=['get_organelle_from_assembly.py','-o',str(out)]
            log='2026-10-26 - INFO: Slimming fromasm_out/initial_assembly_graph.fastg finished!\n'
            with self.assertRaises(ValueError):validate_log(log,command)
            (out/'initial_assembly_graph.fastg').write_text('graph')
            (out/'slimmed_assembly_graph.fastg').write_text('')
            self.assertEqual(validate_log(log,command),'no_target_after_slimming')
            with self.assertRaises(ValueError):validate_log('Starting assembly',command)
            with self.assertRaises(ValueError):validate_log(log+'ERROR: broken\n',command)
            with self.assertRaises(ValueError):validate_log(log,['get_organelle_from_reads.py','-o',str(out)])
            (out/'slimmed_assembly_graph.fastg').write_text('not empty')
            with self.assertRaises(ValueError):validate_log(log,command)


if __name__=='__main__':unittest.main()
