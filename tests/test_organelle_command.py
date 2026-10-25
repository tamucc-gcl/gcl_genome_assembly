from pathlib import Path
import sys
import unittest
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


if __name__=='__main__':unittest.main()
