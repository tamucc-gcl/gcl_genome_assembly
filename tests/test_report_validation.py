import json
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from validate_report import validate


class ReportValidation(unittest.TestCase):
    def fixture(self,root):
        assets=root/'reports/main';assets.mkdir(parents=True)
        (assets/'report_manifest.json').write_text('{}')
        (assets/'report_completeness.tsv').write_text('Feature\tEntity\tStatus\tReason\nQC\tsample\tDISABLED\tOff\n')
        (assets/'artifact_registry.tsv').write_text('type\tid\tpath\tnote\n')
        (root/'assembly_report.md').write_text('# Results\n\n[Scope](#hello-world)\n\n## Hello, world\n')

    def test_disabled_analysis_is_not_a_missing_output(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);self.fixture(root)
            self.assertEqual(validate(root)['status'],'PASS')

    def test_missing_expected_output_fails_even_with_valid_links(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);self.fixture(root)
            (root/'reports/main/report_completeness.tsv').write_text('Feature\tEntity\tStatus\nContact map\tsample\tMISSING\n')
            self.assertEqual(validate(root)['status'],'FAIL')

    def test_missing_file_and_wrong_anchor_fail(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);self.fixture(root)
            with (root/'assembly_report.md').open('a') as h:h.write('[Absent](absent.png)\n[Bad](#unknown)\n')
            result=validate(root)
            self.assertTrue(any('absent.png' in e for e in result['errors']))
            self.assertTrue(any('#unknown' in e for e in result['errors']))

    def test_registered_detail_reports_are_checked_too(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);self.fixture(root)
            (root/'detail.md').write_text('# Detail\n![Missing plot](missing.svg)\n')
            (root/'reports/main/artifact_registry.tsv').write_text('type\tid\tpath\tnote\nassembly_details\trun\tdetail.md\t\n')
            self.assertTrue(any('missing.svg' in e for e in validate(root)['errors']))

    def test_published_symlink_to_work_directory_is_valid(self):
        with tempfile.TemporaryDirectory() as temp, tempfile.TemporaryDirectory() as work:
            root=Path(temp);self.fixture(root)
            target=Path(work)/'figure.svg';target.write_text('<svg />')
            (root/'figure.svg').symlink_to(target)
            with (root/'assembly_report.md').open('a') as h:h.write('<img src="figure.svg" />\n')
            self.assertEqual(validate(root)['status'],'PASS')

    def test_disabled_report_does_not_claim_validation_pass(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);(root/'provenance').mkdir()
            (root/'provenance/run_provenance.json').write_text(json.dumps({'context':{'parameters':{'run_report':False}}}))
            self.assertEqual(validate(root)['status'],'DISABLED')


if __name__=='__main__':unittest.main()
