import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class RealVerification(unittest.TestCase):
    def fixture(self, root, extra=False, pending=False):
        (root/'original.fa').write_text('>parent\nAAAANNCCCC\n>other\nGGGG\n')
        (root/'corrected.fa').write_text('>parent_sub_0_6\nAAAANN\n>parent_sub_6_10\nCCCC\n>other\nGGGG\n'+('>unknown\nA\n' if extra else ''))
        rows=[dict(old_name='parent_sub_0_6',new_name='chr9_1',**{'class':'chromosome'},flags=''),
              dict(old_name='parent_sub_6_10',new_name='chr7_1+chr12_1',**{'class':'composite'},flags='chromosome_assignment_pending' if pending else ''),
              dict(old_name='other',new_name='chr4_1',**{'class':'chromosome'},flags='')]
        if extra:rows.append(dict(old_name='unknown',new_name='unplaced_1',**{'class':'unplaced'},flags=''))
        with (root/'names.tsv').open('w') as handle:
            w=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
        (root/'actions.json').write_text(json.dumps([dict(assembly='sample',scaffold='parent',cut_bp=6,gap_start=4,gap_end=6)]))

    def run_check(self, root):
        script=Path(__file__).resolve().parent/'verify_chimera_real_test.py'
        return subprocess.run([sys.executable,str(script),'--original',str(root/'original.fa'),'--corrected',str(root/'corrected.fa'),
            '--names',str(root/'names.tsv'),'--assembly','sample','--scenario','manual','--actions',str(root/'actions.json'),'--out',str(root/'PASS.json')],capture_output=True,text=True)

    def test_reassigned_composite_is_preserved_and_exactly_verified(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);self.fixture(root)
            result=self.run_check(root);self.assertEqual(result.returncode,0,result.stderr)
            qc=json.loads((root/'PASS.json').read_text());self.assertEqual(qc['pieces'][1]['new_name'],'chr7_1+chr12_1')
            self.assertEqual(qc['output_bp'],14)

    def test_extra_output_record_and_pending_assignment_are_rejected(self):
        for options in (dict(extra=True),dict(pending=True)):
            with self.subTest(options=options),tempfile.TemporaryDirectory() as directory:
                root=Path(directory);self.fixture(root,**options)
                self.assertNotEqual(self.run_check(root).returncode,0)


if __name__=='__main__':unittest.main()
