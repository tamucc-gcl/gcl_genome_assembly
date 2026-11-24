"""Nextflow stages report directories as symlinks; the overview must read those inputs."""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

HELPER=Path(__file__).resolve().parents[1]/'py_scripts/fcs_overview.py'


class OverviewInputs(unittest.TestCase):
    @unittest.skipIf(os.name=='nt','Exercise Linux Nextflow directory symlinks')
    def test_symlinked_stage_directories_produce_both_stage_rows(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);task=root/'task';task.mkdir()
            inputs=[]
            for stage in ('contig','final'):
                report=root/stage/'report';report.mkdir(parents=True)
                (report/'cohort_fcs_sequence_evidence.tsv').write_text(
                    'cohort_assembly\tcoordinate_stage\tsequence\tcat_primary_bp\n'
                    f'A\t{stage}\tc1\t4\n')
                staged=task/'inputs'/stage/'report';staged.parent.mkdir(parents=True)
                staged.symlink_to(report,target_is_directory=True);inputs.append(str(staged))
            result=subprocess.run([sys.executable,str(HELPER),*inputs],cwd=task,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            combined=task/'overview/cohort_fcs_sequence_evidence.tsv'
            self.assertTrue(combined.is_file(),'Overview silently ignored the staged directory symlinks')
            with combined.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
            self.assertEqual({r['coordinate_stage'] for r in rows},{'contig','final'})
            self.assertEqual(len(rows),2)
            scope=json.loads((task/'overview/overview_scope.json').read_text())
            self.assertEqual(scope['records_by_stage'],{'contig':1,'final':1})

    def test_missing_stage_table_fails_instead_of_claiming_success(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            result=subprocess.run([sys.executable,str(HELPER),str(root/'absent')],
                                  cwd=root,capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertFalse((root/'overview/contamination_report.md').exists())


if __name__=='__main__':unittest.main()
