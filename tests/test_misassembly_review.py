import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'py_scripts'))
from misassembly_core import sha,write_table
from misassembly_report import FIELDS


class ReviewIntegrationTests(unittest.TestCase):
    def test_fresh_report_manual_gap_and_added_internal_cut(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);packet=root/'packet';packet.mkdir()
            fasta=root/'assessment.fa';fasta.write_text('>s\nAAAANNCCCC\n>other\nGGGGGGGG\n')
            names=root/'names.tsv';names.write_text('old_name\tnew_name\torient\tlength\tclass\ns\tchr1+chr2\t+\t10\tchromosome\nother\tchr3\t+\t8\tchromosome\n')
            peers=dict(separate=['p1','p2'],continuous=[],conflicting=[])
            option=dict(id='C01G01',scaffold='s',lo=4,hi=6,cut_bp=6,peers=peers,spanning=0,flank_molecules=[5,5],contacts=[])
            event=dict(id='C01',scaffold='s',lo=4,hi=6,left='chr1',right='chr2',recommendation='CUT',reason='Fixture: chromosome separation with local evidence.',peers=peers,cut_options=[option],recommended_cut=option,peer_tests=[])
            data=dict(assembly='synthetic',sample='synthetic',assessment_sha256=sha(fasta),status='assessed',mapping=dict(status='unavailable'),events=[event])
            (packet/'evidence.json').write_text(json.dumps(data))
            def run(script,*args,cwd=root):
                subprocess.run([sys.executable,str(ROOT/'py_scripts'/script),*map(str,args)],cwd=cwd,check=True,capture_output=True)
            run('misassembly_report.py','--packet',packet)
            report=(root/'README.md').read_text(encoding='utf-8')
            self.assertIn('Recommend cutting',report);self.assertIn('Evidence for breaking',report)
            self.assertNotIn('read_support.tsv',report)
            with (root/'review.tsv').open(encoding='utf-8') as handle:rows=list(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual(rows[0]['selected'],'NO');self.assertEqual(rows[0]['cut_bp'],'6')
            rows[0].update(selected='YES',review_disposition='CUT',reviewer='tester',reason='Approve fixture gap.')
            extra=dict(rows[0],id='manual01',scaffold='other',cut_bp='4',action='BREAK_PROBABLE_MISJOIN',gap_start='',gap_end='',localization_status='localized',reason='Manually added cut absent from discovery.')
            actions=root/'reviewed.tsv';write_table(actions,rows+[extra],FIELDS)
            run('break_chimeras.py','--fasta',fasta,'--name-map',names,'--actions',actions,'--assembly','synthetic','--out-fasta',root/'broken.fa','--out-name-map',root/'broken.names.tsv','--audit',root/'audit.tsv','--min-piece-bp',2)
            verification=json.loads((root/'broken.fa.verification.json').read_text())
            self.assertTrue(verification['exact_parent_reconstruction']);self.assertEqual(len(verification['selected_actions']),2)
            self.assertEqual(verification['total_bp'],18)
            second=root/'second';second.mkdir()
            run('misassembly_report.py','--packet',packet,'--actions',actions,'--verification',root/'broken.fa.verification.json',cwd=second)
            self.assertIn('2 actions verified',(second/'README.md').read_text(encoding='utf-8'))


if __name__=='__main__':unittest.main()
