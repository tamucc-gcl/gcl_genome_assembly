import unittest,tempfile,sys,csv,json,hashlib,subprocess
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from chimera_review import generate
class ManualReview(unittest.TestCase):
    def test_packet_defaults_every_cut_to_no_including_former_auto_eligible(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);ctx=root/'context';ctx.mkdir();calls=root/'calls.tsv'
            row=dict(scaffold='s',assembly_sha256='a'*64,transition_lo='5',transition_hi='6')
            with calls.open('w',newline='') as h:
                w=csv.DictWriter(h,fieldnames=list(row),delimiter='\t');w.writeheader();w.writerow(row)
            key=hashlib.sha256(json.dumps(row,sort_keys=True).encode()).hexdigest()[:20]
            (ctx/'decision_measurements.json').write_text(json.dumps({key:dict(verified_gap=True,gap_start=5,gap_end=6,cut_bp=6,auto_eligible=True)}))
            generate('assembly',calls,ctx,root/'review')
            with (root/'review/review.tsv').open() as h:result=list(csv.DictReader(h,delimiter='\t'))
            self.assertEqual(result[0]['selected'],'NO');self.assertEqual(result[0]['reviewer'],'');self.assertEqual(result[0]['reason'],'')
            self.assertIn('No automatic cutting',(root/'review/report.html').read_text())
    def test_standalone_automatic_mode_is_disabled(self):
        script=Path(__file__).resolve().parents[1]/'py_scripts/break_chimeras.py'
        result=subprocess.run([sys.executable,str(script),'--fasta','missing','--name-map','missing','--actions','missing',
            '--assembly','a','--out-fasta','unused','--out-name-map','unused','--audit','unused','--mode','auto'],capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0);self.assertIn('Automated cutting is deferred',result.stderr)
    def test_unselected_bad_coordinates_ignored_but_selected_review_required(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);fa=root/'input.fa';fa.write_text('>s\nAAAANNCCCC\n')
            names=root/'names.tsv';names.write_text('old_name\tnew_name\torient\torder\tlength\tclass\tref_span\tflags\n'+'s\tchr1\t+\t1\t10\tchromosome\t.\t.\n')
            row=dict(selected='NO',id='r',assembly='a',coordinate_stage='pre_finishing',assessment_sha256='stale',scaffold='s',action='UNRESOLVED',cut_bp='',gap_start='',gap_end='',decision_source='review',evidence_packet_id='p',reviewer='',reason='',localization_status='unlocalized')
            path=root/'review.tsv'
            def write():
                with path.open('w',newline='') as h:
                    w=csv.DictWriter(h,fieldnames=list(row),delimiter='\t');w.writeheader();w.writerow(row)
            script=Path(__file__).resolve().parents[1]/'py_scripts/break_chimeras.py'
            args=[sys.executable,str(script),'--fasta',str(fa),'--name-map',str(names),'--actions',str(path),'--assembly','a','--out-fasta',str(root/'out.fa'),'--out-name-map',str(root/'out.tsv'),'--audit',str(root/'audit.tsv'),'--min-piece-bp','4']
            write();subprocess.run(args,check=True,capture_output=True);self.assertEqual((root/'out.fa').read_bytes(),fa.read_bytes())
            row.update(selected='YES');write();result=subprocess.run(args,capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0);self.assertIn('requires reviewer',result.stderr)
