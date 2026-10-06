"""Test assay orchestration with mocked external tools, not synthetic biological validation."""
import csv
import hashlib
import json
import sys
import tempfile
import unittest
import shutil
from pathlib import Path
from unittest.mock import patch, MagicMock
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
import chimera_sequence_context as context
from break_chimeras import fasta


class Orchestration(unittest.TestCase):
    def test_existing_bam_and_library_assay_write_complete_packet_without_source_indices(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);assessment=root/'input.fa';source=root/'source.fa'
            assessment.write_text('>s\n'+'A'*1000+'N'*100+'C'*1000+'\n')
            source.write_text('>source\n'+'A'*1000+'C'*1000+'\n')
            row=dict(assembly='sample',scaffold='s',assembly_sha256=hashlib.sha256(assessment.read_bytes()).hexdigest(),
                     coordinate_stage='pre_finishing',chromosome_member='yes',candidate_verdict='REVIEW',
                     callable='yes',evidence_only='no',transition_lo='1000',transition_hi='1100')
            calls=root/'calls.tsv'
            with calls.open('w') as handle:
                w=csv.DictWriter(handle,fieldnames=list(row),delimiter='\t');w.writeheader();w.writerow(row)
            (root/'last.agp').write_text('s\t1\t1000\t1\tW\tsource\t1\t1000\t+\n'
                's\t1001\t1100\t2\tN\t100\tscaffold\tyes\tproximity_ligation\n'
                's\t1101\t2100\t3\tW\tsource\t1001\t2000\t+\n')
            (root/'pairs').write_text('#chromsize: source 2000\nA_read\tsource\t10\tsource\t1001\t+\t-\tUU\n'
                                      'B_read\tsource\t10\tsource\t1001\t+\t-\tUU\n')
            (root/'libraries.tsv').write_text('library_id\tqname_prefix\nA\tA_\nB\tB_\n')
            (root/'reads.bam').write_text('Mock BAM; not biological data')
            (root/'bam.json').write_text(json.dumps(dict(sha256=row['assembly_sha256'],coordinate_stage='pre_finishing')))
            def fake_run(args,output=None):
                args=list(map(str,args))
                if args[:2]==['samtools','faidx']:
                    records=fasta(args[2])
                    if len(args)==3:
                        Path(args[2]+'.fai').write_text(''.join('%s\t%d\t0\t60\t61\n'%(name,len(seq)) for name,seq in records.items()))
                    else:
                        name,region=args[3].rsplit(':',1);lo,hi=map(int,region.split('-'))
                        Path(output).write_text('>'+name+'\n'+records[name][lo-1:hi]+'\n')
                elif args[:2]==['samtools','idxstats']:
                    Path(output).write_text('s\t2100\t30\t0\n*\t0\t0\t0\n')
                elif args[:2]==['samtools','view']:
                    records=[]
                    for side,start in (('left',1),('right',1101)):
                        for n in range(15):records.append('%s%d\t0\ts\t%d\t60\t1000M\t*\t0\t0\t*\t*\tNM:i:0\n'%(side,n,start))
                    Path(output).write_text(''.join(records))
                elif output:
                    Path(output).write_text('Mock external tool output\n')
            out=root/'packet'
            argv=['context','--fasta',str(assessment),'--calls',str(calls),'--assembly','sample','--sample','individual',
                  '--peers','[]','--motif','CCCTAA','--bam',str(root/'reads.bam'),'--bam-provenance',str(root/'bam.json'),
                  '--agp',str(root/'last.agp'),'--pairs',str(root/'pairs'),'--pairs-source',str(source),
                  '--libraries',str(root/'libraries.tsv'),'--out',str(out)]
            plot=MagicMock();fig=MagicMock();axes=MagicMock()
            plot.subplots.return_value=(fig,axes)
            fig.savefig.side_effect=lambda path,**kwargs:Path(path).write_text('Mock plot, not scientific output')
            matplotlib=MagicMock();matplotlib.pyplot=plot
            # Avoid requiring Windows symlink privileges for a Linux orchestration test.
            with patch.object(context,'run',fake_run),patch.object(sys,'argv',argv),patch.object(Path,'symlink_to',lambda self,target:shutil.copyfile(target,self)),patch.dict(sys.modules,{'matplotlib':matplotlib,'matplotlib.pyplot':plot}):
                context.main()
            measurements=json.loads((out/'decision_measurements.json').read_text())
            value=next(iter(measurements.values()))
            self.assertEqual((value['gap_start'],value['gap_end'],value['cut_bp']),(1000,1100,1100))
            self.assertEqual(value['hifi_spanning_molecules'],0)
            self.assertTrue(value['hifi_informative'])
            self.assertFalse(value['matched_controls_pass'])
            self.assertFalse(value['hic_informative'])
            self.assertEqual({r['library'] for r in value['libraries']},{'A','B'})
            self.assertTrue((out/'candidate_1.controls.png').is_file())
            self.assertTrue((out/'candidate_1.igv.xml').is_file())
            self.assertFalse(Path(str(source)+'.fai').exists())
            self.assertFalse(Path(str(assessment)+'.fai').exists())


if __name__=='__main__':unittest.main()
