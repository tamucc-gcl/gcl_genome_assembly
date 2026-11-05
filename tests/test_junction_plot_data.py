import json
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from misassembly_core import write_table,table
from junction_plot_data import export


class PlacementPlotData(unittest.TestCase):
    def test_retained_counts_missing_gap_and_other_haplotype(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); hifi=root/'junction/hifi';hifi.mkdir(parents=True)
            write_table(hifi/'original_bridges.tsv',[
                dict(event='E1',site='S1',anchor_bp=1000,mapq0=3,mapq30=2),
                dict(event='E1',site='G1',anchor_bp=1000,mapq0=0,mapq30=0)],
                ['event','site','anchor_bp','mapq0','mapq30'])
            base=dict(event='E1',site='S1',mapq=60,aligned=10000,clipped=0,identity=.99,
                      same_haplotype_alternatives=0,other_haplotype_alternatives=1)
            write_table(hifi/'competitive_bridges.tsv',[base,dict(base,same_haplotype_alternatives=1)],list(base))
            site=dict(id='S1',point=12000,credible_bridges=1,unique_arm_bridges=1,
                      minimum_anchor_span_bp=5000,local_reads_long_enough='12')
            event=dict(id='E1',junction_evidence=dict(hifi_status='measured',sites=[site,
                dict(site,id='G1',point=18000,credible_bridges=None,unique_arm_bridges=None)]))
            before=json.dumps(event,sort_keys=True)
            self.assertTrue(export(root,event))
            rows=table(root/'junction/E1.placement_counts.tsv')
            self.assertEqual((rows[0]['credible'],rows[0]['competing'],rows[0]['other_haplotype']),('1','1','2'))
            self.assertEqual((rows[1]['credible'],rows[1]['primary_all']),('','0'))
            self.assertEqual(before,json.dumps(event,sort_keys=True))
            event['junction_evidence']['sites'][0]['credible_bridges']=2
            with self.assertRaisesRegex(ValueError,'disagrees'):export(root,event)

    def test_unavailable_is_not_rendered_as_zero(self):
        self.assertFalse(export('.',dict(id='E1')))


if __name__=='__main__':unittest.main()
