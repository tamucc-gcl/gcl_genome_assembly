import json
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from fcs_category_evidence import enrich,CATEGORIES,table,location

class CategoryEvidence(unittest.TestCase):
    def test_nested_partition_retains_uncertainty_gaps_and_diagnostics(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);tax=root/'tax.rpt'
            tax.write_text('##'+json.dumps([['GX'],{'run-info':{'asserted-div':'anml:fishes'}}])+'\n'+
              '#seq-id\tseq-len\t(xp,lc,co,n)-len\tcvg-by-all\tresult\tdiv\n'+
              'chr~11..90~~1..40\t40\t10,5,0,0\t20\trepeat\tanml:fishes\n'+
              'chr~11..90~~41..80\t40\t0,0,0,0\t15\tlow-coverage\tanml:molluscs\n')
            seq=[dict(cohort_assembly='A',sequence='chr',length_bp='100',plot_taxon='mixed FCS classifications')]
            segments=enrich(seq,{'A':{'taxonomy':str(tax)}},root,'final')
            self.assertEqual(seq[0]['cat_repeat_expected_bp'],40)
            self.assertEqual(seq[0]['cat_repeat_other_bp'],40)
            self.assertEqual(seq[0]['cat_unreported_bp'],20)
            self.assertEqual(sum(seq[0]['cat_'+c+'_bp'] for c in CATEGORIES),100)
            self.assertEqual(sum(r.get('alignment_bp',0) for r in segments if isinstance(r.get('alignment_bp'),int)),35)
            self.assertEqual(segments[1]['result'],'repeat')
            self.assertEqual(len(table(root/'assembly_0001_selected.tsv')),1)
            # Duplicate/overlapping calls must never inflate composition.
            with tax.open('a') as f:f.write('chr~11..90~~1..40\t40\t0,0,0,0\t1\trepeat\tanml:fishes\n')
            with self.assertRaisesRegex(ValueError,'Overlapping'):
                enrich(seq,{'A':{'taxonomy':str(tax)}},root,'final')

    def test_coordinate_like_literal_name_is_not_reinterpreted(self):
        self.assertEqual(location('chr~11..90~~1..40',{'chr~11..90':80}),('chr~11..90',0,40))
        with self.assertRaises(ValueError):location('chr~11..90~~1..81',{'chr':100})

if __name__=='__main__':unittest.main()
