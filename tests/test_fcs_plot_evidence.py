from pathlib import Path
import csv
import json
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from fcs_plot_evidence import build, classification, parent_sequence


class FcsPlotTests(unittest.TestCase):
    def test_inconclusive_is_not_promoted_to_host(self):
        self.assertEqual(classification([dict(result='repeat',div='anml:fishes')]),'unresolved: repeat')
        self.assertEqual(classification([dict(result='primary-div',div='anml:fishes')]),'anml:fishes')
        self.assertEqual(classification([]),'unassigned')
        self.assertEqual(classification([dict(result='primary-div',div='anml:fishes'),
            dict(result='contaminant(div)',div='prok:bacteria')]),'mixed FCS classifications')

    def test_exact_names_and_fcs_segment_suffixes(self):
        lengths={'contig':100,'literal~1..5':20}
        self.assertEqual(parent_sequence('literal~1..5',lengths),'literal~1..5')
        self.assertEqual(parent_sequence('contig~~10..20',lengths),'contig')
        self.assertEqual(parent_sequence('contig~1..20',lengths),'contig')
        with self.assertRaises(ValueError):parent_sequence('other~1..20',lengths)
        with self.assertRaises(ValueError):parent_sequence('contig~~10..200',lengths)

    def test_all_contigs_including_zero_depth_and_unassigned(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            fasta=root/'a.fa'; fasta.write_text('>c1\nACGT\n>c2\nNNNN\n')
            taxonomy=root/'tax.rpt'
            taxonomy.write_text('## metadata\n#seq-id\tseq-len\tresult\tdiv\nc1\t4\tprimary-div\tanml:fishes\n')
            coverage=root/'coverage.txt'
            coverage.write_text('#rname\tstartpos\tendpos\tmeandepth\nc1\t1\t4\t20\nc2\t1\t4\t0\n')
            output=root/'evidence.tsv'
            build(fasta,taxonomy,coverage,'A','final',output)
            with output.open() as handle: data=list(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual(len(data),2)
            self.assertEqual(data[0]['classification_source'],'FCS-GX')
            self.assertEqual(data[1]['plot_taxon'],'unassigned')
            self.assertEqual(data[1]['gc_fraction'],'NA')
            self.assertEqual(data[1]['mean_depth'],'0')
            self.assertEqual(json.loads(output.with_suffix('.json').read_text())['contigs'],2)
            coverage.write_text('#rname\tstartpos\tendpos\tmeandepth\nc1\t1\t3\t20\n')
            with self.assertRaisesRegex(ValueError,'coordinates differ'):
                build(fasta,taxonomy,coverage,'A','final',output)


if __name__=='__main__':unittest.main()
