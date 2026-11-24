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

    def test_nested_gx_segments_validate_each_coordinate_frame(self):
        name='chr1_1~2034221..2314156~~1..149802'
        self.assertEqual(parent_sequence(name,{'chr1_1':3000000}),'chr1_1')
        # The first suffix may already be part of the literal screened ID.
        literal='chr1_1~2034221..2314156'
        self.assertEqual(parent_sequence(name,{literal:279936}),literal)
        self.assertEqual(parent_sequence(name,{name:500}),name)
        with self.assertRaisesRegex(ValueError,'coordinates outside'):
            parent_sequence(name,{'chr1_1':2000000})
        # Inner bounds must fit the 279,936 bp parent segment, not the chromosome.
        with self.assertRaisesRegex(ValueError,'coordinates outside'):
            parent_sequence(literal+'~~1..300000',{'chr1_1':3000000})

    def test_nested_taxonomy_stays_on_exact_contig_with_mixed_evidence(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            fasta=root/'a.fa';fasta.write_text('>chr1_1\n'+'ACGT'*25+'\n')
            taxonomy=root/'tax.rpt'
            taxonomy.write_text('#seq-id\tseq-len\tresult\tdiv\n'
                'chr1_1~11..90~~1..40\t40\tprimary-div\tanml:fishes\n'
                'chr1_1~11..90~~41..80\t40\tcontaminant(div)\tprok:bacteria\n')
            coverage=root/'coverage.txt'
            coverage.write_text('#rname\tstartpos\tendpos\tmeandepth\nchr1_1\t1\t100\t20\n')
            out=root/'evidence.tsv'
            build(fasta,taxonomy,coverage,'A','final',out)
            with out.open() as handle:rows=list(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual(len(rows),1)
            self.assertEqual(rows[0]['sequence'],'chr1_1')
            self.assertEqual(rows[0]['length_bp'],'100')
            self.assertEqual(rows[0]['plot_taxon'],'mixed FCS classifications')

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
