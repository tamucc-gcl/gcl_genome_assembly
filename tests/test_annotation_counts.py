import sys
import tempfile
import unittest
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from annotation_counts import audit


class AnnotationCounts(unittest.TestCase):
    def check(self,features):
        record=SeqRecord(Seq('A'*1000),id='mito',annotations={'molecule_type':'DNA'})
        record.features=features
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'mito.gb';SeqIO.write(record,path,'genbank')
            return audit(path,'sample')

    def test_gene_and_rna_features_are_one_locus_not_two(self):
        loc=SimpleLocation(0,70,strand=1)
        summary,rows=self.check([SeqFeature(loc,type='gene',qualifiers={'product':['tRNA-Phe']}),
                                 SeqFeature(loc,type='tRNA',qualifiers={'product':['tRNA-Phe']})])
        self.assertEqual(summary['annotated_loci'],1)
        self.assertEqual(summary['trna_loci'],1)
        self.assertEqual(rows[0]['annotation_feature_rows'],2)

    def test_duplicate_labels_at_distinct_loci_are_preserved(self):
        summary,_=self.check([SeqFeature(SimpleLocation(lo,lo+60),type='tRNA',qualifiers={'gene':['trnL']}) for lo in (0,100)])
        self.assertEqual(summary['trna_loci'],2)

    def test_circular_compound_location_counts_once(self):
        loc=CompoundLocation([SimpleLocation(900,1000,strand=1),SimpleLocation(0,40,strand=1)])
        summary,_=self.check([SeqFeature(loc,type='gene'),SeqFeature(loc,type='CDS')])
        self.assertEqual(summary['annotated_loci'],1)
        self.assertEqual(summary['protein_coding_loci'],1)

    def test_empty_annotation_is_not_forced_to_expected_gene_total(self):
        summary,rows=self.check([])
        self.assertEqual(summary['annotated_loci'],0)
        self.assertEqual(rows,[])


if __name__=='__main__':unittest.main()
