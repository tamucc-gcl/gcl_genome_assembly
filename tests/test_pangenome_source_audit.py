import csv
import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'py_scripts'))
from pangenome_source_audit import audit


class SourceAuditTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        (self.root/'a.fai').write_text('chr1\t100\nchr2\t80\n')
        self.ids = self.root/'ids.tsv'
        self.ids.write_text('id\tgraph_name\tis_reference\tsequence_count\tbases\n'
                            'sample\tIND_x.1\tFalse\t2\t180\n')
        self.paths = self.root/'paths.tsv'
        self.paths.write_text('assembly\tsequence\tfull_W_coordinate_span_bp\tclip_retained_path_bp\n'
                              'sample\tchr1\t100\t90\n')
        self.specs = [dict(id='sample', fai=str(self.root/'a.fai'), chromosomes=['chr1','chr2'])]

    def run_audit(self, evidence=()):
        return audit(self.specs,self.ids,self.paths,evidence,self.root/'out')

    def table(self, name):
        with (self.root/'out'/name).open() as f:
            return list(csv.DictReader(f,delimiter='\t'))

    def mapping(self, kind, intervals, query='id=IND_x.1|chr2', length=80):
        path = self.root/('map.'+kind+'.gz')
        with gzip.open(path,'wt') as f:
            for lo,hi,mq in intervals:
                f.write(f'{query}\t{length}\t{lo}\t{hi}\t+\tid=_MINIGRAPH_|s1\t100\t0\t{hi-lo}\t{hi-lo}\t{hi-lo}\t{mq}\n')
        return dict(kind=kind,path=str(path))

    def test_missing_mapping_is_unknown_not_measured_zero(self):
        result=self.run_audit()
        self.assertEqual(result['status'],'REVIEW')
        row=self.table('missing_source_paths.tsv')[0]
        self.assertEqual(row['status'],'missing_path_mapping_unavailable')
        self.assertEqual(row['gaf_records'],'')

    def test_mapping_localizes_loss_after_mapping_and_unions_overlap(self):
        result=self.run_audit([self.mapping('gaf',[(0,40,60),(20,60,255)])])
        row=self.table('missing_source_paths.tsv')[0]
        self.assertEqual(row['status'],'mapped_but_absent_from_full_and_clip')
        self.assertEqual(row['gaf_query_covered_bp'],'60')
        self.assertEqual(row['gaf_mapq20_records'],'1')
        self.assertEqual(result['chromosome_scale_paths_requiring_review'],1)

    def test_paf_is_separate_evidence_not_equated_with_alignment_partition(self):
        self.run_audit([self.mapping('paf',[(0,80,50)])])
        row=self.table('missing_source_paths.tsv')[0]
        self.assertEqual(row['status'],'paf_present_but_absent_from_full_and_clip')
        self.assertEqual(row['paf_query_covered_bp'],'80')

    def test_actual_zero_only_when_file_was_read(self):
        self.run_audit([self.mapping('gaf',[])])
        row=self.table('missing_source_paths.tsv')[0]
        self.assertEqual(row['gaf_records'],'0')
        self.assertEqual(row['paf_records'],'')

    def test_gaf_to_paf_discrepancy_is_distinct_from_later_graph_loss(self):
        self.run_audit([self.mapping('gaf',[(0,80,60)]),self.mapping('paf',[])])
        row=self.table('missing_source_paths.tsv')[0]
        self.assertEqual(row['status'],'gaf_present_but_no_paf_query_records')

    def test_synthetic_mapping_queries_are_not_unknown_biological_samples(self):
        result=self.run_audit([self.mapping('paf',[(0,80,60)],query='id=_MINIGRAPH_|s1')])
        self.assertEqual(result['unknown_mapping_name_prefixes'],{})
        self.assertEqual(result['mapping_evidence'][0]['synthetic_query_records'],1)

    def test_fai_and_audited_fasta_must_reconcile(self):
        (self.root/'a.fai').write_text('chr1\t100\nchr2\t79\n')
        with self.assertRaisesRegex(ValueError,'does not reconcile'):
            self.run_audit()

    def test_graph_names_cannot_silently_disappear_from_input_dictionary(self):
        self.paths.write_text(self.paths.read_text().replace('chr1','wrong'))
        with self.assertRaisesRegex(ValueError,'Unknown or duplicate graph source'):
            self.run_audit()

    def test_query_coordinates_must_match_source(self):
        with self.assertRaisesRegex(ValueError,'coordinates disagree'):
            self.run_audit([self.mapping('gaf',[(0,40,60)],length=81)])

    def test_complete_clip_representation_passes(self):
        with self.paths.open('a') as f:f.write('sample\tchr2\t80\t70\n')
        result=self.run_audit()
        self.assertEqual(result['status'],'PASS')
        self.assertEqual(result['source_names_missing_from_both_graphs'],0)

    def test_small_missing_sequences_still_reported_without_chromosome_alarm(self):
        self.specs[0]['chromosomes']=['chr1']
        result=self.run_audit()
        self.assertEqual(result['status'],'PASS')
        self.assertEqual(len(self.table('missing_source_paths.tsv')),1)

    def test_reference_query_absence_is_not_claimed_mapping_failure(self):
        self.ids.write_text(self.ids.read_text().replace('False','True'))
        self.run_audit([self.mapping('gaf',[])])
        row=self.table('missing_source_paths.tsv')[0]
        self.assertEqual(row['status'],'reference_path_missing_mapping_not_diagnostic')


if __name__ == '__main__':unittest.main()
