import csv
import gzip
import io
import json
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'py_scripts'))
from pangenome_build_qc import summarize, union_bp


class BuildQCTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root/'ids').write_text('id\tgraph_name\tbases\nA1\tREF_A1\t10\nA2\tIND_A.2\t12\n')
        (self.root/'clip').write_text('nodes\t1\nedges\t0\nlength\t4\n')
        (self.root/'oldversions').write_text('tool\tversion\ncactus\t9.4.1\n')
        with gzip.open(self.root/'full.gz','wt') as handle:
            handle.write('S\t1\tACGT\nS\t2\t*\tLN:i:6\nW\tREF_A1\t0\tchr1\t0\t4\t>1\n')
        self.files = {'sample-stats.tsv':'chr1\tREF_A1#0\t9\nchr1\tIND_A#2\t8\n',
            'path-stats.tsv':'chr1\tREF_A1#0#chr1\t9\nchr1\tIND_A#2#chr1[0-8]\t8\n',
            'clip-stats.bed':'IND_A#2#chr1\t1\t4\tchr1\nIND_A#2#chr1\t3\t5\tchr1\n'}

    def run_summary(self):
        with tarfile.open(self.root/'stats.tgz','w:gz') as archive:
            for name,text in self.files.items():
                data=text.encode()
                member=tarfile.TarInfo('tool/work/'+name);member.size=len(data)
                archive.addfile(member,io.BytesIO(data))
        return summarize(self.root/'stats.tgz', self.root/'full.gz', self.root/'clip',
                         self.root/'ids', self.root/'out', 'pinned-image', self.root/'oldversions')

    def test_cactus_stats_and_full_graph_accounting(self):
        result=self.run_summary()
        self.assertEqual(result['status'],'PASS')
        self.assertEqual(result['full']['sequence_bp'],10)
        with (self.root/'out/assembly_clipping.tsv').open() as handle:
            rows={r['assembly']:r for r in csv.DictReader(handle,delimiter='\t')}
        self.assertEqual(rows['A2']['reported_clip_interval_bp'],'4')
        self.assertEqual(rows['A2']['retained_path_bp'],'8')
        self.assertEqual(rows['A2']['input_minus_retained_minus_reported_clip_bp'],'0')
        self.assertIn('9.4.1',json.loads((self.root/'out/tool_provenance.json').read_text())['original_build_version_record'])

    def test_missing_bundle_member_fails(self):
        del self.files['sample-stats.tsv']
        with self.assertRaisesRegex(ValueError,'bundle lacks'):
            self.run_summary()

    def test_unmapped_names_are_explicit(self):
        self.files['clip-stats.bed'] += 'unknown#0#chr1\t0\t2\tchr1\n'
        result=self.run_summary()
        self.assertEqual(result['status'],'name_mapping_incomplete')
        self.assertEqual(result['unknown_name_categories'],1)

    def test_invalid_interval_fails(self):
        self.files['clip-stats.bed']='IND_A#2#chr1\t8\t1\n'
        with self.assertRaisesRegex(ValueError,'Invalid clipped interval'):
            self.run_summary()

    def test_union_does_not_double_count_overlapping_intervals(self):
        self.assertEqual(union_bp([(2,5),(0,3),(10,12)]),7)

    def test_path_and_sample_counts_must_agree(self):
        self.files['path-stats.tsv']=self.files['path-stats.tsv'].replace('\t8\n','\t6\n')
        with self.assertRaisesRegex(ValueError,'retained counts disagree'):
            self.run_summary()

    def test_missing_full_walk_is_not_assumed_zero(self):
        self.run_summary()
        with (self.root/'out/full_clip_sequence_paths.tsv').open() as handle:
            rows={r['assembly']:r for r in csv.DictReader(handle,delimiter='\t')}
        self.assertEqual(rows['A2']['full_W_coordinate_span_bp'],'')

    def test_full_composite_missing_from_clip_has_explicit_evidence(self):
        with gzip.open(self.root/'full.gz','at') as handle:
            handle.write('W\tIND_A\t2\tchr4_1+chr14_1\t0\t6\t>2\n')
        self.run_summary()
        with (self.root/'out/full_clip_sequence_paths.tsv').open() as handle:
            row=next(r for r in csv.DictReader(handle,delimiter='\t') if '+' in r['sequence'])
        self.assertEqual(row['full_W_coordinate_span_bp'],'6')
        self.assertEqual(row['clip_retained_path_bp'],'0')
        self.assertIn('chr4_1+chr14_1', (self.root/'out/clipping_report.md').read_text())

    def test_ambiguous_bundle_names_fail(self):
        self.files['other/sample-stats.tsv']=self.files['sample-stats.tsv']
        with self.assertRaisesRegex(ValueError,'Duplicate clipping bundle member'):
            self.run_summary()


if __name__ == '__main__':
    unittest.main()
