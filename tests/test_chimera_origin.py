"""Synthetic origin/support tests; no real genomes or external tools required."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import gzip

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'py_scripts'))
from chimera_origin import (agp_rows, component_interval, exact_projection, exact_span,
                            read_support, revcomp, sha, source_to_query, table)

spec = importlib.util.spec_from_file_location('origin_runner', ROOT/'scripts/comparisons/trace_chimera_origins.py')
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


def hit(length=20, strand='+', cigar='20=', start=100):
    return dict(query='interval_001', query_length=length, query_start=0, query_end=length,
                strand=strand, target='old', target_start=start, target_end=start+length,
                target_length=1000, mapq=60, cigar=cigar)


class OriginTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def test_forward_and_reverse_exact_projection(self):
        self.assertEqual(exact_projection(hit(), 3, 8), (103, 108))
        self.assertEqual(exact_projection(hit(strand='-'), 3, 8), (112, 117))

    def test_insertions_are_not_bridged(self):
        h = hit(cigar='8=2I10=')
        h['target_end'] = 118
        self.assertIsNone(exact_projection(h, 5, 15))
        self.assertEqual(exact_projection(h, 12, 15), (110, 113))

    def test_reverse_indel_projection(self):
        h = hit(strand='-', cigar='8=2I10=')
        h['target_end'] = 118
        self.assertEqual(exact_projection(h, 3, 7), (111, 115))
        self.assertIsNone(exact_projection(h, 7, 15))

    def test_cigar_length_mismatch_rejected(self):
        with self.assertRaisesRegex(ValueError, 'CIGAR disagrees'):
            exact_projection(hit(cigar='19='), 1, 3)

    def test_literal_sequence_required_not_paf_outer_span(self):
        h = hit(length=4, cigar='4=')
        self.assertFalse(exact_span(h, 'ACGT', lambda *a: 'ACGA'))
        self.assertTrue(exact_span(h, 'ANNT', lambda *a: 'ANNT'))
        self.assertTrue(exact_span(hit(length=4, strand='-', cigar='4='), 'ACCC', lambda *a: 'GGGT'))

    def test_verified_reverse_boundary_transform(self):
        self.assertEqual(source_to_query(hit(strand='-'), 102, 107), (13, 18))
        self.assertEqual(source_to_query(hit(strand='-'), 107, 107), (13, 13))

    def test_agp_reverse_component_with_offset(self):
        p = self.root/'a.agp'
        p.write_text('s\t1\t10\t1\tW\tc\t21\t30\t-\n')
        self.assertEqual(component_interval(agp_rows(p)[0], 2, 7), (23, 28))

    def test_agp_unknown_orientation_and_noncontiguity_rejected(self):
        p = self.root/'a.agp'
        for line in ('s\t1\t10\t1\tW\tc\t1\t10\t?\n',
                     's\t2\t10\t1\tW\tc\t1\t9\t+\n'):
            p.write_text(line)
            with self.assertRaises(ValueError):
                agp_rows(p)

    def test_table_rejects_shifted_row(self):
        p = self.root/'a.tsv'
        p.write_text('a\tb\n1\n')
        with self.assertRaises(ValueError):
            table(p)

    @staticmethod
    def sam(name='r', cigar='5000M', mq=60, nm=0, flag=0):
        return '\t'.join([name, str(flag), 's', '1', str(mq), cigar, '*', '0', '0', '*', '*', f'NM:i:{nm}'])+'\n'

    def test_midpoint_spanning_does_not_span_interval_anchors(self):
        result = read_support([self.sam()], 's', 1000, 6000)
        self.assertEqual(result['molecules_passing_chain_screen'], 0)
        self.assertEqual(result['local_read_length_ge_anchor_span'], 1)
        self.assertEqual(result['expected_bridge_count'], 'not_estimated')

    def test_large_deletion_is_not_read_support(self):
        result = read_support([self.sam(cigar='2000M1000D2000M')], 's', 1000, 4000)
        self.assertEqual(result['molecules_bracketing_anchors'], 1)
        self.assertEqual(result['molecules_passing_chain_screen'], 0)

    def test_read_names_deduplicated_and_flags_filtered(self):
        records = [self.sam(), self.sam(), self.sam('supp', flag=2048), self.sam('dup', flag=1024),
                   self.sam('unknown', mq=255), self.sam('errors', nm=500), self.sam('low', mq=1)]
        result = read_support(records, 's', 1000, 4000)
        self.assertEqual(result['molecules_passing_chain_screen'], 1)
        self.assertEqual(result['primary_molecules_observed'], 4)

    def test_partner_end_order_and_half_open_coordinates(self):
        path = self.root/'p.gz'
        with gzip.open(path, 'wt') as handle:
            handle.write('#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\n')
            handle.write('a\ts\t10\tx\t5\t+\t-\tUU\n')  # 9 inside left [0,10)
            handle.write('b\tx\t5\ts\t21\t+\t-\tUU\n')  # 20 inside right [20,30)
            handle.write('c\ts\t11\tx\t5\t+\t-\tUU\n')  # 10 outside
        out = self.root/'partners.tsv'
        runner.partner_summary(path, {'i': dict(scaffold='s', lo=10, hi=20)}, {'s': 30}, 10, out)
        rows = table(out)
        self.assertEqual([(r['side'], r['contacts']) for r in rows], [('left', '1'), ('right', '1')])

    def test_cached_task_ambiguity_not_arbitrarily_selected(self):
        trace = self.root/'trace.tsv'
        trace.write_text('name\thash\nW:HIFIASM (sample)\tab/cdef\n')
        (self.root/'work/ab/cdef001').mkdir(parents=True)
        (self.root/'work/ab/cdef002').mkdir()
        runner.cached_inventory(trace, self.root/'work', self.root/'inventory.tsv')
        self.assertEqual(table(self.root/'inventory.tsv')[0]['status'], 'missing_or_ambiguous_task_directory')

    def test_end_to_end_exact_source_gap_and_missing_stage(self):
        assessment_path, raw = self.root/'assessment.fa', self.root/'raw.fa'
        assessment_path.write_text('>s\nAAAACCCCNNGGGGTTTT\n')
        raw.write_text('>old\nAAAACCCCNNGGGGTTTT\n')
        agp = self.root/'assessment.agp'
        agp.write_text('s\t1\t8\t1\tW\tc1\t1\t8\t+\n'
                       's\t9\t10\t2\tU\t2\tscaffold\tyes\tproximity_ligation\n'
                       's\t11\t18\t3\tW\tc2\t1\t8\t+\n')
        args = SimpleNamespace(results=self.root, samtools='unused', minimap2='unused',
                               flank=4, anchor=2, contact_flank=10, hifi_results=None)
        calls = [dict(assembly_sha256=sha(assessment_path), scaffold='s', transition_lo='6', transition_hi='12')]
        stages = [dict(stage='raw_contigs', fasta=str(raw)),
                  dict(stage='missing', fasta=str(self.root/'absent.fa')),
                  dict(stage='assessment', fasta=str(assessment_path), agp=str(agp))]

        class FakeFasta:
            def __init__(self, path, *unused):
                lines = Path(path).read_text().splitlines()
                self.seq = {lines[0][1:]: lines[1]}
                self.lengths = {k: len(v) for k, v in self.seq.items()}
            def fetch(self, name, lo, hi):
                return self.seq[name][lo:hi]

        def fake_map(args, target, query, prefix):
            if query.name == 'anchors.fa':
                return []
            h = hit(length=14, cigar='14=', start=2)
            return [h]

        with patch.object(runner, 'Fasta', FakeFasta), patch.object(runner, 'map_queries', fake_map):
            result = runner.process_assembly(args, 'sample_hap1', calls, stages, self.root/'out')
        self.assertEqual(result['earliest_exact_stage']['interval_001'], 'raw_contigs')
        ledger = table(self.root/'out/origin_ledger.tsv')
        self.assertEqual([r['status'] for r in ledger], ['exact_sequence_present', 'stage_missing', 'exact_sequence_present'])
        gaps = [r for r in table(self.root/'out/source_boundaries.tsv') if r['kind'] == 'agp_gap']
        self.assertEqual((gaps[0]['assessment_start'], gaps[0]['assessment_end']), ('8', '10'))
        self.assertTrue(all(r['cut_authorized'] == 'no' for r in table(self.root/'out/junction_checks.tsv')))

        # Palindromic full-window matches in both orientations are kept ambiguous.
        def ambiguous_map(args, target, query, prefix):
            if query.name == 'anchors.fa':
                return []
            return [hit(length=14, cigar='14=', start=2),
                    hit(length=14, cigar='14=', start=2, strand='-')]
        with patch.object(runner, 'Fasta', FakeFasta), patch.object(runner, 'map_queries', ambiguous_map):
            runner.process_assembly(args, 'sample_hap1', calls, stages, self.root/'ambiguous')
        self.assertEqual(table(self.root/'ambiguous/origin_ledger.tsv')[0]['status'], 'multiple_exact_placements')

        # A declared gap must actually contain Ns in its paired stage FASTA.
        assessment_path.write_text('>s\nAAAACCCCAAGGGGTTTT\n')
        calls[0]['assembly_sha256'] = sha(assessment_path)
        with patch.object(runner, 'Fasta', FakeFasta), patch.object(runner, 'map_queries', fake_map):
            with self.assertRaisesRegex(ValueError, 'AGP gap disagrees'):
                runner.process_assembly(args, 'sample_hap1', calls, stages, self.root/'bad_gap')

    def test_wrong_assessment_is_rejected_before_mapping(self):
        fasta = self.root/'a.fa'
        fasta.write_text('>s\nACGT\n')
        with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
            runner.process_assembly(SimpleNamespace(), 'sample', [dict(assembly_sha256='0'*64)],
                                    [dict(stage='assessment', fasta=str(fasta))], self.root/'out')

    def test_wrong_hifi_frame_is_rejected(self):
        context = self.root/'assembly/chimeras/sequence_context/sample.sequence_context'
        context.mkdir(parents=True)
        (context/'provenance.json').write_text(json.dumps(dict(assembly='sample', hifi=True,
            sha256='0'*64, coordinate_stage='pre_finishing', intervals={})))
        out = self.root/'out'
        out.mkdir()
        args = SimpleNamespace(hifi_results=self.root)
        with self.assertRaisesRegex(ValueError, 'different assessment'):
            runner.evidence_checks(args, 'sample', '1'*64, None, {}, {}, [], out)


if __name__ == '__main__':
    unittest.main()
