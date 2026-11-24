import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'py_scripts'))
from contamination_cohort_report import build, rows


class CohortReportTests(unittest.TestCase):
    def test_both_sources_keep_frames_and_all_decision_fields(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            action = root / 'A.fcs_gx_report.txt'
            action.write_text('## metadata\n#seq_id\tstart_pos\tend_pos\taction\nc1\t1\t100\tEXCLUDE\n')
            review = root / 'review'
            review.mkdir()
            (review / 'evidence.json').write_text(json.dumps(dict(assembly='A', whole_contig_candidates=0,
                protein_searched_bp=100, megablast_searched_bp=50, nucleotide_status='not-run',
                search_mode='staged', window_bp=5000, minimum_acgt_bp=100,
                sampling=dict(spacing_bp=500000, short_contig_bp=50000, gc_delta=.12, depth_fold=3))))
            (review / 'review.tsv').write_text('action_id\tsource_sha256\tstage\tdecision\tsupport\nA1\tchecksum\tpost_fcs\tkeep\tyes\nA2\tchecksum\tpost_fcs\tkeep\t\n')
            (review / 'sequence_evidence.tsv').write_text('sequence\tlength_bp\nc1\t1000\n')
            spec = dict(stage='contig', assemblies=['A'], fcs_enabled=True, blobtools_enabled=True,
                fcs=[dict(assembly='A', actions=str(action))], blobtools=[dict(assembly='A', review=str(review))])
            metrics = root/'A.fcs_metrics.tsv'
            metrics.write_text('sequence\tlength_bp\tgc_fraction\tmean_depth\tplot_taxon\tclassification_source\nc1\t1000\t0.4\t20\tanml:fishes\tFCS-GX\n')
            spec['fcs_metrics']=[dict(assembly='A',metrics=str(metrics))]
            build(spec, root / 'out')
            summary = rows(root / 'out/report/cohort_summary.tsv')[0]
            self.assertEqual(summary['evidence_intervals'], '1')
            self.assertEqual(summary['whole_contig_candidates'], '0')
            decisions = rows(root / 'out/report/cohort_review.tsv')
            self.assertEqual(len(decisions), 2)
            self.assertEqual(decisions[0]['source_sha256'], 'checksum')
            text = (root / 'out/contamination_report.md').read_text()
            self.assertIn('(fcs/A.fcs_gx_report.txt)', text)
            self.assertIn('1-based inclusive', text)
            self.assertIn('0-based half-open', text)
            self.assertIn('Higher-priority intervals | All signal intervals', text)
            self.assertIn('Protein BLASTX searched: Mb (%)', text)
            self.assertIn('10.00%', text)
            self.assertIn('cohort_gc_coverage.png', text)
            self.assertNotIn('Priority / all evidence', text)
            self.assertNotIn('500,000-bp block', text)
            self.assertIn('FCS-GX divisions', text)
            self.assertEqual(rows(root/'out/report/cohort_fcs_sequence_evidence.tsv')[0]['classification_source'],'FCS-GX')

    def test_disabled_is_not_a_negative_result(self):
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            build(dict(stage='contig', assemblies=['A'], fcs_enabled=False, blobtools_enabled=False), out)
            row = rows(out / 'report/cohort_summary.tsv')[0]
            self.assertEqual(row['fcs_status'], 'not run')
            self.assertEqual(row['whole_contig_candidates'], '')
            self.assertNotIn('![GC and read coverage', (out/'contamination_report.md').read_text())

    def test_enabled_missing_fails(self):
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaisesRegex(ValueError, 'membership differs'):
                build(dict(stage='contig', assemblies=['A'], fcs_enabled=True, blobtools_enabled=False), Path(temp))

    def test_final_adaptor_finding_requires_cleanup_without_claiming_clean(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            gx=root/'A.fcs_gx_report.txt'; gx.write_text('#seq_id\tstart_pos\tend_pos\taction\n')
            adapter=root/'A.fcs_adaptor_report.txt'
            adapter.write_text('#accession\tlength\taction\trange\tname\nc1\t1000\tACTION_TRIM\t1..58\tadaptor\n')
            build(dict(stage='final',assemblies=['A'],fcs_enabled=True,blobtools_enabled=False,
                fcs=[dict(assembly='A',actions=str(gx))],adaptor_reports=[dict(assembly='A',report=str(adapter))]),root/'out')
            row=rows(root/'out/report/cohort_summary.tsv')[0]
            self.assertEqual(row['audit_status'],'REVIEW / CLEANUP REQUIRED')
            self.assertEqual(row['adaptor_actions'],'1')
            self.assertIn('No sequence is edited by this audit',(root/'out/contamination_report.md').read_text())


if __name__ == '__main__':
    unittest.main()
