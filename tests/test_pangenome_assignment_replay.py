import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import Mock, patch
from types import SimpleNamespace

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from pangenome_assignment_replay import (CASES, initialize_cactus_tools, identity_aliases,
    select_queries, split_command, fragment, summarize_split, parse_decisions, report,
    merged, interval_difference, bin_coverage, mapping_intervals, write_intervals, coverage_report)
from pangenome_source_audit import write_table


CONFIG = '''<config><graphmap minMAPQ="5"/><graphmap_split minQueryCoverages="0.75 0.5 0.25"
minQueryCoverageThresholds="100000 1000000" minQueryUniqueness="2" ambiguousName="_AMBIGUOUS_"/></config>'''


class AssignmentReplayTests(unittest.TestCase):
    def test_imported_cactus_helpers_explicitly_use_existing_container_tools(self):
        setup=Mock()
        with patch.dict(sys.modules,{'cactus.shared.common':SimpleNamespace(setupBinaries=setup)}):
            initialize_cactus_tools()
        options=setup.call_args.args[0]
        self.assertEqual(options.binariesMode,'local')
        self.assertFalse(options.latest)
        self.assertIsNone(options.workDir)

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.config=self.root/'config.xml';self.config.write_text(CONFIG)
        self.identities=self.root/'identities.tsv'
        self.identities.write_text('id\tgraph_name\nisland\tIND_island.1\n')

    def row(self,seq='chr1',status='represented_in_clip',large='yes'):
        return dict(assembly='island',sequence=seq,source_bp='100',chromosome_scale=large,
                    status=status,filtered_query_covered_bp='90')

    def command(self,case):
        return split_command('graph','paf',self.config,'REF',['chr1','chr2'],case,'split_','log')

    def test_baseline_matches_pinned_coverage_and_uniqueness_options(self):
        cmd=self.command('baseline')
        self.assertEqual(cmd[cmd.index('-Q')+1],'2')
        self.assertEqual([cmd[i+1] for i,x in enumerate(cmd) if x=='-n'],['0.75','0.5','0.25'])
        self.assertEqual(cmd[cmd.index('-r')+1],'id=REF|')
        self.assertNotIn('-G',cmd)

    def test_relaxed_control_changes_only_uniqueness(self):
        expected=self.command('baseline');expected[expected.index('-Q')+1]='1'
        self.assertEqual(self.command('whole_relaxed'),expected)

    def test_regional_cases_keep_uniqueness_and_use_generic_chunk_sizes(self):
        for case,size in [('regional_1Mb','1000000'),('regional_3Mb','3000000')]:
            self.assertEqual(self.command(case),self.command('baseline')+['-u',size,'-s'])

    def test_unknown_cases_fail(self):
        with self.assertRaisesRegex(ValueError,'Unknown assignment'):self.command('custom')

    def test_missing_and_composite_paths_selected_without_sample_specific_rules(self):
        rows=[self.row('unplaced','missing_path_mapping_unavailable','no'),self.row('chr8+chr10','represented_in_clip','no')]
        self.assertEqual(select_queries(rows),rows)

    def test_bounded_controls_are_selected_per_assembly(self):
        rows=[self.row('chr'+str(i)) for i in range(10)]
        self.assertEqual(len(select_queries(rows)),2)

    def test_low_retention_paths_are_selected_beyond_normal_control_limit(self):
        rows=[self.row('chr'+str(i)) for i in range(5)]
        rows[-1]['clip_retained_path_bp']='49'
        self.assertIn(rows[-1],select_queries(rows))

    def test_coverage_uses_interval_sets_and_distinguishes_loss_from_gain(self):
        left=merged([(0,10),(5,20),(20,30),(50,70)])
        right=merged([(10,40),(60,80)])
        self.assertEqual(left,[(0,30),(50,70)])
        self.assertEqual(interval_difference(left,right),[(0,10),(50,60)])
        self.assertEqual(interval_difference(right,left),[(30,40),(70,80)])
        self.assertEqual(bin_coverage([(900000,1100000),(2000000,2000003)],2000003),
                         [100000,100000,3])

    def test_original_mapping_reader_unions_overlap_and_ignores_unselected_queries(self):
        paf=self.root/'raw.paf'
        paf.write_text('id=IND_island.1|chr1\t100\t0\t50\t+\ts\t100\t0\t50\t50\t50\t60\n'
                       'id=IND_island.1|chr1\t100\t20\t80\t+\ts\t100\t0\t60\t60\t60\t255\n'
                       'id=_MINIGRAPH_|s\t100\t0\t90\t+\ts\t100\t0\t90\t90\t90\t60\n')
        spans,counts=mapping_intervals(paf,{'IND_island.1':'island'},
                                      {('island','chr1'):100,('island','chr2'):80})
        self.assertEqual(spans[('island','chr1')],[(0,80)])
        self.assertEqual(spans[('island','chr2')],[])
        self.assertEqual(counts[('island','chr1')],2)
        paf.write_text(paf.read_text().replace('\t100\t0\t50','\t101\t0\t50',1))
        with self.assertRaisesRegex(ValueError,'coordinates disagree'):
            mapping_intervals(paf,{'IND_island.1':'island'},{('island','chr1'):100})

    def test_full_report_includes_present_composites_and_exact_stage_losses(self):
        source=self.row('chr8+chr10')
        source.update(clip_retained_path_bp='45',retention_review='yes')
        selection=self.root/'selection.tsv';write_table(selection,[source],list(source))
        key=('island','chr8+chr10')
        coverage=self.root/'coverage.tsv'
        write_intervals(coverage,{'raw_gaf':{key:[(0,100)]},'raw_paf':{key:[(0,90)]},
                                  'filtered_paf':{key:[(50,90)]}})
        dirs=[]
        for case in CASES:
            directory=self.root/case;directory.mkdir();dirs.append(directory)
            row=dict(case=case,assembly=key[0],sequence=key[1],source_bp='100',
                chromosome_scale='yes',graph_status='represented_in_clip',
                filtered_query_covered_bp='40',assigned_chromosomes='chr10',
                assigned_span_bp='100',ambiguous_span_bp='0',surviving_query_covered_bp='30',
                retention_review='yes')
            write_table(directory/'assignments.tsv',[row],list(row))
            write_table(directory/'decisions.tsv',[],['assembly','sequence'])
            write_intervals(directory/'surviving_intervals.tsv.gz',{case:{key:[(60,90)]}})
            write_table(directory/'regions.tsv',[dict(case=case,assembly=key[0],sequence=key[1],
                chromosome='chr10',start_0based=0,end_exclusive=100)],
                ['case','assembly','sequence','chromosome','start_0based','end_exclusive'])
        out=self.root/'summary'
        report(dirs,out,coverage,selection)
        text=(out/'assignment_report.md').read_text(encoding='utf-8')
        self.assertIn('island / chr8+chr10',text)
        self.assertIn('coverage_changes.tsv.gz',text)
        from pangenome_source_audit import read_table
        trace={r['stage']:r for r in read_table(out/'coverage_trace.tsv')}
        self.assertEqual(trace['raw_paf']['lost_covered_bp'],'10')
        self.assertEqual(trace['filtered_paf']['lost_covered_bp'],'50')
        self.assertEqual(trace['baseline']['lost_covered_bp'],'10')
        self.assertEqual(trace['baseline']['covered_bp'],'30')
        changes=read_table(out/'coverage_changes.tsv.gz')
        self.assertTrue(any(r['stage']=='filtered_paf' and r['start_0based']=='0'
                            and r['end_exclusive']=='50' for r in changes))
        self.assertIn('| regional_1Mb | 0 | 0 |',text)  # composites are not normal controls

    def test_fragment_coordinates_are_one_based_in_tool_zero_based_in_report(self):
        aliases={'IND_island.1':'island'};lengths={('island','chr1'):100}
        self.assertEqual(fragment('id=IND_island.1|chr1:11-60',aliases,lengths),(('island','chr1'),10,60))
        self.assertEqual(fragment('id=IND_island.1|chr1',aliases,lengths),(('island','chr1'),0,100))
        self.assertIsNone(fragment('id=_MINIGRAPH_|s12',aliases,lengths))

    def test_invalid_fragment_cannot_be_silently_reported(self):
        with self.assertRaisesRegex(ValueError,'Invalid replay fragment'):
            fragment('id=IND_island.1|chr1:0-60',{'IND_island.1':'island'},{('island','chr1'):100})

    def test_surviving_alignments_are_not_equated_to_assigned_whole_scaffold(self):
        row=self.row(status='mapped_but_absent_from_full_and_clip')
        selection=self.root/'selection.tsv';write_table(selection,[row],list(row))
        (self.root/'split_chr1.fa_contigs').write_text('id=IND_island.1|chr1\nid=_MINIGRAPH_|s1\n')
        (self.root/'split_chr1.paf').write_text('id=IND_island.1|chr1\t100\t0\t40\t+\ts1\t40\t0\t40\t40\t40\t60\n')
        rows,regions=summarize_split(self.root,selection,self.identities,'whole_relaxed','_AMBIGUOUS_')
        self.assertEqual(rows[0]['assigned_span_bp'],100)
        self.assertEqual(rows[0]['surviving_query_covered_bp'],40)

    def test_chunk_alignment_offsets_reconcile_to_source_coordinates(self):
        row=self.row(status='mapped_but_absent_from_full_and_clip')
        selection=self.root/'selection.tsv';write_table(selection,[row],list(row))
        (self.root/'split_chr1.fa_contigs').write_text('id=IND_island.1|chr1:11-60\n')
        (self.root/'split_chr1.paf').write_text('id=IND_island.1|chr1:11-60\t50\t0\t50\t+\ts1\t50\t0\t50\t50\t50\t60\n')
        rows,regions=summarize_split(self.root,selection,self.identities,'regional_1Mb','_AMBIGUOUS_')
        self.assertEqual(regions[0]['start_0based'],10)
        self.assertEqual(regions[0]['end_exclusive'],60)
        self.assertEqual(rows[0]['surviving_query_covered_bp'],50)

    def test_ambiguous_span_is_not_counted_as_surviving_query(self):
        row=self.row(status='mapped_but_absent_from_full_and_clip')
        selection=self.root/'selection.tsv';write_table(selection,[row],list(row))
        (self.root/'split__AMBIGUOUS_.fa_contigs').write_text('id=IND_island.1|chr1\n')
        (self.root/'split__AMBIGUOUS_.paf').write_text('id=IND_island.1|chr1\t100\t0\t100\t+\ts1\t100\t0\t100\t100\t100\t60\n')
        rows,_=summarize_split(self.root,selection,self.identities,'baseline','_AMBIGUOUS_')
        self.assertEqual(rows[0]['ambiguous_span_bp'],100)
        self.assertEqual(rows[0]['surviving_query_covered_bp'],0)

    def test_report_requires_all_four_cases(self):
        directory=self.root/'baseline';directory.mkdir()
        row=dict(case='baseline',assembly='a',sequence='chr1')
        write_table(directory/'assignments.tsv',[row],list(row))
        with self.assertRaisesRegex(ValueError,'Incomplete assignment'):report([directory],self.root/'summary')

    def test_native_log_identifies_uniqueness_rejection_separately_from_coverage(self):
        path=self.root/'assignment.txt'
        path.write_text('Query contig is ambiguous: id=IND_island.1|chr1  len=100 cov=0.51 (vs 0.25) uf=1.1 (vs 2)\n')
        row=parse_decisions(path,{'IND_island.1':'island'})[0]
        self.assertEqual(row['reason'],'chromosome assignment not unique enough')

    def test_native_log_parses_infinite_uniqueness_and_assigned_reference(self):
        path=self.root/'assignment.txt'
        path.write_text('Assigned ref-contig to chr1: id=IND_island.1|chr1  len=100 cov=1 (vs 0.25) uf= infinity (vs 2)\n')
        row=parse_decisions(path,{'IND_island.1':'island'})[0]
        self.assertEqual(row['decision'],'assigned')
        self.assertEqual(row['uniqueness_ratio'],'infinity')


if __name__=='__main__':unittest.main()
