import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
import pangenome_regional_trace as mod
from pangenome_regional_patch import verified_source, OLD, NEW


class RegionalTraceTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.root=Path(self.tmp.name)
        self.identities=self.root/'ids.tsv';self.identities.write_text('id\tgraph_name\nsample\tIND_sample.1\n')
        self.sources=[dict(assembly='sample',sequence='piece',source_bp=8000000)]

    def test_instrumentation_preserves_every_original_line(self):
        _,text=verified_source(Path(__file__).resolve().parents[1]/'vendor/cactus-gfa-tools-regional')
        for original in (text,text.replace(OLD,NEW)):
            observed=mod.instrument(original)
            # Every original line must appear in order, with additions only.
            it=iter(observed.splitlines())
            for line in original.splitlines():
                self.assertTrue(any(x==line for x in it),line)
            for stage in mod.FIELDS:self.assertIn('TRACE\\t'+stage,observed)

    def test_context_mismatch_fails(self):
        with self.assertRaisesRegex(ValueError,'exactly once'):mod.instrument('wrong source')

    def test_byte_equivalence_is_required(self):
        a=self.root/'a';b=self.root/'b';a.mkdir();b.mkdir()
        for p in (a,b):(p/'split_A.paf').write_text('same\n')
        self.assertEqual(mod.verify_equivalent(a,b),1)
        (b/'split_A.paf').write_text('different\n')
        with self.assertRaisesRegex(ValueError,'output bytes'):mod.verify_equivalent(a,b)

    def test_missing_output_file_fails(self):
        a=self.root/'a';b=self.root/'b';a.mkdir();b.mkdir()
        (a/'split_A.paf').write_text('same')
        with self.assertRaisesRegex(ValueError,'output bytes'):mod.verify_equivalent(a,b)

    def event(self,stage,*values):
        return '\t'.join(('TRACE',stage,'id=IND_sample.1|piece',*map(str,values)))+'\n'

    def parse(self,lines):
        p=self.root/'log';p.write_text('native human log\n'+''.join(lines))
        return mod.parse_trace(p,self.identities,self.sources,self.root/'events.tsv.gz')

    def setup_events(self):
        return [self.event('before_overlap','A',3,900000,2500000,2400000),
                self.event('after_overlap','A',3,900000,2500000,2400000),
                self.event('before_overlap','B',1,5000000,5000000,4800000),
                self.event('after_overlap','B',1,5000000,5000000,4800000),
                self.event('dominant','B',4800000,1000000)]

    def final(self,passed=1):
        return self.event('final','B',4800000,8000000,0,0,0,.6,.25,passed)

    def test_long_run_with_short_intervals_is_length_rejection(self):
        rows=self.parse(self.setup_events()+[
            self.event('candidate','A',0,3000000,3,900000,1000000,2400000,.8,.25,0,1),self.final()])
        a=next(r for r in rows if r['chromosome']=='A')
        self.assertIn('longest interval too short',a['reason'])
        self.assertEqual(a['largest_candidate_span_bp'],3000000)
        self.assertEqual(a['length_failed_runs'],1);self.assertEqual(a['density_failed_runs'],0)
        with gzip.open(self.root/'events.tsv.gz','rt') as f:self.assertIn('candidate',f.read())

    def test_density_is_separate_from_length(self):
        rows=self.parse(self.setup_events()+[
            self.event('candidate','A',0,3000000,3,2000000,1000000,300000,.1,.25,1,0),self.final()])
        a=next(r for r in rows if r['chromosome']=='A')
        self.assertIn('density too low',a['reason']);self.assertEqual(a['length_failed_runs'],0)

    def test_failure_counts_can_overlap(self):
        rows=self.parse(self.setup_events()+[
            self.event('candidate','A',0,3000000,3,900000,1000000,300000,.1,.25,0,0),self.final()])
        a=next(r for r in rows if r['chromosome']=='A')
        self.assertEqual((a['candidate_runs'],a['length_failed_runs'],a['density_failed_runs']),(1,1,1))

    def test_local_pass_can_be_overruled_by_whole_query_gate(self):
        rows=self.parse(self.setup_events()+[
            self.event('candidate','A',0,3000000,3,2000000,1000000,2400000,.8,.25,1,1),self.final(0)])
        self.assertTrue(all(r['reason']=='Whole-query adjusted coverage failed' for r in rows))
        self.assertEqual(next(r for r in rows if r['chromosome']=='A')['locally_accepted_runs'],1)

    def test_overlap_removal_is_distinguishable(self):
        rows=self.parse([self.event('before_overlap','A',1,500000,500000,400000)])
        self.assertEqual(rows[0]['reason'],'All intervals removed by overlap resolution')

    def test_missing_final_trace_fails(self):
        with self.assertRaisesRegex(ValueError,'Incomplete'):self.parse(self.setup_events())

    def test_unknown_query_fails(self):
        with self.assertRaisesRegex(ValueError,'Unexpected'):
            self.parse([self.event('before_overlap','A',1,10,10,10).replace('IND_sample.1','UNKNOWN')])

    def test_malformed_trace_fails(self):
        with self.assertRaisesRegex(ValueError,'Malformed'):self.parse(['TRACE\tcandidate\tbad\n'])

    def report_inputs(self):
        inputs=[]
        base=self.parse(self.setup_events()+[
            self.event('candidate','A',0,3000000,3,900000,1000000,2400000,.8,.25,0,1),self.final()])
        for i,policy in enumerate(('native','chromosome_isolated','disabled')):
            p=self.root/policy;p.mkdir();inputs.append(p)
            (p/'trace_case.json').write_text(json.dumps(dict(policy=policy,version=mod.TRACE_VERSION)))
            checks=[dict(variant=v,byte_identical='yes',coordinate_integrity='PASS') for v in ('native','fixed')]
            mod.write_table(p/'observation_equivalence.tsv',checks,list(checks[0]))
            rows=[dict(**r,variant=v,final_assigned_span_bp=0,outcome='Not assigned') for v in ('native','fixed') for r in base]
            mod.write_table(p/'chromosome_gates.tsv',rows,list(rows[0]))
            source=dict(**self.sources[0],chromosome_scale='yes',clip_retained_path_bp=100,filtered_query_covered_bp=100+i)
            mod.write_table(p/'selection.tsv',[source],list(source))
        return inputs

    def test_report_allows_expected_filter_coverage_differences(self):
        mod.report(self.report_inputs(),self.root/'report')
        self.assertIn('longest interval too short',(self.root/'report/regional_gate_report.md').read_text())

    def test_report_refuses_changed_source_identity(self):
        inputs=self.report_inputs();p=inputs[0]/'selection.tsv';p.write_text(p.read_text().replace('piece','another'))
        with self.assertRaisesRegex(ValueError,'source selections'):mod.report(inputs,self.root/'report')

    def test_report_refuses_missing_policies(self):
        with self.assertRaisesRegex(ValueError,'Incomplete trace policies'):mod.report(self.report_inputs()[:2],self.root/'report')

    def test_report_refuses_changed_observation_outputs(self):
        inputs=self.report_inputs();p=inputs[0]/'observation_equivalence.tsv';p.write_text(p.read_text().replace('yes','no'))
        with self.assertRaisesRegex(ValueError,'validation failed'):mod.report(inputs,self.root/'report')


if __name__=='__main__':unittest.main()
