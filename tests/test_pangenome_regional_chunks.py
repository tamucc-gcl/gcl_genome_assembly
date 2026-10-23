import csv
import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
import pangenome_regional_chunks as mod
from pangenome_regional_patch import coordinate_audit


class RegionalChunksTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.root=Path(self.tmp.name)

    def test_only_chunk_option_changes(self):
        config=self.root/'config.xml';config.write_text('<config><graphmap minMAPQ="5"/><graphmap_split minQueryCoverages=".75 .5 .25" minQueryCoverageThresholds="100000 1000000" minQueryUniqueness="2"/></config>')
        original=mod.chunk_command('graph','paf',config,'ref',['A','B'],1000000,'out','log')
        for chunk in (100000,25000):
            cmd=mod.chunk_command('graph','paf',config,'ref',['A','B'],chunk,'out','log')
            diffs=[i for i,(a,b) in enumerate(zip(original,cmd)) if a!=b]
            self.assertEqual(diffs,[original.index('-u')+1]);self.assertEqual(cmd[diffs[0]],str(chunk))

    def test_unplanned_setting_fails(self):
        with self.assertRaisesRegex(ValueError,'Unsupported'):mod.chunk_command('g','p','c','r',[],10,'o','l')

    def fixtures(self):
        sources=[dict(assembly='sample',sequence='A+B',source_bp=8000000,chromosome_scale='yes',
                      status='represented_in_clip',clip_retained_path_bp=4000000),
                 dict(assembly='sample',sequence='A',source_bp=100,chromosome_scale='yes',
                      status='represented_in_clip',clip_retained_path_bp=100)]
        bases=[];comparisons=[]
        for policy in mod.POLICIES:
            baseline=self.root/('trace_'+policy);baseline.mkdir();bases.append(baseline)
            (baseline/'trace_case.json').write_text(json.dumps(dict(policy=policy)))
            checks=[dict(variant=v,byte_identical='yes',coordinate_integrity='PASS') for v in mod.VARIANTS]
            mod.write_table(baseline/'observation_equivalence.tsv',checks,list(checks[0]))
            mod.write_table(baseline/'selection.tsv',sources,list(sources[0]))
            for variant in mod.VARIANTS:
                for chunk in mod.CHUNKS:
                    if chunk==1000000:p=baseline;prefix=variant+'_'
                    else:
                        p=self.root/f'{policy}_{variant}_{chunk}';p.mkdir();comparisons.append(p);prefix=''
                        (p/'case.json').write_text(json.dumps(dict(policy=policy,variant=variant,chunk_bp=chunk,ambiguous_name='_AMBIGUOUS_')))
                        mod.write_table(p/'selection.tsv',sources,list(sources[0]))
                    spans=[('A+B','B',4000000,8000000)]
                    recovered=4000000 if chunk==25000 and variant=='fixed' else 100000 if chunk==25000 else 0
                    if recovered:spans.append(('A+B','A',0,recovered))
                    if recovered<4000000:spans.append(('A+B','_AMBIGUOUS_',recovered,4000000))
                    spans += [('A','A',0,90),('A','B',90,100)] if chunk==25000 and variant=='fixed' else [('A','A',0,100)]
                    regions=[dict(case='test',assembly='sample',sequence=s,chromosome=c,start_0based=lo,end_exclusive=hi) for s,c,lo,hi in spans]
                    assignments=[];mapped=[]
                    for source in sources:
                        s=source['sequence'];named=[r for r in regions if r['sequence']==s and r['chromosome']!='_AMBIGUOUS_']
                        bp=sum(r['end_exclusive']-r['start_0based'] for r in named)
                        assignments.append(dict(assembly='sample',sequence=s,source_bp=source['source_bp'],
                            surviving_query_covered_bp=bp,assigned_span_bp=bp,ambiguous_span_bp=source['source_bp']-bp))
                        for r in named:mapped.append(dict(assembly='sample',sequence=s,chromosome=r['chromosome'],mapped_union_bp=r['end_exclusive']-r['start_0based']))
                    mod.write_table(p/(prefix+'assignments.tsv'),assignments,list(assignments[0]))
                    mod.write_table(p/(prefix+'regions.tsv'),regions,list(regions[0]))
                    audits=coordinate_audit(regions,sources)
                    mod.write_table(p/(prefix+'coordinate_audit.tsv'),audits,list(audits[0]))
                    if chunk!=1000000:mod.write_table(p/'chromosome_mapping.tsv',mapped,list(mapped[0]))
        return bases,comparisons

    def test_complete_matrix_and_recovery_report(self):
        mod.report(*self.fixtures(),self.root/'report')
        rows=mod.read_table(self.root/'report/chunk_source_comparison.tsv')
        self.assertEqual(len(rows),36)
        r=next(r for r in rows if r['policy']=='chromosome_isolated' and r['variant']=='fixed' and r['chunk_bp']=='25000' and r['sequence']=='A+B')
        self.assertEqual(int(r['mapped_delta_bp']),4000000);self.assertEqual(r['chromosomes'],'A,B')

    def test_small_secondary_and_control_side_effects_are_flagged(self):
        mod.report(*self.fixtures(),self.root/'report')
        flags=mod.read_table(self.root/'report/chunk_review_flags.tsv')
        self.assertTrue(any(r['review_flag']=='new chromosome on ordinary control' for r in flags))
        self.assertTrue(any(r['review_flag']=='new chromosome represented only by regions below 1 Mb' for r in flags))

    def test_exact_source_interval_changes(self):
        mod.report(*self.fixtures(),self.root/'report')
        with gzip.open(self.root/'report/chunk_region_changes.tsv.gz','rt') as f:rows=list(csv.DictReader(f,delimiter='\t'))
        self.assertTrue(any(r['sequence']=='A+B' and r['variant']=='fixed' and r['chunk_bp']=='25000'
                            and r['change']=='added' and r['start_0based']=='0' and r['end_exclusive']=='4000000' for r in rows))

    def test_missing_matrix_case_fails(self):
        bases,comparisons=self.fixtures()
        with self.assertRaisesRegex(ValueError,'Incomplete 18-case'):mod.load_cases(bases,comparisons[:-1])

    def test_duplicate_matrix_case_fails(self):
        bases,comparisons=self.fixtures()
        with self.assertRaisesRegex(ValueError,'duplicate'):mod.load_cases(bases,comparisons+[comparisons[0]])

    def test_unverified_cached_baseline_fails(self):
        bases,comparisons=self.fixtures();p=bases[0]/'observation_equivalence.tsv';p.write_text(p.read_text().replace('yes','no'))
        with self.assertRaisesRegex(ValueError,'Unverified'):mod.load_cases(bases,comparisons)

    def test_overlapping_assignments_fail(self):
        bases,comparisons=self.fixtures();p=comparisons[0]/'regions.tsv'
        rows=mod.read_table(p);rows.append(rows[0].copy());mod.write_table(p,rows,list(rows[0]))
        with self.assertRaisesRegex(ValueError,'Overlapping'):mod.load_cases(bases,comparisons)

    def test_missing_coordinate_audit_fails(self):
        bases,comparisons=self.fixtures();p=comparisons[0]/'coordinate_audit.tsv';rows=mod.read_table(p)
        mod.write_table(p,rows[:1],list(rows[0]))
        with self.assertRaisesRegex(ValueError,'Incomplete coordinate'):mod.load_cases(bases,comparisons)

    def test_changed_source_selection_fails(self):
        bases,comparisons=self.fixtures();p=comparisons[0]/'selection.tsv';p.write_text(p.read_text().replace('8000000','9000000'))
        with self.assertRaisesRegex(ValueError,'source selection'):mod.load_cases(bases,comparisons)

    def test_loss_of_mapping_is_explicit(self):
        bases,comparisons=self.fixtures();p=comparisons[0]/'assignments.tsv';rows=mod.read_table(p);rows[1]['surviving_query_covered_bp']='50'
        mod.write_table(p,rows,list(rows[0]));mod.report(bases,comparisons,self.root/'report')
        summaries=mod.read_table(self.root/'report/chunk_matrix_summary.tsv')
        self.assertTrue(any(r['controls_losing']=='1' for r in summaries))


if __name__=='__main__':unittest.main()
