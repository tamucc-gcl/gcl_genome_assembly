import gzip
import json
import shutil
import subprocess
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
import pangenome_cactus_policy as mod
from pangenome_sequence_audit import scan_graph, resolve_sequence, audit as graph_audit


class CactusPolicyTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.root=Path(self.tmp.name)
        self.sources=[dict(query='id=S.1|a',sample='S.1',sequence='a',source_bp=100)]
        mod.write_table(self.root/'lengths.tsv',self.sources,list(self.sources[0]))
        self.data=dict(source_lengths=str(self.root/'lengths.tsv'),ambiguous_name='_AMBIGUOUS_',evidence=str(self.root))

    def test_existing_chunk_flag_is_replaced_once(self):
        cmd=mod.regional_command(dict(corrected_binary='fixed',chunk_bp=25000),['-g','g','-u','1000000','-s','-p','p'])
        self.assertEqual(cmd,['fixed','-g','g','-p','p','-u','25000','-s'])

    def test_fragment_coordinates(self):
        self.assertEqual(mod.source_fragment('id=S.1|a:11-50',{'id=S.1|a':100}),('id=S.1|a',10,50))
        with self.assertRaisesRegex(ValueError,'outside'):mod.source_fragment('id=S.1|a:11-101',{'id=S.1|a':100})

    def test_unknown_biological_source_fails(self):
        with self.assertRaisesRegex(ValueError,'Unknown'):mod.source_fragment('id=other|a',{})

    def test_literal_source_identifier_takes_precedence(self):
        self.assertEqual(mod.source_fragment('id=S.1|a:1-5',{'id=S.1|a:1-5':100}),('id=S.1|a:1-5',0,100))

    def test_overlap_is_rejected_and_gaps_explicit(self):
        self.assertEqual(mod.complement([(10,30),(40,90)],100),[(0,10),(30,40),(90,100)])
        with self.assertRaisesRegex(ValueError,'Overlapping'):mod.complement([(0,60),(50,100)],100)

    def test_split_audit_reports_exclusions(self):
        (self.root/'split_chrA.fa_contigs').write_text('id=S.1|a:11-80\n')
        (self.root/'split__AMBIGUOUS_.fa_contigs').write_text('id=S.1|a:81-90\n')
        paf=self.root/'input.paf';paf.write_text('')
        data=mod.audit_split(self.data,self.root/'split_',paf,'split_test')
        self.assertEqual((data['assigned_bp'],data['ambiguous_bp'],data['unreported_bp']),(70,10,20))
        with gzip.open(self.root/'split_test.excluded.tsv.gz','rt') as f:text=f.read()
        self.assertIn('outside_partition',text);self.assertIn('ambiguous_not_in_chromfile',text)

    def test_mapping_outside_assigned_fragment_fails(self):
        (self.root/'split_chrA.fa_contigs').write_text('id=S.1|a:11-80\n')
        (self.root/'split_chrA.paf').write_text('id=S.1|a:11-80\t70\t0\t71\t+\tn\t100\t0\t71\t71\t71\t60\n')
        paf=self.root/'in.paf';paf.write_text('')
        with self.assertRaisesRegex(ValueError,'outside'):mod.audit_split(self.data,self.root/'split_',paf,'split_test')

    def test_fragment_paf_name_matches_sanitized_fasta_without_changing_coordinates(self):
        name='id=S.1|a:11-80'
        regions=self.root/'split_A.fa_contigs';regions.write_text(name+'\n')
        paf=self.root/'split_A.paf'
        original=f'{name}\t70\t5\t65\t+\tnode\t100\t9\t69\t60\t60\t60\tcg:Z:60M\n'
        paf.write_text(original)
        result=mod.normalize_split_pafs(self.data,self.root/'split_')
        self.assertEqual(paf.read_text(),original.replace(name,'id=S.1|a_sub_10_80',1))
        self.assertEqual(regions.read_text(),name+'\n')
        self.assertEqual(result['renamed_mapping_records'],1)
        self.assertEqual(resolve_sequence('S.1','a_sub_10_80',{'id=S.1|a':100}),('id=S.1|a',10,80))

    def test_fragment_alias_collision_fails(self):
        sources=self.sources+[dict(query='id=S.1|a_sub_10_80',sample='S.1',sequence='a_sub_10_80',source_bp=70)]
        mod.write_table(self.root/'lengths.tsv',sources,list(sources[0]))
        (self.root/'split_A.fa_contigs').write_text('id=S.1|a:11-80\n')
        with self.assertRaisesRegex(ValueError,'collides'):mod.normalize_split_pafs(self.data,self.root/'split_')

    @unittest.skipUnless(shutil.which('sed'),'GNU sed regression runs in Linux integration checks')
    def test_pinned_upstream_rename_leaves_colon_paf_unchanged(self):
        sys.path.insert(0,str(Path(__file__).parent/'fixtures'))
        from cactus_fragment_rename import get_faidx_subpath_rename_cmd
        raw='id=S.1|a:11-80\t70\t0\t70\t+\tnode\t100\t0\t70\t70\t70\t60\n'
        result=subprocess.run(get_faidx_subpath_rename_cmd(),input=raw,text=True,stdout=subprocess.PIPE,check=True)
        self.assertEqual(result.stdout,raw)
        fixed=raw.replace('a:11-80','a_sub_10_80')
        result=subprocess.run(get_faidx_subpath_rename_cmd()+['-e','s/ /\t/g'],input=fixed,text=True,stdout=subprocess.PIPE,check=True)
        self.assertEqual(result.stdout,fixed)

    def graph(self,walks,segments='S\t1\tACGT\n'):
        p=self.root/'graph.gfa';p.write_text(segments+walks);return p

    def test_graph_nodes_counted_once_across_multiple_walks(self):
        sources=self.sources+[dict(query='id=T|a',sample='T',sequence='a',source_bp=100)]
        p=self.graph('W\tS\t1\ta\t0\t4\t>1\nW\tT\t0\ta\t0\t4\t>1\n')
        stats,rows=scan_graph(p,sources,'FULL',self.root)
        self.assertEqual(stats['distinct_node_bp'],4);self.assertEqual(stats['biological_source_span_bp'],8)

    def test_duplicate_nodes_fail(self):
        p=self.graph('W\tS\t1\ta\t0\t4\t>1\n','S\t1\tACGT\nS\t1\tACGT\n')
        with self.assertRaisesRegex(ValueError,'Duplicate graph'):scan_graph(p,self.sources,'FULL',self.root)

    def test_overlapping_walks_fail(self):
        p=self.graph('W\tS\t1\ta\t0\t70\t>1\nW\tS\t1\ta#1\t60\t100\t>1\n')
        with self.assertRaisesRegex(ValueError,'Overlapping source'):scan_graph(p,self.sources,'FULL',self.root)

    def test_source_coordinate_bounds_fail(self):
        p=self.graph('W\tS\t1\ta\t0\t101\t>1\n')
        with self.assertRaisesRegex(ValueError,'exceeds'):scan_graph(p,self.sources,'FULL',self.root)

    def test_fragment_walk_offsets(self):
        p=self.graph('W\tS\t1\ta[40-80]\t0\t40\t>1\n')
        stats,rows=scan_graph(p,self.sources,'FULL',self.root);self.assertEqual(rows[0]['source_span_bp'],40)
        self.assertEqual(resolve_sequence('S.1','a_sub_40_80',{'id=S.1|a':100}),('id=S.1|a',40,80))

    def test_unknown_walk_sample_fails(self):
        p=self.graph('W\tunknown\t0\ta\t0\t4\t>1\n')
        with self.assertRaisesRegex(ValueError,'Unrecognized'):scan_graph(p,self.sources,'FULL',self.root)

    def test_source_stage_differences_are_explicit(self):
        policy=self.root/'policy';policy.mkdir()
        mod.write_table(policy/'source_lengths.tsv',self.sources,list(self.sources[0]))
        mod.dump(policy/'sequence_accounting.json',dict(policy=mod.POLICY))
        mod.dump(policy/'policy.json',dict(ambiguous_name='_AMBIGUOUS_'))
        row=dict(**self.sources[0],assigned_bp=30,ambiguous_bp=0,unreported_bp=70)
        mod.write_table(policy/'split_test.sources.tsv',[row],list(row))
        with gzip.open(policy/'split_test.regions.tsv.gz','wt') as f:
            f.write('query\tchromosome\tstart_0based\tend_exclusive\nid=S.1|a\tA\t0\t30\n')
        full=self.root/'full.gfa';full.write_text('S\t1\tACGT\nW\tS\t1\ta\t20\t50\t>1\n')
        clip=self.root/'clip.gfa';clip.write_text('S\t1\tACGT\nW\tS\t1\ta\t10\t40\t>1\n')
        graph_audit(full,clip,policy,self.root/'out')
        report=json.loads((self.root/'out/sequence_graph_audit.json').read_text())
        self.assertEqual(report['representation_review'],'required')
        self.assertEqual(report['coordinate_discrepancies'],dict(assigned_missing_from_FULL_bp=20,FULL_outside_assigned_bp=20,CLIP_outside_FULL_bp=10))

    def test_no_hidden_gap_in_source_accounting(self):
        stats=mod.audit_split(self.data,self.root/'split_',self._empty_paf(),'split_none')
        self.assertEqual(stats['assigned_bp'],0);self.assertEqual(stats['unreported_bp'],100)

    def _empty_paf(self):
        p=self.root/'empty.paf';p.write_text('');return p


if __name__=='__main__':unittest.main()
