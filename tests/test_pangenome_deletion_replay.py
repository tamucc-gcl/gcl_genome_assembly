import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from pangenome_source_audit import read_table, write_table
from pangenome_deletion_replay import node_owners, lane_paf, restore_names, deletion_command, report, run_filter, POLICIES, ASSIGNMENTS


class DeletionReplayTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.gfa=self.root/'graph.gfa'
        self.gfa.write_text('S\ta\tACGT\tSN:Z:chrA\tSO:i:0\tSR:i:0\n'
            'S\tb\tACGT\tSN:Z:chrB\tSO:i:0\tSR:i:0\n'
            'S\tx\tAC\tSN:Z:other\tSO:i:0\tSR:i:1\n'
            'S\ty\tAC\tSN:Z:other\tSO:i:2\tSR:i:1\n'
            'L\ta\t+\tb\t+\t0M\nL\ta\t+\tx\t+\t0M\nL\tx\t+\ty\t+\t0M\n')
        self.owners=self.root/'owners.tsv.gz'

    def test_reference_nodes_are_boundaries_not_chromosome_mergers(self):
        node_owners(self.gfa,self.owners)
        rows={r['node']:r for r in read_table(self.owners)}
        self.assertEqual(rows['a']['reference_chromosome'],'chrA')
        self.assertEqual(rows['b']['reference_chromosome'],'chrB')
        self.assertEqual(rows['y']['reference_chromosome'],'chrA')

    def test_conflicting_nonreference_component_is_not_guessed(self):
        with self.gfa.open('a') as f:f.write('L\ty\t+\tb\t+\t0M\n')
        node_owners(self.gfa,self.owners)
        rows={r['node']:r for r in read_table(self.owners)}
        self.assertEqual(rows['x']['status'],'conflicting')
        paf=self.root/'in.paf';paf.write_text(self.paf('x'))
        with self.assertRaisesRegex(ValueError,'Cannot safely assign'):
            lane_paf(paf,self.owners,self.root/'out.paf',self.root/'lanes.tsv')

    def test_missing_reference_tags_fail_early(self):
        self.gfa.write_text('S\tx\tACGT\n')
        with self.assertRaisesRegex(ValueError,'rank-zero'):node_owners(self.gfa,self.owners)

    def test_unknown_and_unanchored_mapping_nodes_fail(self):
        with self.gfa.open('a') as f:f.write('S\tisolated\tAC\n')
        node_owners(self.gfa,self.owners)
        for target in ('unknown','isolated'):
            paf=self.root/'in.paf';paf.write_text(self.paf(target))
            with self.assertRaisesRegex(ValueError,'Cannot safely assign'):
                lane_paf(paf,self.owners,self.root/'out.paf',self.root/'lanes.tsv')

    @staticmethod
    def paf(node,start=0,end=2):
        return f'id=sample|chr8+chr10\t100\t{start}\t{end}\t+\tid=_MINIGRAPH_|{node}\t4\t0\t2\t2\t2\t60\n'

    def test_lanes_preserve_coordinates_lengths_order_and_restore_every_field(self):
        node_owners(self.gfa,self.owners)
        source=self.root/'in.paf';original=self.paf('a')+self.paf('b',50,52)+self.paf('x',10,12)
        source.write_text(original)
        lanes=self.root/'lanes.paf';ledger=self.root/'lanes.tsv'
        names=lane_paf(source,self.owners,lanes,ledger)
        rows=[line.split('\t') for line in lanes.read_text().splitlines()]
        self.assertEqual(rows[0][0],rows[2][0]);self.assertNotEqual(rows[0][0],rows[1][0])
        self.assertEqual(rows[1][1:4],['100','50','52'])
        restored=self.root/'restored.paf';restore_names(lanes,restored,names)
        self.assertEqual(restored.read_text(),original)
        self.assertEqual(sum(int(r['records']) for r in read_table(ledger)),3)

    def test_native_thresholds_remain_exactly_configured(self):
        cfg=self.root/'config.xml';cfg.write_text('<config><graphmap delFilter="10000000" delFilterThreshold="0.01" delFilterQuerySizeThreshold="2"/></config>')
        cmd=deletion_command(self.root,'source.paf',cfg,6)
        self.assertEqual(cmd[-10:],['-d','10000000','-v','-p','-t','6','-m','0.01','-s','2'])

    def test_nonreproducing_native_execution_stops_before_assignment(self):
        p=self.root/'prepared';p.mkdir()
        (p/'preparation.json').write_text(json.dumps(dict(expected_native_sha256='mismatch')))
        (p/'pre_deletion.paf').write_text(self.paf('a'))
        (p/'effective_config.xml').write_text('<config><graphmap delFilter="10000000"/></config>')
        modules={'cactus.shared.configWrapper':SimpleNamespace(ConfigWrapper=object),
                 'cactus.refmap.cactus_graphmap':SimpleNamespace(filter_paf=object)}
        def native(cmd,stdout,stderr,check): stdout.write(self.paf('a'))
        with patch.dict(sys.modules,modules),patch('pangenome_deletion_replay.check_tools'),patch('pangenome_deletion_replay.subprocess.run',native):
            with self.assertRaisesRegex(ValueError,'does not exactly reproduce'):
                run_filter(p,'native',self.root/'native',1)

    def test_matrix_report_requires_native_reproduction_and_all_six_cases(self):
        filters=[];assignments=[]
        source=dict(assembly='island',sequence='chr8+chr10',source_bp='100',chromosome_scale='yes',
                    status='represented_in_clip',clip_retained_path_bp='50')
        for policy in POLICIES:
            d=self.root/policy;d.mkdir();filters.append(d)
            (d/'preparation.json').write_text(json.dumps(dict(policy=policy,native_exact_match=policy=='native')))
            write_table(d/'selection.tsv',[source],list(source))
            rows=[dict(policy=policy,assembly='island',sequence='chr8+chr10',stage=s,source_bp=100,covered_bp=50,records=1)
                  for s in ('pre_deletion','post_deletion','post_paf_filter')]
            write_table(d/'filter_coverage.tsv',rows,list(rows[0]))
            for case in ASSIGNMENTS:
                p=d/case;p.mkdir();assignments.append(p)
                (p/'case.json').write_text(json.dumps(dict(policy=policy,case=case)))
                row=dict(case=case,assembly='island',sequence='chr8+chr10',surviving_query_covered_bp=50,assigned_chromosomes='chr10')
                write_table(p/'assignments.tsv',[row],list(row))
        out=self.root/'summary';report(filters,assignments,out)
        self.assertIn('byte for byte',(out/'deletion_filter_report.md').read_text())
        self.assertEqual(len(read_table(out/'assignment_comparison.tsv')),6)
        with self.assertRaisesRegex(ValueError,'Incomplete or duplicated'):
            report(filters,assignments[:-1],self.root/'incomplete')
        (filters[0]/'preparation.json').write_text(json.dumps(dict(policy='native',native_exact_match=False)))
        with self.assertRaisesRegex(ValueError,'non-reproducing'):
            report(filters,assignments,self.root/'mismatch')


if __name__=='__main__':unittest.main()
