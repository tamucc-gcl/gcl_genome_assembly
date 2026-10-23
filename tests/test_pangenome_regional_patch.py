import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
import pangenome_regional_patch as mod
from pangenome_source_audit import write_table


class RegionalPatchTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.sources=[dict(assembly='sample',sequence='piece',source_bp=100)]

    def region(self,lo,hi,chrom='A'):
        return dict(assembly='sample',sequence='piece',chromosome=chrom,start_0based=lo,end_exclusive=hi)

    def test_source_is_pinned_and_patch_is_one_line(self):
        manifest,text=mod.verified_source(Path(__file__).resolve().parents[1]/'vendor/cactus-gfa-tools-regional')
        self.assertEqual(manifest['commit'],mod.COMMIT)
        self.assertEqual(text.count(mod.OLD),1)
        before=text.splitlines();after=text.replace(mod.OLD,mod.NEW).splitlines()
        self.assertEqual(sum(a!=b for a,b in zip(before,after)),1)

    def test_adjacent_chromosomes_are_not_duplicates(self):
        row=mod.coordinate_audit([self.region(0,40),self.region(40,100,'B')],self.sources)[0]
        self.assertEqual(row['assigned_union_bp'],100)
        self.assertEqual(row['assigned_sum_bp'],100)

    def test_overlaps_fail_even_between_different_chromosomes(self):
        for chrom in ('A','B','_AMBIGUOUS_'):
            with self.subTest(chrom=chrom),self.assertRaisesRegex(ValueError,'Overlapping'):
                mod.coordinate_audit([self.region(0,60),self.region(50,100,chrom)],self.sources)

    def test_duplicate_regions_fail(self):
        with self.assertRaisesRegex(ValueError,'duplicate'):
            mod.coordinate_audit([self.region(0,50)]*2,self.sources)

    def test_out_of_bounds_or_empty_regions_fail(self):
        for lo,hi in ((-1,50),(0,101),(40,40),(60,40)):
            with self.subTest(lo=lo,hi=hi),self.assertRaisesRegex(ValueError,'outside'):
                mod.coordinate_audit([self.region(lo,hi)],self.sources)

    def test_duplicate_sources_fail(self):
        with self.assertRaisesRegex(ValueError,'Duplicate source'):
            mod.coordinate_audit([],self.sources*2)

    def test_ambiguous_and_unreported_are_explicit(self):
        row=mod.coordinate_audit([self.region(0,30),self.region(40,60,'_AMBIGUOUS_')],self.sources)[0]
        self.assertEqual((row['assigned_union_bp'],row['ambiguous_bp'],row['unreported_bp']),(30,20,50))

    def test_no_assignment_does_not_claim_source_coverage(self):
        self.assertEqual(mod.coordinate_audit([],self.sources)[0]['unreported_bp'],100)

    def paf(self,name='id=IND_sample.1|piece:41-100',length=60,lo=0,hi=60):
        return f'{name}\t{length}\t{lo}\t{hi}\t+\ts1\t100\t0\t60\t60\t60\t60\n'

    def membership(self,text,chrom='B'):
        identities=self.root/'identities.tsv'
        identities.write_text('id\tgraph_name\nsample\tIND_sample.1\n')
        (self.root/('split_'+chrom+'.paf')).write_text(text)
        mod.audit_paf_membership(self.root,[self.region(40,100,'B')],self.sources,identities)

    def test_mapping_coordinates_are_relative_to_the_assigned_fragment(self):
        self.membership(self.paf())

    def test_mapping_in_another_chromosome_bin_fails(self):
        with self.assertRaisesRegex(ValueError,'does not belong'):self.membership(self.paf(),chrom='A')

    def test_mapping_outside_fragment_fails(self):
        for text in (self.paf(hi=61),self.paf(length=100),self.paf(lo=-1)):
            with self.subTest(text=text),self.assertRaisesRegex(ValueError,'outside'):self.membership(text)

    def test_unknown_biological_query_fails(self):
        with self.assertRaisesRegex(ValueError,'Unknown'):self.membership(self.paf(name='unknown'))

    def test_binary_provenance_mismatch_fails_before_assignment(self):
        binary=self.root/'binary';binary.write_bytes(b'wrong')
        info=self.root/'build.json';info.write_text(json.dumps(dict(patch=mod.PATCH,patched_binary_sha256='wrong')))
        with self.assertRaisesRegex(ValueError,'provenance'):
            mod.assign(self.root,binary,info,[],self.root/'output')

    def matrix(self):
        original=[];fixed=[];filters=[]
        for policy in ('native','chromosome_isolated','disabled'):
            for case,dirs in (('regional_1Mb',original),(mod.CASE,fixed)):
                p=self.root/(policy+'_'+case);p.mkdir();dirs.append(p)
                (p/'case.json').write_text(json.dumps(dict(policy=policy,case=case,coordinate_integrity='PASS')))
                row=dict(assembly='sample',sequence='piece',filtered_query_covered_bp=100,
                    surviving_query_covered_bp=100,assigned_chromosomes='A,B',ambiguous_span_bp=0)
                write_table(p/'assignments.tsv',[row],list(row))
                audit=mod.coordinate_audit([self.region(0,100)],self.sources)
                write_table(p/'coordinate_audit.tsv',audit,list(audit[0]))
            p=self.root/('filter_'+policy);p.mkdir();filters.append(p)
            row=dict(**self.sources[0],chromosome_scale='yes',status='represented_in_clip',clip_retained_path_bp=100)
            write_table(p/'selection.tsv',[row],list(row))
        return original,fixed,filters

    def test_complete_report_and_coordinate_accounting(self):
        mod.report(*self.matrix(),self.root/'report')
        self.assertIn('3 source/policy audits passed',(self.root/'report/regional_assignment_report.md').read_text())

    def test_missing_policy_fails(self):
        old,new,filters=self.matrix()
        with self.assertRaisesRegex(ValueError,'policy matrix'):mod.report(old,new[:-1],filters,self.root/'report')

    def test_missing_source_audit_fails(self):
        old,new,filters=self.matrix()
        (new[0]/'coordinate_audit.tsv').write_text('assembly\tsequence\n')
        with self.assertRaisesRegex(ValueError,'audit matrix'):mod.report(old,new,filters,self.root/'report')

    def test_invalid_source_accounting_fails(self):
        old,new,filters=self.matrix()
        p=new[0]/'coordinate_audit.tsv';p.write_text(p.read_text().replace('\tPASS','\tFAIL'))
        with self.assertRaisesRegex(ValueError,'accounting'):mod.report(old,new,filters,self.root/'report')


if __name__=='__main__':unittest.main()
