import csv
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from psmc_plan import plan,manifest_rows
from psmc_results import parse,export,resample,pattern_groups
from psmc_report_index import index
from psmc_input import bed,genotype,genotype_confidence


def fitted(iteration=2):
    return ('MM\tpattern:2+2, n:4, skip:1,\nRD\t'+str(iteration)+'\nRI\t0\nTR\t0.1\t0.02\n'+
            ''.join(f'RS\t{i}\t{i*.5}\t2\t0.2\t0.1\t0.1\n' for i in range(4))+'//\n')


class PsmcContracts(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
    def tearDown(self):self.tmp.cleanup()
    def row(self,sample='a',hap='hap1',ploidy=2,taxid='1'):
        return dict(sample=sample,haplotype=hap,ploidy=ploidy,taxid=taxid,completion='finalized',
                    representation='hic_phased',fasta=hap+'.fa',species='Species')
    def refs(self):return [dict(taxid='1',id='ref',coordinate_stage='final',sha256='a'*64)]
    def test_pair_identity_not_graph_admission(self):
        p=plan(dict(assemblies=[self.row(),self.row(hap='hap2')]),self.refs(),[])
        self.assertEqual(p['jobs'][0]['mode'],'assemblies')
        self.assertIsNone(p['jobs'][0]['mutation_rate'])
    def test_individuals_never_paired(self):
        p=plan(dict(assemblies=[self.row(),self.row(sample='b',hap='hap2')]),self.refs(),[])
        self.assertEqual(len(p['skipped']),2)
    def test_single_diploid_full_reads(self):
        p=plan(dict(assemblies=[self.row(hap='primary')]),self.refs(),[dict(sample='a',fastq='complete.fq')])
        self.assertEqual(p['jobs'][0]['mode'],'hifi')
    def test_true_haploid_rejected(self):
        p=plan(dict(assemblies=[self.row(hap='primary',ploidy=1)]),self.refs(),[dict(sample='a',fastq='reads')])
        self.assertIn('biologically diploid',p['skipped'][0]['reason'])
    def test_no_hifi_reasoned_skip(self):
        self.assertEqual(len(plan(dict(assemblies=[self.row(hap='primary')]),self.refs(),[])['skipped']),1)
    def test_no_long_reference_scaffolds_skips_before_calling(self):
        refs=self.refs();refs[0]['sequences']=[dict(name='short',length=1000)]
        p=plan(dict(assemblies=[self.row(),self.row(hap='hap2')]),refs,[])
        self.assertIn('minimum length',p['skipped'][0]['reason'])
    def test_missing_expected_haplotype(self):
        row=self.row(hap='hap2');row['completion']='missing'
        self.assertIn('incomplete',plan(dict(assemblies=[self.row(),row]),self.refs(),[])['skipped'][0]['reason'])
    def test_scaling_positive_paired(self):
        e=dict(assemblies=[self.row(),self.row(hap='hap2')])
        for mu,g in [(0,1),(-1,1),(float('nan'),1),(1,None)]:
            with self.assertRaises(ValueError):plan(e,self.refs(),[],mutation_rate=mu,generation_time=g)
        self.assertEqual(plan(e,self.refs(),[],mutation_rate=1e-8,generation_time=.5)['jobs'][0]['generation_time'],.5)
    def test_unknown_override_not_silently_ignored(self):
        with self.assertRaises(ValueError):plan(dict(assemblies=[self.row()]),self.refs(),[],{'typo':{}})
    def test_reference_final_hash_required(self):
        r=self.refs();r[0]['coordinate_stage']='scaffold'
        with self.assertRaises(ValueError):plan(dict(assemblies=[self.row(),self.row(hap='hap2')]),r,[])
    def test_cross_species_identity_and_scaling(self):
        with self.assertRaises(ValueError):plan(dict(assemblies=[self.row(),self.row(hap='hap2',taxid='2')]),self.refs(),[])
        with self.assertRaises(ValueError):plan(dict(assemblies=[self.row(),self.row(sample='b',taxid='2')]),self.refs(),[],mutation_rate=1e-8,generation_time=1)
    def test_manifest_relative_intervals(self):
        (self.root/'scope.bed').write_text('c\t0\t3\n')
        p=self.root/'input.tsv';p.write_text('sample\tinclude_bed\na\tscope.bed\n')
        self.assertEqual(manifest_rows(p)['a']['include_bed'],str((self.root/'scope.bed').resolve()))
    def test_masks_boundaries_and_merge(self):
        p=self.root/'scope.bed';p.write_text('c\t0\t3\nc\t2\t5\n')
        self.assertEqual(bed(p,{'c':5}),{'c':[(0,5)]})
        p.write_text('c\t0\t6\n')
        with self.assertRaises(ValueError):bed(p,{'c':5})
    def test_genotypes_missing_and_complex(self):
        self.assertEqual(genotype(('A','G'),(0,1)),'R')
        self.assertEqual(genotype(('A','G'),(1,1)),'G')
        self.assertIsNone(genotype(('A','AG'),(0,1)))
        self.assertIsNone(genotype(('A','G'),(None,1)))
    def test_homozygous_likelihood_confidence_is_not_missing_gq(self):
        self.assertEqual(genotype_confidence(dict(GT=(0,0),PL=(0,54,213))),(54.,'PL_separation'))
        self.assertEqual(genotype_confidence(dict(GT=(0,1),PL=(54,0,213))),(54.,'PL_separation'))
        self.assertEqual(genotype_confidence(dict(GT=(0,0),GQ=40)),(40.,'GQ'))
        self.assertIsNone(genotype_confidence(dict(GT=(0,0)))[0])
    def test_complete_and_truncated_iterations(self):
        p=self.root/'observed.psmc';p.write_text(fitted())
        self.assertEqual(parse(p)['blocks'][0]['iteration'],2)
        p.write_text(fitted().replace('//\n',''))
        with self.assertRaisesRegex(ValueError,'Truncated'):parse(p)
    def test_nonfinite_and_interval_count_rejected(self):
        p=self.root/'observed.psmc'
        for text in [fitted().replace('0.1\t0.02','nan\t0.02'),fitted().replace('RS\t3\t1.5\t2\t0.2\t0.1\t0.1\n','')]:
            p.write_text(text)
            with self.assertRaises(ValueError):parse(p)
    def test_explicit_fit_set_no_stale_bootstraps(self):
        p=self.root/'observed.psmc';p.write_text(fitted())
        with self.assertRaises(ValueError):export([p],self.root/'out',2,100,expected_bootstraps=1)
    def test_scaling_nondefault_bin_fractional_generation(self):
        p=self.root/'observed.psmc';p.write_text(fitted())
        rows=export([p],self.root/'out',2,50,1e-8,.5)
        self.assertEqual(len(rows),2) # terminal free-parameter block excluded, as upstream
        self.assertAlmostEqual(rows[1]['time'],25000)
        self.assertAlmostEqual(rows[1]['population'],100000)
        rows=export([p],self.root/'unscaled',2,50)
        self.assertAlmostEqual(rows[1]['time'],.001)
        self.assertAlmostEqual(rows[1]['population'],.004)
    def test_pattern_validation(self):
        self.assertEqual(pattern_groups('4+2*3+6'),[4,3,3,6])
        for p in ['0+4','4; rm -rf','256','2*0']:
            with self.assertRaises(ValueError):pattern_groups(p)
    def test_seeded_bootstrap_reproducible(self):
        source=self.root/'split.fa';source.write_text('>a\nKKTTNN\n>b\nKTTTTN\n>c\nKKKKTT\n')
        import os
        before=Path.cwd()
        try:
            os.chdir(self.root)
            a=resample(source,'one.fa','sample',1,1729)
            b=resample(source,'two.fa','sample',1,1729)
            c=resample(source,'three.fa','sample',2,1729)
            self.assertEqual(a,b);self.assertNotEqual(a['derived_seed'],c['derived_seed'])
            self.assertEqual(Path('one.fa').read_bytes(),Path('two.fa').read_bytes())
        finally:os.chdir(before)
    def test_missing_individual_report_rejected(self):
        p=self.root/'plan.json';p.write_text(json.dumps(dict(jobs=[dict(sample='a')],skipped=[])))
        with self.assertRaises(ValueError):index(p,[],self.root/'out')


if __name__=='__main__':unittest.main()
