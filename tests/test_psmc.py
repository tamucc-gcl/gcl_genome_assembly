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
from psmc_reference import verify as verify_reference


def fitted(iteration=2):
    return ('MM\tpattern:2+2, n:4, skip:1,\nRD\t'+str(iteration)+'\nRI\t0\nTR\t0.1\t0.02\n'+
            ''.join(f'RS\t{i}\t{i*.5}\t2\t0.2\t0.1\t0.1\n' for i in range(4))+'//\n')


class PsmcContracts(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
    def tearDown(self):self.tmp.cleanup()
    def row(self,sample='a',hap='hap1',ploidy=2,taxid='1'):
        return dict(id=sample+'_'+hap,sample=sample,haplotype=hap,ploidy=ploidy,taxid=taxid,completion='finalized',
                    representation='hic_phased',fasta=hap+'.fa',species='Species')
    def refs(self):return [dict(sample='a',taxid='1',id='a_hap1',fasta='hap1.fa',
        coordinate_stage='final',sha256='a'*64,sequences=[dict(name='long',length=2000000)])]
    def test_pair_identity_not_graph_admission(self):
        e=dict(assemblies=[self.row(),self.row(hap='hap2')],
            linear_references=[dict(sample='other',id='other_hap1')],cohorts=[dict(ready=False)])
        p=plan(e,self.refs(),[])
        self.assertEqual(p['jobs'][0]['mode'],'assemblies')
        self.assertEqual(p['jobs'][0]['reference']['id'],'a_hap1')
        self.assertIsNone(p['jobs'][0]['mutation_rate'])
    def test_individuals_never_paired(self):
        p=plan(dict(assemblies=[self.row(),self.row(sample='b',hap='hap2')]),self.refs(),[])
        self.assertEqual(len(p['skipped']),2)
    def test_single_diploid_full_reads(self):
        refs=self.refs();refs[0].update(id='a_primary',fasta='primary.fa')
        p=plan(dict(assemblies=[self.row(hap='primary')]),refs,[dict(sample='a',fastq='complete.fq')])
        self.assertEqual(p['jobs'][0]['mode'],'hifi')
    def test_hifi_override_on_paired_individual_keeps_own_hap1_reference(self):
        e=dict(assemblies=[self.row(),self.row(hap='hap2')])
        p=plan(e,self.refs(),[dict(sample='a',fastq='full.fq')],{'a':dict(mode='hifi')})
        self.assertEqual(p['jobs'][0]['mode'],'hifi')
        self.assertEqual(p['jobs'][0]['reference']['id'],'a_hap1')
        self.assertEqual(p['jobs'][0]['reads'],'full.fq')
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
    def test_each_individual_uses_own_reference_despite_shared_species(self):
        rows=[self.row(sample=s,hap=h) for s in ['a','b'] for h in ['hap1','hap2']]
        for r in rows:r['fasta']=r['id']+'.fa'
        refs=[dict(sample=s,taxid='1',id=s+'_hap1',fasta=s+'_hap1.fa',
            coordinate_stage='final',sha256=c*64,sequences=[dict(name='long',length=2000000)]) for s,c in [('a','a'),('b','b')]]
        p=plan(dict(assemblies=rows),refs,[])
        self.assertEqual([j['reference']['id'] for j in p['jobs']],['a_hap1','b_hap1'])
        self.assertNotEqual(p['jobs'][0]['reference']['sha256'],p['jobs'][1]['reference']['sha256'])
        p=plan(dict(assemblies=rows),refs[:1],[])
        self.assertEqual(p['skipped'][0]['sample'],'b') # no species-reference fallback
    def test_foreign_reference_rejected_even_if_mislabeled_as_own(self):
        refs=self.refs();refs[0].update(id='other_hap1',fasta='other.fa')
        with self.assertRaisesRegex(ValueError,'same individual'):
            plan(dict(assemblies=[self.row(),self.row(hap='hap2')]),refs,[])
    def test_individual_reference_checksum_and_fai_validation(self):
        fa=self.root/'own.fa';fa.write_text('>chr1\nACGTNNACGT\n')
        fai=self.root/'own.fa.fai';fai.write_text('chr1\t10\t6\t10\t11\n')
        selection=dict(sample='a',id='a_hap1',taxid='1',selection='own_hap1')
        r=verify_reference(selection,fa,fai)
        self.assertEqual(r['sha256'],hashlib.sha256(fa.read_bytes()).hexdigest())
        self.assertEqual(r['sequences'],[dict(name='chr1',length=10)])
        fai.write_text('chr1\t9\t6\t10\t11\n')
        with self.assertRaisesRegex(ValueError,'disagree'):verify_reference(selection,fa,fai)
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
    def test_cohort_keeps_exported_curve_values_and_reports_unfitted_individuals(self):
        p=self.root/'plan.json';p.write_text(json.dumps(dict(jobs=[dict(sample='a'),dict(sample='b')],skipped=[])))
        dirs=[]
        for sample,status in [('a','completed'),('b','not_fitted')]:
            d=self.root/sample;d.mkdir();dirs.append(d)
            (d/'summary.json').write_text(json.dumps(dict(sample=sample,taxid='1',mode='assemblies',status=status,
                reference='ref',callable_bp=100,heterozygous_snps=2,retained_callable_bp=80,retained_scaffolds=1,
                scaling='unscaled',report='report.md')))
        source='fit\tkind\ttime\tpopulation\tscaling\nobserved\tobserved\t0.001234\t0.005678\tunscaled\n'
        (dirs[0]/'curves.tsv').write_text(source)
        index(p,dirs,self.root/'out')
        with (self.root/'out/combined_curves.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
        self.assertEqual(len(rows),1)
        self.assertEqual((rows[0]['time'],rows[0]['population']),('0.001234','0.005678'))
        self.assertEqual(rows[0]['sample'],'a')
        self.assertEqual((dirs[0]/'curves.tsv').read_text(),source)
        with (self.root/'out/summary.tsv').open() as f:summary=list(csv.DictReader(f,delimiter='\t'))
        self.assertEqual(summary[1]['status'],'not_fitted')
    def test_all_skipped_has_no_phantom_cohort_figure(self):
        p=self.root/'plan.json';p.write_text(json.dumps(dict(jobs=[],skipped=[dict(sample='a',taxid='1',reason='haploid')])))
        index(p,[],self.root/'out')
        self.assertFalse((self.root/'out/combined_curves.tsv').exists())
        self.assertNotIn('combined_history.png',(self.root/'out/README.md').read_text())


if __name__=='__main__':unittest.main()
