import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

spec=importlib.util.spec_from_file_location('core_audit',Path(__file__).resolve().parents[1]/'scripts/comparisons/audit_core_integration.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)


class IntegrationAuditTests(unittest.TestCase):
    def setup_run(self,root):
        baseline=root/'baseline';results=root/'results';out=root/'audit'
        baseline.mkdir();results.mkdir();out.mkdir()
        rows=[]
        for hap in ['hap1','hap2']:
            fai=root/(hap+'.fai');fai.write_text('chr1\t1000\t0\t1000\t1001\n')
            rows.append(dict(id=hap,sample='individual',taxid='1',completion='finalized',haplotype=hap,
                             ploidy=2,fasta=str(root/(hap+'.fasta')),fai=str(fai)))
        manifest=dict(assemblies=rows,cohorts=[dict(taxid='1',ready=True)])
        for base in [baseline,results]:
            (base/'assembly_eligibility.json').write_text(json.dumps(manifest))
            graph=base/'pangenome/1/representation';graph.mkdir(parents=True)
            for name in ['clipping_audit.json','sequence_graph_audit.json','source_audit.json']:
                (graph/name).write_text(json.dumps(dict(status='REVIEW' if name=='source_audit.json' else 'PASS')))
            ref=base/'references/1';ref.mkdir(parents=True)
            (ref/'reference.json').write_text(json.dumps(dict(id='hap1',sha256='abc',sequences=[dict(name='chr1',length=1000)],
                                                            coordinate_stage='final',sequence_count=1,bases=1000)))
        comparison=results/'pairwise_alignments';comparison.mkdir()
        (comparison/'comparison_plan.tsv').write_text('ref_id\tquery_id\nhap1\thap2\n')
        (comparison/'hap1_vs_hap2.plot_status.tsv').write_text('plot\tstatus\ndotplot\tplotted\nribbon\tplotted\n')
        for suffix in ['_dotplot.png','_riparian.png']:(comparison/('hap1_vs_hap2'+suffix)).write_bytes(b'fixture')
        for name in ['assembly_run_summary.md','assembly_report.md','reports/assembly_qc_report.html',
                     'provenance/run_provenance.json','provenance/software_versions.tsv','provenance/dependency_catalog.tsv',
                     'provenance/software_version_issues.tsv']:
            p=results/name;p.parent.mkdir(exist_ok=True);p.write_text('fixture\n')
        trace=root/'trace.tsv';trace.write_text('name\tstatus\texit\tworkdir\nPANGENOME:CACTUS_PANGENOME (1)\tCACHED\t0\t-\n')
        return baseline,results,trace,out

    def test_accepts_cached_core_and_existing_review_disposition(self):
        with tempfile.TemporaryDirectory() as tmp:
            result=mod.audit(*self.setup_run(Path(tmp)))
            self.assertEqual(result['status'],'PASS');self.assertEqual(result['comparisons'],1)

    def test_rejects_ignored_failure_and_graph_recomputation(self):
        with tempfile.TemporaryDirectory() as tmp:
            args=self.setup_run(Path(tmp));args[2].write_text('name\tstatus\texit\tworkdir\nPANGENOME:CACTUS_PANGENOME (1)\tCOMPLETED\t0\t-\nOPTIONAL\tFAILED\t1\t-\n')
            result=mod.audit(*args)
            self.assertEqual(result['status'],'FAIL');self.assertTrue(any('cached' in e for e in result['errors']))

    def test_rejects_final_input_or_graph_admission_change(self):
        with tempfile.TemporaryDirectory() as tmp:
            args=self.setup_run(Path(tmp));path=args[1]/'assembly_eligibility.json';data=json.loads(path.read_text())
            data['assemblies'][0]['fasta']='changed.fasta';data['cohorts'][0]['ready']=False;path.write_text(json.dumps(data))
            result=mod.audit(*args);self.assertEqual(result['status'],'FAIL')
            self.assertTrue(any('admission' in e for e in result['errors']));self.assertTrue(any('fasta' in e for e in result['errors']))


if __name__=='__main__':unittest.main()
