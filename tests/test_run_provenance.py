import csv
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from run_provenance import collect,task_context,environment_inventory


class ProvenanceTests(unittest.TestCase):
    def fixture(self,root,values):
        staged=root/'staged';staged.mkdir();specs=[]
        for i,text in enumerate(values,1):
            task=root/('task'+str(i));task.mkdir()
            (task/'.command.run').write_text('# NEXTFLOW TASK: WORKFLOW:TOOL ('+str(i)+')\n')
            (task/'.command.sh').write_text('tool --version\n')
            target=task/'versions.tsv';target.write_text(text);specs.append(dict(path=str(target)))
            dest=staged/('version'+str(i));dest.mkdir();(dest/'versions.tsv').write_text(text)
        return specs,staged

    def test_current_schemas_headers_and_multiple_versions_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);specs,staged=self.fixture(root,['tool\tversion\nminimap2\t2.29\n','process\ttool\tversion\nMAP\tminimap2\t2.30\n'])
            collect(specs,staged,dict(catalog_sources=[]),root/'out')
            with (root/'out/software_versions.tsv').open() as h:rows=list(csv.DictReader(h,delimiter='\t'))
            self.assertEqual([r['version'] for r in rows],['2.29','2.30'])
            self.assertTrue(all(r['tool']=='minimap2' for r in rows));self.assertEqual(rows[1]['process'],'MAP')

    def test_empty_and_unlabeled_sources_are_explicit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);specs,staged=self.fixture(root,['','0.10.3\n','tool\t\n'])
            collect(specs,staged,dict(catalog_sources=[]),root/'out')
            data=json.loads((root/'out/run_provenance.json').read_text())
            self.assertEqual(data['version_records'],0);self.assertEqual(data['issues'],5)

    def test_empty_channel_creates_complete_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'staged').mkdir()
            collect([],root/'staged',dict(catalog_sources=[]),root/'out')
            self.assertEqual(len((root/'out/software_versions.tsv').read_text().splitlines()),1)
            self.assertEqual(json.loads((root/'out/run_provenance.json').read_text())['sources'],[])

    def test_resolved_packages_and_container_fingerprint(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);env=root/'conda env';meta=env/'conda-meta';meta.mkdir(parents=True)
            (meta/'tool.json').write_text(json.dumps(dict(name='samtools',version='1.22',build='abc',url='https://example.org/pkg',sha256='123')))
            image=root/'tool image.sif';image.write_bytes(b'fixture container')
            task=root/'task';task.mkdir()
            (task/'.command.run').write_text(f"# NEXTFLOW TASK: NESTED\nsource /opt/conda/bin/activate '{env}'\nsingularity exec '{image}' /bin/bash .command.sh\n")
            (task/'.command.sh').write_text('samtools --version')
            nested=task/'representation';nested.mkdir()
            context=task_context(nested);data=environment_inventory([context,context])
            self.assertEqual(context['process'],'NESTED');self.assertEqual(len(data['conda']),1)
            self.assertEqual(data['conda'][str(env)]['packages'][0]['version'],'1.22')
            self.assertEqual(data['containers'][str(image)]['sha256'],hashlib.sha256(image.read_bytes()).hexdigest())

    def test_missing_context_not_fabricated(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(task_context(tmp)['context_status'],'unavailable')
            data=environment_inventory([dict(conda_prefixes=[str(Path(tmp)/'absent')],container_files=[])])
            self.assertEqual(next(iter(data['conda'].values()))['status'],'unavailable')


if __name__=='__main__':unittest.main()
