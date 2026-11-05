import json
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/comparisons'))
from prepare_integration_params import prepare


class IntegrationSettings(unittest.TestCase):
    def test_baseline_settings_are_explicit_and_retired_or_run_paths_are_not_restored(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'provenance').mkdir()
            data=dict(taxid=123,inspector_min=10,outdir='old',sample_sheet='old.csv',pangenome_private_analysis=True)
            (root/'provenance/run_provenance.json').write_text(json.dumps(dict(context=dict(parameters=data))))
            config=root/'config';config.write_text('manifest {\n    version = "core"\n}\nparams {\n    taxid = null\n    inspector_min = 1\n    outdir = "new"\n    sample_sheet = null\n}\n')
            values=prepare(root,config,root/'params.json')
            self.assertEqual(values,dict(taxid=123,inspector_min=10))


if __name__=='__main__':unittest.main()
