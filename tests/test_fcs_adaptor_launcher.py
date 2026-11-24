from pathlib import Path
import json
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from run_fcs_adaptor import launch
from setup_fcs_adaptor_host_checks import setup


class LauncherTests(unittest.TestCase):
    def test_shellless_image_is_given_application_directly(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            runtime=setup(root)
            fasta=root/'input with spaces.fasta';fasta.write_text('>c1\nACGT\n')
            with patch('run_fcs_adaptor.shutil.which',return_value=str(runtime)):
                launch(fasta,root/'output','docker://ncbi/fcs-adaptor:0.5.5','euk','singularity',
                       root/'screen.log',root/'launch.json')
            metadata=json.loads((root/'launch.json').read_text())
            self.assertIn('/app/fcs/bin/av_screen_x',metadata['command'])
            self.assertNotIn('/bin/bash',metadata['command'])
            self.assertNotIn('/bin/sh',metadata['command'])
            self.assertEqual(metadata['command'][-1],'/sample-volume/input with spaces.fasta')
            self.assertEqual(len(metadata['source_sha256']),64)
            self.assertTrue((root/'output/fcs_adaptor_report.txt').is_file())

    def test_missing_runtime_fails_before_screening(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);fasta=root/'a.fa';fasta.write_text('>c1\nACGT\n')
            with patch('run_fcs_adaptor.shutil.which',return_value=None):
                with self.assertRaisesRegex(RuntimeError,'required on the host'):
                    launch(fasta,root/'out','docker://ncbi/fcs-adaptor:0.5.5','euk','singularity',
                           root/'log',root/'manifest')

    def test_existing_local_image_is_used_directly(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);runtime=setup(root)
            fasta=root/'input.fasta';fasta.write_text('>c1\nACGT\n')
            image=root/'cached.img';image.write_bytes(b'fixture image')
            with patch('run_fcs_adaptor.shutil.which',return_value=str(runtime)):
                launch(fasta,root/'out',str(image),'euk','singularity',root/'log',root/'manifest')
            metadata=json.loads((root/'manifest').read_text())
            self.assertEqual(metadata['image'],str(image.resolve()))
            self.assertIn(str(image.resolve()),metadata['command'])


if __name__=='__main__':unittest.main()
