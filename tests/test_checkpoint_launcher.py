import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


class CheckpointLauncher(unittest.TestCase):
    def test_no_project_biology_and_explicit_overrides_survive_launch_directory(self):
        repo=Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'gcl_genome_assembly').symlink_to(repo,target_is_directory=True)
            tools=root/'tools';tools.mkdir()
            startup=root/'startup.sh';startup.write_text('module() { return 0; }\n')
            tool=tools/'nextflow'
            tool.write_text('#!/usr/bin/env python3\nimport json,os,sys\nfrom pathlib import Path\nPath(os.environ["CAPTURE"]).write_text(json.dumps(sys.argv[1:]))\n')
            tool.chmod(0o755);capture=root/'args.json'
            env=dict(os.environ,SLURM_SUBMIT_DIR=str(root),BASH_ENV=str(startup),CAPTURE=str(capture),
                PATH=str(tools)+os.pathsep+os.environ['PATH'])
            script=repo/'run_assembly_checkpoint.sbatch'
            result=subprocess.run(['bash',str(script),'data/samples.csv','NONE','validate',
                '--outdir','output folder'],env=env,text=True,capture_output=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            args=json.loads(capture.read_text())
            for forbidden in ['--taxid','--yahs_round1_enzyme','--inspector_contig_min_contig_bp','--hic_readsets']:
                self.assertNotIn(forbidden,args)
            self.assertEqual(args[args.index('--outdir')+1],str(root/'output folder'))
            self.assertEqual(args[args.index('--publish_dir_mode')+1],'symlink')
            result=subprocess.run(['bash',str(script),'samples.csv','hic.csv','validate',
                '--taxid','123','--ploidy','4'],env=env,text=True,capture_output=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            args=json.loads(capture.read_text())
            self.assertEqual(args[args.index('--taxid')+1],'123')
            self.assertEqual(args[args.index('--ploidy')+1],'4')
            self.assertEqual(args[args.index('--hic_readsets')+1],str(root/'hic.csv'))


if __name__=='__main__':unittest.main()
