"""Launcher bootstrap/packaging checks; controlled tools do not validate YaHS biology."""
import json
import os
from pathlib import Path
import subprocess
import tarfile
import tempfile
import unittest


class ReplayLauncherTests(unittest.TestCase):
    def test_module_exports_conda_and_retained_bin_has_distinct_name(self):
        repo = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root/'gcl_genome_assembly').symlink_to(repo, target_is_directory=True)
            tools = root/'tools';tools.mkdir()
            conda_tools = root/'conda-tools';conda_tools.mkdir()
            conda = conda_tools/'conda'
            conda.write_text('#!/bin/bash\necho /fixture/conda\n');conda.chmod(0o755)
            startup = root/'startup.sh'
            startup.write_text('module() {\n'
                               '  [[ "$1 $2" == "load miniconda3" ]] || return 1\n'
                               '  export PATH="$FIXTURE_CONDA_TOOLS:$PATH"\n}\n')
            nextflow = tools/'nextflow'
            nextflow.write_text('''#!/bin/bash
set -euo pipefail
conda info --base >/dev/null
while [[ $# -gt 0 ]]; do
  case "$1" in
    --replay_links) links="$2"; shift;;
    --outdir) out="$2"; shift;;
  esac
  shift
done
[[ "${links##*/}" == retained_hic_links.bin ]]
[[ -s "$links" ]]
mkdir -p "$out"
printf '{}\n' > replay_check.json
''')
            nextflow.chmod(0o755)
            fasta = root/'input.fa';fasta.write_text('>s\nACGT\n')
            links = root/'fixture_scaffolds_final.bin';links.write_bytes(b'links')
            agp = root/'expected.agp';agp.write_text('s\t1\t4\t1\tW\ts\t1\t4\t+\n')
            out = root/'output'
            env = dict(os.environ, SLURM_SUBMIT_DIR=str(root), BASH_ENV=str(startup),
                       FIXTURE_CONDA_TOOLS=str(conda_tools), PATH=str(tools)+':'+os.environ['PATH'])
            result = subprocess.run(['bash', str(repo/'scripts/comparisons/run_scaffolding_replay.sbatch'),
                'fixture', str(fasta), str(fasta), str(links), str(links), str(agp), str(out)],
                env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            self.assertEqual(result.returncode, 0, result.stdout)
            self.assertEqual((out/'retained_hic_links.bin').resolve(), links)
            self.assertTrue(json.loads((out/'input_identity.json').read_text())['links_match'])
            with tarfile.open(str(out)+'.review.tar.gz') as archive:
                self.assertIn('replay_check.json', archive.getnames())
                self.assertFalse(any(n.endswith('.bin') for n in archive.getnames()))


if __name__ == '__main__':
    unittest.main()
