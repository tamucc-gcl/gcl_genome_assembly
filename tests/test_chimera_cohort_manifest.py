"""Full-cohort target selection and explicit unavailable-read evidence."""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class CohortManifest(unittest.TestCase):
    def test_all_samples_are_targets_and_missing_bams_are_explicit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);assessment=root/'assessment';harm=assessment/'assembly/harmonization'
            harm.mkdir(parents=True)
            ids=['one_hap1','two_hap1']
            (harm/'1.chromosome_sets.tsv').write_text('id\n'+'\n'.join(ids)+'\n')
            (harm/'1.reference_id.txt').write_text(ids[0])
            (harm/'1.chimera_candidates.tsv').write_text('assembly\n')
            for assembly in ids:
                sample=assembly.rsplit('_hap',1)[0]
                files=[harm/(assembly+'.harmonized_name_map.tsv'),harm/(assembly+'.ref.paf.gz'),
                    assessment/'assembly/scaffold/yahs_round2'/(assembly+'_round2_scaffolds.fa'),
                    assessment/'assembly/scaffold/yahs_round2'/(assembly+'_round2_scaffolds_final.agp'),
                    assessment/'assembly/contig/hifiasm'/(sample+'.hap1.p_ctg.gfa'),
                    assessment/'bam/hic/scaffold/filtered'/(assembly+'.pairs.gz'),
                    assessment/'assembly/scaffold/misassembly_correction'/(assembly+'_corrected.fasta'),
                    assessment/'bam/hic/scaffold/raw'/(assembly+'.readsets.tsv')]
                for f in files:f.parent.mkdir(parents=True,exist_ok=True);f.write_text('fixture\n')
            script=Path(__file__).resolve().parents[1]/'scripts/comparisons/prepare_chimera_real_test.py'
            args=[sys.executable,str(script),'--assessment',str(assessment),'--bam-root',str(root/'bams'),
                  '--species','1','--sample','all','--out',str(root/'inputs')]
            result=subprocess.run(args,capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('Missing retained HiFi BAMs',result.stderr)
            subprocess.run(args+['--allow-missing-bams'],check=True,capture_output=True)
            manifest=json.loads((root/'inputs/manifest.json').read_text())
            self.assertEqual({r['meta']['id'] for r in manifest['targets']},set(ids))
            self.assertEqual(json.loads((root/'inputs/reviewed-actions.json').read_text()),[])
            with (root/'inputs/evidence-availability.tsv').open() as handle:
                availability=list(csv.DictReader(handle,delimiter='\t'))
            self.assertTrue(all(r['hifi'].startswith('Unavailable') for r in availability))
