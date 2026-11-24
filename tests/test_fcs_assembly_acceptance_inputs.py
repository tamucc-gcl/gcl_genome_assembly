"""Resolve exact retained task inputs; never substitute post-cleanup/final frames."""
import argparse
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('prepare_fcs',ROOT/'scripts/comparisons/prepare_fcs_assembly_acceptance.py')
prepare=importlib.util.module_from_spec(spec)
spec.loader.exec_module(prepare)


class RetainedInputs(unittest.TestCase):
    def test_configure_binds_original_screen_input_and_actual_final_output(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            trace=root/'trace.tsv'
            assembly='Sde-CTlk_104_hap1'
            sample='Sde-CTlk_104'
            rows=['hash\tname\tstatus\texit']
            paths=[]
            for key,process,identifier in [('aa','FCS_GX_SCREEN',assembly),
                                            ('bb','FINALIZE_ASSEMBLY',assembly),
                                            ('cc','BAM_TO_FASTQ',sample)]:
                path=root/'work'/key/'123456full';path.mkdir(parents=True)
                paths.append(path)
                rows.append(f'{key}/123456\tPREFIX:{process} ({identifier})\tCACHED\t0')
            trace.write_text('\n'.join(rows)+'\n')
            (paths[0]/'original.fa').write_text('>c\naCgT\n')
            (paths[0]/'.command.sh').write_text('/app/bin/run_gx \\\n --fasta original.fa \\\n --tax-id 373251\n')
            (paths[1]/f'{assembly}.fasta').write_text('>chr1\nACGT\n')
            (paths[2]/f'{sample}.fastq.gz').write_bytes(b'retained read placeholder')
            db=root/'db';db.mkdir();(db/'all.gxi').write_bytes(b'db')
            args=argparse.Namespace(project=root,trace=trace,assembly=assembly,taxid=373251,
                source_fasta=None,final_fasta=None,hifi_fastq=None,gxdb_dir=db,out=root/'run')
            prepare.configure(args)
            params=json.loads((args.out/'params.json').read_text())
            self.assertEqual(Path(params['source_fasta']),paths[0]/'original.fa')
            self.assertEqual(Path(params['final_fasta']),paths[1]/f'{assembly}.fasta')
            self.assertFalse(params['run_contamination_evidence'])
            self.assertEqual(params['sample'],sample)
            # Ambiguous successful histories must fail rather than selecting a task.
            with trace.open('a') as handle:
                handle.write(f'dd/654321\tPREFIX:FCS_GX_SCREEN ({assembly})\tCOMPLETED\t0\n')
            with self.assertRaisesRegex(ValueError,'Expected one successful'):
                prepare.configure(args)


if __name__=='__main__':unittest.main()
