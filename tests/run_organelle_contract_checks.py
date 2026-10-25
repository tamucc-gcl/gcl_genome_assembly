"""Bounded production-workflow tests with native assemblers/database downloaders doubled."""
from pathlib import Path
import os
import subprocess
repo=Path(__file__).resolve().parents[1]
checks=[(c,True,'PASS: recovery states') for c in
        ['circular','linear','incomplete','no_recovery','retry','fungus_retry','multi_target','multi_species']]+[
    ('execution_failure',False,'17'),
    ('logged_error',False,'GetOrganelle reported an error'),
    ('incomplete_log',False,'completion marker missing'),
    ('database_failure',False,'17'),
    ('mitos_database_failure',False,'17'),
    ('annotation_failure',False,'19'),
    ('missing_bed',False,'annotation is incomplete'),
    ('multiple_records',False,'Current MITOS adapter requires one FASTA record'),
    ('missing_targets',False,'Missing required organelle targets'),
    ('missing_gcode',False,'Missing required genetic code'),
    ('duplicate_targets',False,'Invalid or duplicate organelle targets'),
    ('missing_reference',False,'Missing required mitochondrial reference'),
    ('hifi_missing_gcode',False,'Missing required genetic code'),
]
out=repo/'comparisons/organelle-contract-validation';out.mkdir(parents=True,exist_ok=True)
for case,success,message in checks:
    result=subprocess.run(['bash',str(repo/'tests/run_chimera_channel_reproduction.sh'),
        'organelle_contract_checks.nf',case],cwd=repo,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (out/(case+'.log')).write_text(result.stdout)
    failed=result.returncode!=0 or 'ERROR ~' in result.stdout
    assert failed != success and message in result.stdout,f'{case}: unexpected outcome; see {out/(case+".log")}'
    print('PASS',case,flush=True)
for case,success,message in [('taxonomy_success',True,'PASS: taxonomy resolution'),
                             ('unknown_taxid',False,'did not resolve'),
                             ('empty_lineage',False,'did not resolve'),
                             ('taxonomy_failure',False,'17')]:
    result=subprocess.run(['bash',str(repo/'tests/run_chimera_channel_reproduction.sh'),
        'taxonomy_contract_checks.nf',case],cwd=repo,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (out/(case+'.log')).write_text(result.stdout)
    failed=result.returncode!=0 or 'ERROR ~' in result.stdout
    assert failed != success and message in result.stdout,f'{case}: unexpected outcome; see {out/(case+".log")}'
    print('PASS',case,flush=True)
