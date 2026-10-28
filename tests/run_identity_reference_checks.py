"""Actual name resolution/reference modules; remote reference CLI is a boundary double."""
import subprocess
from pathlib import Path

repo=Path(__file__).resolve().parents[1]
checks=[
    ('identity_mixed',True,'PASS: identity resolution'),
    ('identity_known_only',True,'PASS: identity resolution'),
    ('identity_unknown',False,'Unknown scientific name'),
    ('identity_ambiguous',False,'Ambiguous scientific name'),
    ('identity_missing_node',False,'no node for matched taxid'),
    ('identity_missing',False,'Missing organism identity'),
    ('identity_invalid_taxid',False,'Invalid taxid'),
    ('reference_success',True,'PASS: validated mitochondrial reference pair'),
    ('reference_missing',False,'no unambiguous reference available'),
    ('reference_empty',False,'Empty or mismatched reference files'),
    ('reference_mismatch',False,'Empty or mismatched reference files'),
    ('reference_multiple',False,'no unambiguous reference available'),
    ('reference_header_only',False,'Invalid reference record'),
    ('reference_short',False,'Invalid reference record'),
    ('reference_no_organism',False,'Invalid reference record'),
    ('reference_unavailable_executable',False,'exit status (127)'),
    ('reference_execution_failure',False,'exit status (17)'),
]
out=repo/'comparisons/identity-reference-validation';out.mkdir(exist_ok=True)
for case,success,message in checks:
    r=subprocess.run(['bash',str(repo/'tests/run_chimera_channel_reproduction.sh'),
                      'identity_reference_checks.nf',case],cwd=repo,text=True,
                     stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (out/(case+'.log')).write_text(r.stdout)
    assert (r.returncode==0)==success and message in r.stdout,f'{case}: see {out/(case+".log")}'
    if case=='reference_execution_failure':assert 'Execution is retried (1)' in r.stdout
    if case.startswith('reference_') and case not in ['reference_execution_failure','reference_success']:
        assert 'Execution is retried' not in r.stdout
    print('PASS',case,flush=True)
