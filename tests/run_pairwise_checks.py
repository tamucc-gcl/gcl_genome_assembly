"""Production visualization workflows; external alignment/QUAST/tidk are controlled doubles."""
from pathlib import Path
import subprocess

repo=Path(__file__).resolve().parents[1]
out=repo/'comparisons/final-viz-validation';out.mkdir(exist_ok=True)
checks=[(case,True,'PASS: species-aware pairing') for case in
        ['within_sample','all_species','no_pairs','disabled','no_alignments']]+[
    ('empty',True,'PASS: empty cohort'),
    ('invalid_mode',False,'Invalid pairwise_alignment_mode'),
    ('duplicate_id',False,'Duplicate final comparison identity'),
    ('conflicting_species',False,'conflicting species'),
    ('alignment_failure',False,'exit status (17)'),
    ('malformed_paf',False,'PAF record has fewer'),
    ('wrong_dictionary',False,'final assembly indexes'),
]
for case,success,message in checks:
    r=subprocess.run(['bash',str(repo/'tests/run_chimera_channel_reproduction.sh'),'pairwise_checks.nf',case],
                     cwd=repo,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (out/(case+'.log')).write_text(r.stdout)
    assert (r.returncode==0)==success and message in r.stdout,f'{case}: see {out/(case+".log")}'
    if case in ['malformed_paf','wrong_dictionary']:assert 'Execution is retried' not in r.stdout
    if case=='alignment_failure':assert 'Execution is retried (1)' in r.stdout
    print('PASS',case,flush=True)
