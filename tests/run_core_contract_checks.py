#!/usr/bin/env python3
"""Local Nextflow 23.10 checks. External tool doubles do not validate graph biology."""
import subprocess
import os
from pathlib import Path

repo = Path(__file__).resolve().parents[1]
checks = [
    ('input_reference_checks.nf', 'success', True, 'PASS: production shared reference process'),
    ('input_reference_checks.nf', 'duplicate_sample', False, 'Duplicate sample_id'),
    ('input_reference_checks.nf', 'duplicate_file', False, 'Read file assigned twice'),
    ('input_reference_checks.nf', 'unknown_hic_sample', False, 'unknown sample'),
    ('input_reference_checks.nf', 'no_valid_samples', False, 'No assembly-capable valid samples'),
    ('core_contract_checks.nf', 'success', True, 'PASS: validated graph roles'),
    ('core_contract_checks.nf', 'assignment_inputs', True, 'PASS: cached assignment inputs retain species-specific'),
    ('core_contract_checks.nf', 'graph_failure', False, 'exit status (17)'),
    ('core_contract_checks.nf', 'fragment_header_failure', False, 'Installed Cactus fragment FASTA convention differs'),
    ('core_contract_checks.nf', 'cross_species', False, 'Cross-species graph input'),
    ('core_contract_checks.nf', 'missing_graph_role', False, 'Graph output audit failed'),
    ('core_contract_checks.nf', 'organelle_failure', False, 'exit status (17)'),
]
for fixture, case, success, message in checks:
    result = subprocess.run(['bash', str(repo / 'tests/run_chimera_channel_reproduction.sh'), fixture, case],
                            cwd=repo, env=dict(os.environ, CORE_INPUT_CASE=case),
                            text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    log = repo / 'comparisons' / ('core-check-' + fixture.removesuffix('.nf') + '-' + case + '.log')
    log.write_text(result.stdout)
    assert (result.returncode == 0) == success, f'{case}: unexpected exit {result.returncode}; see {log}'
    assert message in result.stdout, f'{case}: expected diagnostic missing; see {log}'
    if success:
        assert 'ERROR ~' not in result.stdout, f'{case}: Nextflow error; see {log}'
    print(f'PASS {fixture}: {case}', flush=True)
