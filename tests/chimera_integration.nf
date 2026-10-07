nextflow.enable.dsl=2

params.fixture = null
params.scenario = 'unresolved'
params.outdir = 'chimera-integration-output'
params.publish_dir_mode = 'copy'
params.chimera_break = false
params.chimera_min_piece_bp = 4

// Set defaults before including modules, which capture parameter values.
include { CHIMERA_REVIEW } from '../modules/chimera_review.nf'
include { BREAK_CHIMERAS } from '../modules/break_chimeras.nf'

process VERIFY_RESULT {
    publishDir "${params.outdir}/verification", mode: 'copy'
    input:
    tuple val(meta), path(fasta), path(names)
    val(scenario)
    output:
    path('PASS.json')
    script:
    """
    python3 - '${fasta}' '${names}' '${scenario}' <<'PY'
import csv, json, sys
from pathlib import Path
records = {}
for line in Path(sys.argv[1]).read_text().splitlines():
    if line.startswith('>'):
        name = line[1:]; records[name] = ''
    else:
        records[name] += line
with open(sys.argv[2]) as handle:
    rows = list(csv.DictReader(handle, delimiter='\t'))
assert set(records) == {r['old_name'] for r in rows}
assert sum(map(len, records.values())) == 14
assert records['other'] == 'GGGG'
if sys.argv[3] == 'unresolved':
    assert records['composite'] == 'AAAANNCCCC' and len(records) == 2
else:
    assert records['composite_sub_0_6'] == 'AAAANN'
    assert records['composite_sub_6_10'] == 'CCCC'
    assert len(records) == 3
    for row in rows:
        if row['old_name'].startswith('composite_sub_'):
            assert row['new_name'] == row['old_name']
            assert row['class'] == 'unplaced' and row['orient'] == '+'
Path('PASS.json').write_text(json.dumps({'scenario': sys.argv[3], 'status': 'PASS'}))
PY
    """
}

workflow {
    if (!params.fixture) error 'Supply --fixture directory'
    def root = file(params.fixture)
    def scripts = file("${projectDir}/../py_scripts")
    def meta = [id:'synthetic', sample:'synthetic', taxid:'synthetic']
    CHIMERA_REVIEW(Channel.of(tuple(meta, file("${root}/calls.tsv"), file("${root}/context"))),
                      file("${scripts}/chimera_review.py"))
    selected = Channel.of(tuple(meta, file("${root}/manual-actions.tsv")))
    BREAK_CHIMERAS(selected.map { m, actions -> tuple(m, file("${root}/assessment.fa"), file("${root}/names.tsv"), actions) },
                  Channel.value([file("${scripts}/break_chimeras.py"), file("${scripts}/chimera_actions.py")]))
    VERIFY_RESULT(BREAK_CHIMERAS.out.assemblies, params.scenario)
}
