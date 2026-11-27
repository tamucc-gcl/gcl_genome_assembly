"""Combine current-run stage tables for the cross-stage FCS figure."""
import csv
import argparse
import json
from pathlib import Path
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('directories',nargs='+',type=Path,help='Exact staged stage-report directories')
args=parser.parse_args()
root=Path('overview');root.mkdir(exist_ok=True)
rows=[];seen=set()
# Path.rglob does not descend into directory symlinks on supported task Python
# versions. Nextflow stages report directories that way: open each exact input.
paths=[directory/'cohort_fcs_sequence_evidence.tsv' for directory in args.directories]
for path in paths:
    with path.open() as f:current=list(csv.DictReader(f,delimiter='\t'))
    for row in current:
        if 'cat_primary_bp' not in row:
            raise ValueError('Stage table lacks FCS category evidence: '+str(path))
        key=(row['cohort_assembly'],row['coordinate_stage'],row['sequence'])
        if key in seen:raise ValueError('Duplicate FCS stage record: '+str(key))
        seen.add(key);rows.append(row)
if rows:
    with (root/'cohort_fcs_sequence_evidence.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
stages=sorted({r['coordinate_stage'] for r in rows})
(root/'overview_scope.json').write_text(json.dumps(dict(input_tables=[str(p) for p in paths],
    records=len(rows),assemblies=sorted({r['cohort_assembly'] for r in rows}),stages=stages,
    records_by_stage={stage:sum(r['coordinate_stage']==stage for r in rows) for stage in stages}),indent=2)+'\n')
(root/'contamination_report.md').write_text('# FCS evidence across assembly stages\n\n'+(
    '![Mixed-record category composition](cohort_fcs_composition.png)\n\n'
    'One row per assembly; available stage panels share row order and category colours. Bars partition '
    'mixed-record base pairs, not the complete assembly. Repeat / low GX coverage remain inconclusive. '
    'Taxonomy flags are not cleaning recommendations.\n\n'
    '[PDF](cohort_fcs_composition.pdf) · [All record evidence](cohort_fcs_sequence_evidence.tsv)\n\n'
    if rows else 'FCS category evidence was not supplied.\n\n')+
    ' · '.join(f'[{stage.title()} reports and GC/coverage](../../assembly/{stage}/decontam/contamination_report.md)'
               for stage in ('contig','scaffold','final') if stage in stages)+'\n',encoding='utf-8')
