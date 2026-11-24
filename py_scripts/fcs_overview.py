"""Combine current-run stage tables for the cross-stage FCS figure."""
import csv
from pathlib import Path
root=Path('overview');root.mkdir(exist_ok=True)
rows=[];seen=set()
for path in sorted(Path('inputs').rglob('cohort_fcs_sequence_evidence.tsv')):
    with path.open() as f:current=list(csv.DictReader(f,delimiter='\t'))
    for row in current:
        if 'cat_primary_bp' not in row:continue
        key=(row['cohort_assembly'],row['coordinate_stage'],row['sequence'])
        if key in seen:raise ValueError('Duplicate FCS stage record: '+str(key))
        seen.add(key);rows.append(row)
if rows:
    with (root/'cohort_fcs_sequence_evidence.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
(root/'contamination_report.md').write_text('# FCS evidence across assembly stages\n\n'+(
    '![Mixed-record category composition](cohort_fcs_composition.png)\n\n'
    'One row per assembly; available stage panels share row order and category colours. Bars partition '
    'mixed-record base pairs, not the complete assembly. Repeat / low GX coverage remain inconclusive. '
    'Taxonomy flags are not cleaning recommendations.\n\n'
    '[PDF](cohort_fcs_composition.pdf) · [All record evidence](cohort_fcs_sequence_evidence.tsv)\n\n'
    if rows else 'FCS category evidence was not supplied.\n\n')+
    '[Contig reports and GC/coverage](../../contig/decontam/contamination_report.md) · '
    '[Final reports and GC/coverage](../../final/decontam/contamination_report.md)\n',encoding='utf-8')
