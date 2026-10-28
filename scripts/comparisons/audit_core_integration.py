"""Removable acceptance audit of one normal cached main.nf run."""
import argparse
import csv
import json
from itertools import combinations
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'py_scripts'))
from run_provenance import sha256,task_context,environment_inventory


def table(path):
    with Path(path).open() as handle:return list(csv.DictReader(handle,delimiter='\t'))


def audit(baseline,results,trace,out):
    baseline=Path(baseline);results=Path(results);out=Path(out)
    old=json.loads((baseline/'assembly_eligibility.json').read_text())
    new=json.loads((results/'assembly_eligibility.json').read_text())
    before={r['id']:r for r in old['assemblies'] if r['completion']=='finalized'}
    after={r['id']:r for r in new['assemblies'] if r['completion']=='finalized'}
    errors=[];notes=[]
    if {str(c['taxid']) for c in old['cohorts'] if c['ready']}!={str(c['taxid']) for c in new['cohorts'] if c['ready']}:
        errors.append('Graph cohort admission changed')
    if set(before)!=set(after):errors.append('Finalized assembly membership changed')
    for aid in set(before)&set(after):
        for key in ['sample','taxid','haplotype','ploidy','n_hap','fasta','fai','graph_member','graph_name']:
            if before[aid].get(key)!=after[aid].get(key):errors.append(f'{aid}: final {key} changed')
        if sha256(before[aid]['fai'])!=sha256(after[aid]['fai']):errors.append(aid+': final index changed')
        relative=Path('assembly/final')/(aid+'.name_map.tsv')
        if (baseline/relative).is_file():
            if not (results/relative).is_file() or sha256(baseline/relative)!=sha256(results/relative):errors.append(aid+': final name map changed')
    tasks=table(trace)
    if not tasks or any(r['status'] not in ['COMPLETED','CACHED'] or r['exit'] not in ['0','-',''] for r in tasks):
        errors.append('Trace contains failed, ignored, incomplete or missing tasks')
    graphs=[r for r in tasks if ':CACTUS_PANGENOME ' in r['name']]
    expected_graphs=[c for c in new['cohorts'] if c['ready']]
    if len(graphs)!=len(expected_graphs) or any(r['status']!='CACHED' for r in graphs):errors.append('Accepted Cactus construction did not remain cached')
    for cohort in expected_graphs:
        taxid=str(cohort['taxid'])
        relative=Path('references')/taxid/'reference.json'
        a=json.loads((baseline/relative).read_text());b=json.loads((results/relative).read_text())
        for key in ['id','sha256','sequences','coordinate_stage','sequence_count','bases']:
            if a[key]!=b[key]:errors.append(taxid+': reference changed: '+key)
        for filename in ['clipping_audit.json','sequence_graph_audit.json','source_audit.json']:
            relative=Path('pangenome')/taxid/'representation'/filename
            a=json.loads((baseline/relative).read_text());b=json.loads((results/relative).read_text())
            if filename!='source_audit.json' and b['status']!='PASS':errors.append(taxid+': '+filename+' failed')
            if filename=='source_audit.json':
                for key in ['status','source_names_missing_from_both_graphs','chromosome_scale_paths_requiring_review','unknown_mapping_name_prefixes']:
                    if a.get(key)!=b.get(key):errors.append(taxid+': source-audit disposition changed: '+key)
    plan=table(results/'pairwise_alignments/comparison_plan.tsv')
    expected_pairs={frozenset([a['id'],b['id']]) for a,b in combinations(after.values(),2)
                    if a['sample']==b['sample'] and a['taxid']==b['taxid']}
    if len(plan)!=len(expected_pairs) or {frozenset([r['ref_id'],r['query_id']]) for r in plan}!=expected_pairs:
        errors.append('Within-sample comparison membership is incomplete or incorrect')
    for row in plan:
        prefix=row['ref_id']+'_vs_'+row['query_id']
        status=table(results/'pairwise_alignments'/(prefix+'.plot_status.tsv'))
        if len(status)!=2 or {r['plot'] for r in status}!={'dotplot','ribbon'} or any(
                r['status'] not in ['plotted','no_sequences_above_threshold','no_alignments_passing_plot_filters'] for r in status):
            errors.append(prefix+': figure status missing/invalid')
        for suffix in ['_dotplot.png','_riparian.png']:
            p=results/'pairwise_alignments'/(prefix+suffix)
            if not p.is_file() or p.stat().st_size==0:errors.append(prefix+': figure missing/empty')
    for required in ['assembly_run_summary.md','assembly_report.md','reports/assembly_qc_report.html',
                     'provenance/run_provenance.json','provenance/software_versions.tsv','provenance/dependency_catalog.tsv']:
        if not (results/required).is_file():errors.append('Missing required output: '+required)
    contexts=[task_context(r['workdir']) for r in tasks if r.get('workdir') not in [None,'','-']]
    env=environment_inventory(contexts)
    (out/'execution_environments.json').write_text(json.dumps(dict(tasks=contexts,environments=env),indent=2)+'\n')
    unresolved=sum(c['context_status']!='recorded' for c in contexts)
    notes.append(f'{unresolved} task wrappers unavailable; environment snapshots never substitute invented versions.')
    issues=table(results/'provenance/software_version_issues.tsv')
    notes.append(f'{len(issues)} explicit version-record issues; inspect provenance/software_version_issues.tsv.')
    result=dict(status='FAIL' if errors else 'PASS',finalized=len(after),comparisons=len(plan),
        completed=sum(r['status']=='COMPLETED' for r in tasks),cached=sum(r['status']=='CACHED' for r in tasks),
        graph_construction_cached=len(graphs),errors=errors,notes=notes)
    (out/'core_integration_check.json').write_text(json.dumps(result,indent=2)+'\n')
    text=['# Combined core integration check','',result['status'],'',
          f"{result['finalized']} final assemblies; {result['comparisons']} within-sample pairs; {result['completed']} completed / {result['cached']} cached tasks.",'',
          'Final input/name-map consistency, unchanged source-audit dispositions and Cactus cache reuse were checked.',
          'This check does not independently hash the complete graph or certify biological completeness.','']
    text+=['- '+item for item in errors+notes]
    (out/'core_integration_check.md').write_text('\n'.join(text)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['baseline','results','trace','out']:p.add_argument('--'+name,required=True)
    a=p.parse_args()
    try:result=audit(a.baseline,a.results,a.trace,a.out)
    except (OSError,ValueError,KeyError) as exc:
        result=dict(status='FAIL',errors=[str(exc)])
        Path(a.out).mkdir(parents=True,exist_ok=True)
        (Path(a.out)/'core_integration_check.json').write_text(json.dumps(result,indent=2)+'\n')
        (Path(a.out)/'core_integration_check.md').write_text('# Combined core integration check\n\nFAIL\n\n'+str(exc)+'\n')
    raise SystemExit(0 if result['status']=='PASS' else 1)
