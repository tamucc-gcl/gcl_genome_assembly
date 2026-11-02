"""Execute the real R summary renderer with graph reporting enabled and disabled."""
import argparse
import csv
import json
import os
from pathlib import Path
import re
import subprocess


def run(rscript,root,name,graph,actual=None):
    out=root/name;out.mkdir(parents=True,exist_ok=True)
    identities=['fixture_hap1','fixture_hap2']
    if actual:
        identities=[r['id'] for r in json.loads((actual/'assembly_eligibility.json').read_text())['assemblies'] if r['completion']=='finalized']
    manifest=out/'manifest.tsv'
    with manifest.open('w',newline='') as h:
        w=csv.writer(h,delimiter='\t');w.writerow(['type','id','id2','filename','subdir'])
        for aid in identities:w.writerow(['assembly',aid,'.',aid+'.fasta','assembly/final'])
    qc=actual/'qc/assembly/assembly_qc_metrics.csv' if actual else out/'qc.csv'
    if not actual:
        with qc.open('w',newline='') as h:
            w=csv.writer(h);w.writerow(['sample_id','stage','analysis','metric','hap1','hap2','both'])
            for metric,value in [('Total length (>= 0 bp)',1000000),('N50',500000),('# contigs (>= 0 bp)',2)]:
                w.writerow(['fixture','final','quast',metric,value,value,'NA'])
    versions=out/'software_versions.tsv'
    versions.write_text('tool\tversion\tprocess\tsource\nminimap2\t2.28\tONE\tfixture\nminimap2\t2.30\tTWO\tfixture\n')
    fragment=out/'pangenome.md';fragment.write_text('## Pangenome\n\nGRAPH_REPORT_TEST_CONTENT\n')
    command=[rscript,str(Path(__file__).resolve().parents[1]/'r_scripts/generate_summary_report.R'),
             '--manifest',str(manifest),'--compiled_qc',str(qc),'--telomere_summary','NO_TELOMERES',
             '--pairwise_summary','NO_PAIRWISE','--mito_stats','NO_MITO_STATS','--teloclip_stats','NO_TELOCLIP',
             '--pangenome_report',str(fragment) if graph else 'NO_PANGENOME','--versions',str(versions),
             '--output',str(out/'report.md')]
    env=os.environ.copy();env['TZ']='UTC';env['PATH']=str(Path(rscript).parent)+os.pathsep+env['PATH']
    result=subprocess.run(command,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,env=env)
    (out/'renderer.log').write_text(result.stdout)
    if result.returncode:raise AssertionError(name+': '+result.stdout[-3500:])
    text=(out/'report.md').read_text()
    assert ('GRAPH_REPORT_TEST_CONTENT' in text)==graph
    assert bool(re.search(r'^## \d+\. Pangenome$',text,re.M))==graph
    assert '| minimap2 | 2.28 |' in text and '| minimap2 | 2.30 |' in text
    print('PASS',name,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('rscript');p.add_argument('--actual-root',type=Path)
    a=p.parse_args();root=Path(__file__).resolve().parents[1]/'comparisons/summary-renderer-validation';root.mkdir(exist_ok=True)
    run(a.rscript,root,'graph_enabled',True)
    run(a.rscript,root,'graph_disabled',False)
    if a.actual_root:run(a.rscript,root,'real_qc_replay',True,a.actual_root.resolve())
