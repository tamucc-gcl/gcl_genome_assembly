#!/usr/bin/env python3
"""Production junction evidence; consumes staged inputs, never discovers old cache paths."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
from misassembly_junction import labels,hifi,hic,windows,table,save,sha,COMMANDS,credible_bridge,refine_advice


def attach_summaries(data, root):
    bridge=table(root/'hifi'/'competitive_bridges.tsv') if (root/'hifi'/'competitive_bridges.tsv').exists() else []
    reach={r['site']:r for r in table(root/'hifi'/'unique_anchor_reach.tsv')} if (root/'hifi'/'unique_anchor_reach.tsv').exists() else {}
    islands=table(root/'labels'/'opposite_side_islands.tsv')
    contacts=table(root/'hic'/'unique_arm_contacts.tsv') if (root/'hic'/'unique_arm_contacts.tsv').exists() else []
    for event in data['events']:
        sites=[]
        for site in event['junction_sites']:
            rows=[r for r in bridge if r['site']==site['id']]
            credible=[r for r in rows if credible_bridge(r)]
            unique=[r for r in credible if min(int(r['unique_left_bp']),int(r['unique_right_bp']))>=1000]
            rr=reach.get(site['id'],{})
            unique_islands=[r for r in islands if r['site']==site['id'] and r['eligible']=='True' and r['sample']!=data['sample']]
            # Gap-spanning chains use a separate, N-aware assay in the base assessment.
            # This follow-up tests uninterrupted sequence and must not turn an N gap
            # into apparent absence of physical molecule support.
            applicable=bool(rr) and site['kind']!='gap_option'
            counts=len(credible) if applicable else None
            sites.append(dict(id=site['id'],point=site['point'],credible_bridges=counts,
                unique_arm_bridges=len(unique) if applicable else None,
                minimum_anchor_span_bp=rr.get('minimum_unique_anchor_span_bp') or None,
                local_read_p95=rr.get('local_read_p95') or None,
                local_reads_long_enough=rr.get('local_reads_long_enough') or None,
                largest_independent_opposite_island_bp=max((int(r['largest_unique_block_bp']) for r in unique_islands),default=0),
                contacts=[r for r in contacts if r['site']==site['id'] and r['library']!='POOLED' and r['anchor_bp']=='100000']))
        event['junction_evidence']=dict(sites=sites,
            interpretation='Positive unique-arm bridges support continuity. Unreachable anchors and missing controls cannot establish a deficit. Small opposite-side islands are evidence to inspect, not a veto on a manual cut.',
            uniqueness_scope='Unique among retained emitted assembly alignments; not proven single-copy.',
            hifi_status='measured' if reach else 'unavailable')
        event['junction_evidence']['plots']=[p.name for p in root.glob(event['id']+'.*.png')]
        event['junction_evidence']['artifacts']=[str(p.relative_to(root)) for p in root.rglob('*') if p.is_file()]
        by_site={r['id']:r for r in sites}
        for target in [event]+event.get('alternatives',[])+event.get('cut_options',[]):
            target['junction_support']=by_site.get(target['id'] if target is not event else event['id']+'T',{})
        # Re-evaluate site advice with direct positive molecule counterevidence.
        if 'peers' in event:refine_advice(event)
        event.pop('junction_sites',None)


def make_plots(job, root):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    for event in job['events']:
        for record in json.loads((root/'hic'/'summary.json').read_text())['maps']:
            if record['event']!=event['event'] or record['library']=='POOLED' or record['resolution'] not in (10000,25000):continue
            with np.load(root/'hic'/record['file']) as archive:
                matrix=archive['matrix'];lo=int(archive['lo']);res=int(archive['resolution'])
            fig,ax=plt.subplots(figsize=(8,7),facecolor='white')
            hi=lo+len(matrix)*res
            im=ax.imshow(matrix+1,origin='lower',extent=[lo/1e6,hi/1e6,lo/1e6,hi/1e6],
                cmap='magma',norm=LogNorm(vmin=1,vmax=max(2,int(matrix.max())+1)),interpolation='nearest')
            # Names are on axes, keeping many nearby markers legible.
            for site in event['sites']:
                point=site['point']/1e6
                ax.axvline(point,color='cyan',lw=.5,alpha=.5);ax.axhline(point,color='cyan',lw=.5,alpha=.5)
            ax.set(xlabel='Original scaffold coordinate (Mb)',ylabel='Original scaffold coordinate (Mb)',
                title=f"{job['assembly']} / {event['event']} / {record['library']} / {res//1000} kb")
            fig.colorbar(im,ax=ax,label='Retained pairs + 1 (log scale)')
            fig.tight_layout();fig.savefig(root/f"{event['event']}.{record['library']}.{res}.png",dpi=150,facecolor='white');plt.close(fig)


def main():
    p=argparse.ArgumentParser()
    for key in ('packet','fasta','mapping','catalog','alignments','competitive','out'):p.add_argument('--'+key,required=True)
    for key in ('agp','pairs','source','libraries'):p.add_argument('--'+key)
    p.add_argument('--threads',type=int,default=8);p.add_argument('--plots',choices=('true','false'),default='true')
    a=p.parse_args();out=Path(a.out);shutil.copytree(a.packet,out)
    evidence=out/'evidence.json';data=json.loads(evidence.read_text());root=out/'junction';root.mkdir()
    if data['assessment_sha256']!=sha(a.fasta):raise ValueError('Junction FASTA checksum mismatch')
    data['junction_policy']='junction-evidence-v1'
    if not data['events']:
        save(root/'status.json',dict(status='not_needed',reason='No detected events'));save(evidence,data);return
    mapping=json.loads((Path(a.mapping)/'mapping.json').read_text())
    if mapping['assessment_sha256']!=data['assessment_sha256']:raise ValueError('Junction mapping checksum mismatch')
    sizes={};name=None
    with open(a.fasta) as f:
        for line in f:
            if line.startswith('>'):name=line[1:].split()[0];sizes[name]=0
            else:sizes[name]+=len(line.strip())
    events=[]
    for e in data['events']:
        sites,lo,hi=windows(e,sizes[e['scaffold']]);e['junction_sites']=sites
        events.append(dict(event=e['id'],sites=sites,scaffold=e['scaffold'],left=e['left'],right=e['right'],lo=lo,hi=hi))
    pairs=json.loads(a.alignments)
    expected={(r['a'],r['b'],r['sha256']) for r in data['alignments']}
    if {(r['a'],r['b'],sha(r['path'])) for r in pairs}!=expected:raise ValueError('Junction cohort alignment checksums differ')
    job=dict(assembly=data['assembly'],sample=data['sample'],events=events,sizes=sizes,
        catalog=a.catalog,evidence=str(evidence),alignments=pairs,paths=dict(fasta=a.fasta),competitive_fastas=json.loads(a.competitive),resolutions=(5000,10000,25000),include_pooled=False)
    for key in ('agp','pairs','source','libraries'):
        if getattr(a,key):job['paths'][key]=getattr(a,key)
    if a.pairs:
        for key in ('agp','source'):
            if sha(job['paths'][key])!=data['hic_audit'][key+'_sha256']:raise ValueError('Junction Hi-C source mismatch')
        job['hic_expected']=data['hic_audit']['audit']['input_pairs']
    for branch in ('labels','hifi','hic'):(root/branch).mkdir()
    labels(job,root/'labels',a.threads)
    mask=root/'labels'/'masks.json'
    if mapping['status'] in ('mapped','reused'):
        job['bam']=str((Path(a.mapping)/mapping['bam']).resolve())
        job['index']=str((Path(a.mapping)/mapping['index']).resolve())
        subprocess.run(['samtools','quickcheck',job['bam']],check=True)
        hifi(job,root/'hifi',a.threads,mask)
    else:save(root/'hifi'/'summary.json',dict(status='unavailable'))
    if a.pairs:
        hic(job,root/'hic',a.threads,mask)
        if a.plots=='true':make_plots(job,root)
    else:save(root/'hic'/'summary.json',dict(status='unavailable',maps=[]))
    attach_summaries(data,root)
    if a.plots=='true':
        from misassembly_plot import plot_event
        for event in data['events']:
            if 'peers' in event:
                plot_event(event,out)
                plot_event(event,out,detail=True)
    save(root/'status.json',dict(status='complete',coordinate_stage='pre_finishing',assessment_sha256=data['assessment_sha256'],commands=COMMANDS))
    save(evidence,data)


if __name__=='__main__':main()
