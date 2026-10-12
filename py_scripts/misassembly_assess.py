#!/usr/bin/env python3
"""Measure cohort events, localize supported cut options, and retain review evidence."""
import argparse
import json
import math
import subprocess
import tempfile
import sys
import time
from collections import defaultdict
from pathlib import Path
from misassembly_core import (sha,table,write_table,peer_test,summarize_peers,support_grid,gap_spanning,contact_result,recommend,POLICY,matched_control_sites,read_depth,consensus_regions,coarse_tracks)
from chimera_controls import agp,scan_contacts,literal_gaps
from chimera_graph_evidence import measure_graph

_START=time.monotonic()


def progress(message):
    print(f'[assessment +{time.monotonic()-_START:.1f}s] {message}',file=sys.stderr,flush=True)


def draw_chromosome_track(axis,blocks,row,colors):
    """One collection per chromosome and peer, preserving every interval and gap."""
    groups=defaultdict(list)
    for block in blocks:
        if block.get('chrom'):
            groups[block['chrom']].append((block['lo']/1e6,(block['hi']-block['lo'])/1e6))
    for chrom,intervals in groups.items():
        axis.broken_barh(intervals,(row-.35,.7),facecolors=colors[chrom],edgecolors='none',linewidth=0)
    return len(groups)


def run(args,output=None):
    args=list(map(str,args))
    started=time.monotonic();progress('Starting '+ ' '.join(args[:2]))
    if output:
        with open(output,'w') as handle:subprocess.run(args,stdout=handle,check=True)
    else:subprocess.run(args,check=True)
    progress('Finished '+ ' '.join(args[:2])+f' in {time.monotonic()-started:.1f}s')


def sequence(fasta,chrom,lo,hi):
    text=subprocess.check_output(['samtools','faidx',str(fasta),'%s:%d-%d'%(chrom,lo+1,hi)],text=True)
    return ''.join(text.splitlines()[1:])


def peer_tests(event,lo,hi):
    return [peer_test(dict(id=t['peer'],sample=t['sample'],eligible=t['eligible']),t['blocks'],lo,hi) for t in event['peer_tracks']]


def assessment_gaps(fasta,current_gaps):
    """Recover N runs inherited inside components, without duplicating AGP gaps."""
    from bisect import bisect_right
    minimum=min([10]+[g['hi']-g['lo'] for g in current_gaps])
    literal=literal_gaps(fasta,minimum=minimum)
    by_scaffold=defaultdict(list)
    for gap in literal:by_scaffold[gap['scaffold']].append(gap)
    starts={sc:[g['lo'] for g in gaps] for sc,gaps in by_scaffold.items()}
    current_runs=set()
    for gap in current_gaps:
        index=bisect_right(starts.get(gap['scaffold'],[]),gap['lo'])-1
        enclosing=by_scaffold[gap['scaffold']][index] if index>=0 else None
        if enclosing is None or gap['hi']>enclosing['hi']:
            raise ValueError('Current AGP gap does not match literal N sequence')
        current_runs.add((enclosing['scaffold'],enclosing['lo'],enclosing['hi']))
    return [dict(g,origin='current_agp' if (g['scaffold'],g['lo'],g['hi']) in current_runs else 'inherited_sequence_gap') for g in literal]


def chromosome_colors(labels):
    """Stable chromosome colors across assemblies; no ten-chromosome wraparound."""
    import colorsys
    import hashlib
    result={}
    for label in labels:
        number=label.removeprefix('chr').split('_')[0]
        seed=int(number) if number.isdigit() else int(hashlib.sha256(label.encode()).hexdigest()[:8],16)
        result[label]=colorsys.hsv_to_rgb((seed*.61803398875)%1,.60,.78)
    return result


def contact_profile(matrix,resolution=100000,flank_bins=3):
    """Binned cross-boundary support; not the calibrated local gap assay."""
    import numpy as np
    rows=[]
    for i in range(flank_bins,len(matrix)-flank_bins+1):
        left=matrix[i-flank_bins:i,i-flank_bins:i];right=matrix[i:i+flank_bins,i:i+flank_bins]
        left_within=float((left.sum()+np.trace(left))/2)
        right_within=float((right.sum()+np.trace(right))/2)
        cross=float(matrix[i-flank_bins:i,i:i+flank_bins].sum())
        ratio=cross/math.sqrt(left_within*right_within) if min(left_within,right_within)>0 else None
        rows.append(dict(position=i*resolution,ratio=ratio,cross=cross,left_within=left_within,right_within=right_within))
    return rows


def shade_chromosomes(axes,regions,colors):
    """Only shade assigned regions; do not fill ambiguous/unmapped intervals."""
    for axis in axes:
        for region in regions:
            axis.axvspan(region['lo']/1e6,region['hi']/1e6,facecolor=colors[region['chrom']],alpha=.10,zorder=0)


def plot_event(event,out):
    started=time.monotonic();progress('Starting candidate plot '+event['id'])
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    grid=event.get('hifi',{});tracks=event['peer_tracks']
    regions=event.get('chromosome_regions',[])
    labels=sorted({b['chrom'] for t in tracks for b in t['blocks'] if b.get('chrom')}|{r['chrom'] for r in regions})
    colors=chromosome_colors(labels)
    fig,axes=plt.subplots(5,1,figsize=(12,max(12,8+.32*len(tracks))),sharex=True,
                         gridspec_kw={'height_ratios':[max(2,len(tracks)*.25),1.5,1,1.3,1]})
    collections=sum(draw_chromosome_track(axes[0],t['blocks'],i,colors) for i,t in enumerate(tracks))
    progress(f"{event['id']}: {sum(sum(bool(b.get('chrom')) for b in t['blocks']) for t in tracks)} intervals in {collections} drawing collections")
    axes[0].set_yticks(range(len(tracks)),[t['peer'] for t in tracks]);axes[0].invert_yaxis()
    fig.legend(handles=[Patch(color=colors[c],label=c) for c in labels],ncol=min(8,max(1,len(labels))),
               loc='upper center',bbox_to_anchor=(.5,.965),fontsize=9,frameon=False)
    axes[0].set_title('Chromosome assignments; blank regions have no usable assignment',fontsize=10)
    for name,color,label,style in [('hifi','#2166ac','MAPQ ≥20','-'),('hifi_strict','#b35806','MAPQ ≥30','-'),
                                  ('hifi_all_primary','gray','All primary placements (includes ambiguous)',':')]:
        probes=event.get(name+'_profile',event.get(name,{})).get('probes',[])
        if probes:axes[1].plot([p['position']/1e6 for p in probes],[p['molecules'] for p in probes],color=color,label=label,ls=style,lw=1)
    axes[1].axhline(2,color='gray',ls=':');axes[1].set_ylabel('HiFi spanning molecules')
    if grid.get('probes'):axes[1].legend(fontsize=8)
    else:axes[1].text(.5,.5,'HiFi evidence unavailable',ha='center',transform=axes[1].transAxes)
    for name,color,label in [('depth_all_primary','gray','All primary placements'),('depth','#238b45','MAPQ ≥20')]:
        if event.get(name):axes[2].plot([(p['start']+p['end'])/2e6 for p in event[name]],[p['depth'] for p in event[name]],color=color,label=label,lw=1)
    if event.get('depth'):axes[2].legend(fontsize=8)
    axes[2].set_ylabel('Mean HiFi depth')
    profiles=event.get('contact_profiles',{})
    for library,rows in sorted(profiles.items()):
        axes[3].plot([r['position']/1e6 for r in rows],[float('nan') if r['ratio'] is None else r['ratio'] for r in rows],label=library,lw=1)
    axes[3].set_ylabel('Hi-C contact support\nCross / √(within L × R)')
    axes[3].set_title('Binned Hi-C support · 100 kb bins, 300 kb flanks · libraries shown separately',fontsize=10)
    if profiles:axes[3].legend(fontsize=8)
    else:axes[3].text(.5,.5,'Hi-C profiles unavailable',ha='center',transform=axes[3].transAxes)
    telomeres=event.get('telomeres',[])
    if telomeres:axes[4].plot([r['start']/1e6 for r in telomeres],[r['count'] for r in telomeres],color='#762a83',lw=.8)
    else:axes[4].text(.5,.5,'Telomere measurements unavailable',ha='center',transform=axes[4].transAxes)
    axes[4].set_ylabel('Telomere motif count\nper 1 kb');axes[4].set_xlabel('Original scaffold coordinate (Mb)')
    shade_chromosomes(axes[1:],regions,colors)
    for ax in axes:
        ax.set_facecolor('white')
        ax.axvspan(event['lo']/1e6,event['hi']/1e6,color='orange',alpha=.12,zorder=0)
        for cut in event['cut_options']:ax.axvline(cut['cut_bp']/1e6,color='#b2182b',ls='--')
        ax.grid(alpha=.15)
    for ax in axes[1:]:ax.set_ylim(bottom=0)
    axes[-1].set_xlim(max(0,event['lo']-500000)/1e6,(event['hi']+500000)/1e6)
    fig.suptitle(event['id']+' · '+event['scaffold']+' · '+event['left']+' → '+event['right'],y=.995)
    fig.tight_layout(rect=(0,0,1,.92));fig.savefig(out/'plots'/(event['id']+'.png'),dpi=150,facecolor='white');plt.close(fig)
    progress(f"Finished candidate plot {event['id']} in {time.monotonic()-started:.1f}s")


def plot_scaffold(chrom,matrix,events,out,resolution=100000):
    """The standalone scaffold figure contains only the pooled contact matrix."""
    started=time.monotonic();progress('Starting contact map '+chrom)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    fig,axis=plt.subplots(figsize=(10,9))
    if matrix is not None:
        extent=(0,matrix.shape[0]*resolution/1e6,matrix.shape[0]*resolution/1e6,0)
        image=axis.imshow(np.log1p(matrix),extent=extent,cmap='magma',aspect='equal')
        fig.colorbar(image,ax=axis,label='log(1 + Hi-C pairs); libraries pooled')
        for event in events:
            if event['scaffold']!=chrom:continue
            position=(event['lo']+event['hi'])/2e6
            axis.axvline(position,color='cyan',ls='--',lw=.8);axis.axhline(position,color='cyan',ls='--',lw=.8)
            axis.annotate(event['id'],(position,1),xycoords=('data','axes fraction'),fontsize=9,rotation=90,va='bottom')
    else:axis.text(.5,.5,'Hi-C mappings unavailable',ha='center',transform=axis.transAxes)
    axis.set_title(chrom+' · contact map',pad=28)
    axis.set_ylabel('Scaffold position (Mb)');axis.set_xlabel('Scaffold position (Mb)')
    fig.tight_layout();fig.savefig(out/'plots'/(chrom+'.contacts.png'),dpi=150,facecolor='white');plt.close(fig)
    progress(f'Finished contact map {chrom} in {time.monotonic()-started:.1f}s')


def main():
    p=argparse.ArgumentParser()
    for key in ('fasta','discovery','mapping','out'):p.add_argument('--'+key,required=True)
    for key in ('agp','pairs','source','libraries','graph'):p.add_argument('--'+key)
    p.add_argument('--motif',default='CCCTAA');p.add_argument('--threads',type=int,default=8);p.add_argument('--plots',choices=('true','false'),default='true')
    a=p.parse_args();out=Path(a.out);out.mkdir()
    data=json.loads(Path(a.discovery).read_text());mapping=json.loads(Path(a.mapping).read_text())
    if mapping.get('bam'):
        mapping['bam']=str((Path(a.mapping).resolve().parent/mapping['bam']).resolve())
        mapping['index']=str((Path(a.mapping).resolve().parent/mapping['index']).resolve())
    if data['assessment_sha256']!=sha(a.fasta) or mapping['assessment_sha256']!=data['assessment_sha256']:raise ValueError('Assessment sources do not match')
    events=data['events'];data['mapping']=mapping
    progress(f"{data['assembly']}: {len(events)} candidate events")
    data['assessment_fasta']=str(Path(a.fasta).resolve())
    if events and a.plots=='true':(out/'plots').mkdir()
    data['tools']={tool:subprocess.check_output([tool,'--version'],text=True).splitlines()[0] for tool in ('samtools','minimap2','tidk')}
    support_rows=[];contact_rows=[];block_rows=[];telomere_rows=[]
    with tempfile.TemporaryDirectory(prefix='misassembly-',dir='.') as temp:
        work=Path(temp).resolve();fa=work/'assessment.fa';fa.symlink_to(Path(a.fasta).resolve());run(['samtools','faidx',fa])
        sizes={f[0]:int(f[1]) for f in (l.split('\t') for l in Path(str(fa)+'.fai').read_text().splitlines())}
        gaps=[];placements={}
        if a.agp:
            gaps,placements,lengths=agp(a.agp)
            if lengths!=sizes:raise ValueError('Current AGP and FASTA dictionaries differ')
        gaps=assessment_gaps(a.fasta,gaps) if events else []
        intervals={};queries=work/'queries.fa'
        with queries.open('w') as handle:
            for event in events:
                lo,hi=event['lo'],event['hi'];sc=event['scaffold'];key=event['id']
                if not 0<=lo<=hi<=sizes[sc]:raise ValueError('Discovery transition outside FASTA')
                event['peer_tests']=peer_tests(event,lo,hi);event['peers']=summarize_peers(event['peer_tests'],data['sample'],(event['left'],event['right']))
                for trial in event['peer_tests']:trial['sister']=trial['sample']==data['sample']
                event['cut_options']=[]
                for gap in gaps:
                    if gap['scaffold']!=sc or gap['lo']<lo-250000 or gap['hi']>hi+250000:continue
                    if set(sequence(fa,sc,gap['lo'],gap['hi']).upper())!={'N'}:raise ValueError('AGP gap is not literal N sequence')
                    option=dict(scaffold=sc,lo=gap['lo'],hi=gap['hi'],cut_bp=gap['hi'],origin=gap['origin'],id=key+'G%02d'%(len(event['cut_options'])+1))
                    option['peer_tests']=peer_tests(event,option['lo'],option['hi']);option['peers']=summarize_peers(option['peer_tests'],data['sample'],(event['left'],event['right']))
                    # A nearby gap is an alternative hypothesis. It must not jump
                    # over a reversal or a third chromosome to become a recommendation.
                    bridge_lo=min(lo,option['lo']);bridge_hi=max(hi,option['hi'])
                    bridge=[]
                    for t in event['peer_tracks']:
                        if not t['eligible'] or t['sample']==data['sample']:continue
                        blocks=[b for b in coarse_tracks(t['blocks'],1000) if b['lo']<bridge_hi and b['hi']>bridge_lo]
                        bases=sum(max(0,min(bridge_hi,b['hi'])-max(bridge_lo,b['lo']))*b['aligned_bp']/max(1,b['hi']-b['lo']) for b in blocks)
                        order=[]
                        for b in blocks:
                            if not order or order[-1]!=b['chrom']:order.append(b['chrom'])
                        valid=order in ([event['left']],[event['right']],[event['left'],event['right']])
                        bridge.append(dict(peer=t['peer'],sample=t['sample'],assigned_bp=bases,order=order,consistent=valid,observable=bases>=.5*max(1,bridge_hi-bridge_lo)))
                    option['bridge']=bridge
                    option['bridge_consistent']=(lo<=option['lo'] and option['hi']<=hi or len({b['sample'] for b in bridge if b['observable'] and b['consistent']})>=2) and not any(b['observable'] and not b['consistent'] for b in bridge)
                    event['cut_options'].append(option)
                for target in [event]+event['cut_options']:
                    flank=500000 if target is event else 100000
                    target['start']=max(0,target['lo']-flank);target['end']=min(sizes[sc],target['hi']+flank)
                    intervals[target['id']]=dict(target,scaffold=sc)
                    handle.write('>'+target['id']+'\n'+sequence(fa,sc,target['start'],target['end'])+'\n')
        progress('Building matched controls and chromosome context')
        controls,control_links=matched_control_sites(intervals,gaps,sizes)
        # Gap controls require independent same-chromosome context. Read-supported
        # continuous controls require both that context and a measured read chain.
        for control in controls.values():
            tests=[peer_test(dict(id=peer['peer'],sample=peer['sample'],eligible=peer['eligible']),peer['scaffolds'].get(control['scaffold'],[]),control['lo'],control['hi']) for peer in data.get('coarse_context',[])]
            control['context_individuals']=summarize_peers(tests,data['sample'],('__none__','__none__'))['continuous']
        for event in events:
            progress('Starting HiFi assays '+event['id'])
            if mapping['status'] in ('mapped','reused'):
                for target in [event]+event['cut_options']:
                    sam=work/(target['id']+'.sam');run(['samtools','view','-h',mapping['bam'],'%s:%d-%d'%(event['scaffold'],target['start']+1,target['end'])],sam)
                    if target is event:
                        target['depth']=read_depth(sam,target['start'],target['end']);target['depth_all_primary']=read_depth(sam,target['start'],target['end'],mapq=0)
                        target['hifi']=support_grid(sam,event['lo'],event['hi'],mapq=20);target['hifi_strict']=support_grid(sam,event['lo'],event['hi'],mapq=30);target['hifi_all_primary']=support_grid(sam,event['lo'],event['hi'],mapq=0)
                        for tier in ('hifi','hifi_strict','hifi_all_primary'):
                            support_rows.extend(dict(event=event['id'],mapq=target[tier]['mapq'],**site) for site in target[tier]['probes'])
                            target[tier+'_profile']=support_grid(sam,target['start'],target['end'],mapq=target[tier]['mapq'])
                    else:
                        target['spanning']=gap_spanning(sam,target['lo'],target['hi'])
                        target['flank_molecules']=[support_grid(sam,target['lo']-1000,target['lo']-1000)['minimum'],support_grid(sam,target['hi']+1000,target['hi']+1000)['minimum']]
                # Failed reads cannot justify a cut; direct evidence is explicitly unavailable.
            else:event['hifi']=dict(state='unavailable')
            progress('Finished HiFi assays '+event['id'])
        progress('Starting read-supported control validation')
        usable_controls={}
        for key,control in controls.items():
            if len(control['context_individuals'])<2:continue
            if control['role']=='continuous_control':
                if mapping['status'] not in ('mapped','reused'):continue
                sam=work/(key+'.sam');run(['samtools','view','-h',mapping['bam'],'%s:%d-%d'%(control['scaffold'],control['lo']-2000+1,control['hi']+2000)],sam)
                if support_grid(sam,control['lo'],control['hi'])['minimum']<2:continue
            usable_controls[key]=control
        matrices={};library_matrices={};profiles={}
        if a.pairs and events:
            if not all((a.source,a.libraries,a.agp)):raise ValueError('Hi-C evidence requires its source FASTA, last AGP and library identities')
            source=work/'source.fa';source.symlink_to(Path(a.source).resolve());run(['samtools','faidx',source])
            source_sizes={f[0]:int(f[1]) for f in (l.split('\t') for l in Path(str(source)+'.fai').read_text().splitlines())}
            if any(name not in source_sizes or v[4]>source_sizes[name] for name,values in placements.items() for v in values):raise ValueError('AGP components outside Hi-C source FASTA')
            libraries=table(a.libraries)
            if a.plots=='true':
                import numpy as np
                matrices={sc:np.zeros((math.ceil(sizes[sc]/100000),)*2,dtype=np.uint32) for sc in {e['scaffold'] for e in events}}
                library_matrices={library:{sc:np.zeros_like(matrix) for sc,matrix in matrices.items()} for library in sorted({r['library_id'] for r in libraries})}
            def capture(library,x,y):
                if x[0]==y[0] and x[0] in matrices:
                    i,j=x[1]//100000,y[1]//100000
                    for matrix in (matrices[x[0]],library_matrices[library][x[0]]):
                        matrix[i,j]+=1
                        if i!=j:matrix[j,i]+=1
            progress(f'Starting Hi-C scan: {len(intervals)} locations, {len(usable_controls)} controls')
            counts,totals,audit=scan_contacts(a.pairs,placements,dict(intervals,**usable_controls),libraries,on_pair=capture)
            progress('Finished Hi-C scan')
            profiles={sc:{library:contact_profile(values[sc]) for library,values in library_matrices.items()} for sc in matrices}
            data['hic_audit']=dict(audit=audit,libraries=totals,source_sha256=sha(a.source),agp_sha256=sha(a.agp))
            for event in events:
                for target in [event]+event['cut_options']:
                    target['contacts']=[]
                    population='continuous_control' if target is event else 'gap_control'
                    for library in sorted({r['library_id'] for r in libraries}):
                        result=contact_result(counts[library,target['id']],[dict(counts[library,key],id=key) for key in control_links[target['id']] if key in usable_controls])
                        result.update(library=library,control_population=population);target['contacts'].append(result)
                        contact_rows.append(dict(event=event['id'],location=target['id'],library=library,**result['raw_counts'],ratio=result['ratio'],matched_controls=result['matched_controls'],calibrated=result['calibrated'],contact_loss=result['loss']))
        if a.graph and events:
            graph=measure_graph(a.graph,work,queries,intervals,run,a.threads)
            for event in events:
                for target in [event]+event['cut_options']:target.update(graph.get(target['id'],{}))
        if events and a.plots=='true':
            all_fa=work/'scaffolds.fa'
            with all_fa.open('w') as handle:
                for sc in sorted({e['scaffold'] for e in events}):handle.write('>'+sc+'\n'+sequence(fa,sc,0,sizes[sc])+'\n')
            run(['tidk','search','--string',a.motif,'--window','1000','--output','scaffolds','--dir',work/'tidk',all_fa])
            for path in sorted((work/'tidk').glob('*.tsv')):
                for r in table(path):
                    if 'id' in r and 'window' in r:telomere_rows.append(dict(scaffold=r['id'],start=max(0,int(r['window'])-1000),count=int(r['forward_repeat_number'])+int(r['reverse_repeat_number'])))
        for event in events:
            event['recommendation'],event['reason']=recommend(event)
            for t in event['peer_tracks']:block_rows.extend(dict(event=event['id'],peer=t['peer'],**b) for b in t['blocks'])
            if a.plots=='true':
                event['chromosome_regions']=consensus_regions(event['peer_tracks'],data['sample'])
                event['contact_profiles']={library:[r for r in rows if event['lo']-500000<=r['position']<=event['hi']+500000] for library,rows in profiles.get(event['scaffold'],{}).items()}
                event['telomeres']=[r for r in telomere_rows if r['scaffold']==event['scaffold'] and event['lo']-500000<=r['start']<=event['hi']+500000]
                plot_event(event,out)
        if a.plots=='true':
            for sc in sorted({e['scaffold'] for e in events}):plot_scaffold(sc,matrices.get(sc),events,out)
        data['control_registry']=usable_controls
    # Publish only evidence used by the report and decisions. No SAMs, duplicate FASTAs,
    # speculative IGV sessions, recovered historical joins or diagnostic HTML files.
    data.pop('coarse_context',None)
    progress('Writing evidence JSON and tables')
    (out/'evidence.json').write_text(json.dumps(data,indent=2),encoding='utf-8')
    if support_rows:write_table(out/'read_support.tsv',support_rows,['event','mapq','position','molecules'])
    if contact_rows:write_table(out/'contact_support.tsv',contact_rows,['event','location','library','cross','left_within','right_within','ratio','matched_controls','calibrated','contact_loss'])
    if block_rows:write_table(out/'chromosome_blocks.tsv',block_rows,['event','peer','lo','hi','chrom','status'])
    if telomere_rows:write_table(out/'telomeres.tsv',telomere_rows,['scaffold','start','count'])
    progress('Assessment complete')


if __name__=='__main__':main()
