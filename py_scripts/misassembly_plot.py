"""Render review figures from retained measurements, without re-running assays."""
import sys,time
from collections import defaultdict
_START=time.monotonic()
def progress(message):
    print(f"[plot +{time.monotonic()-_START:.1f}s] {message}",file=sys.stderr,flush=True)

def draw_chromosome_track(axis,blocks,row,colors):
    """One collection per chromosome and peer, preserving every interval and gap."""
    groups=defaultdict(list)
    for block in blocks:
        if block.get('chrom'):
            groups[block['chrom']].append((block['lo']/1e6,(block['hi']-block['lo'])/1e6))
    for chrom,intervals in groups.items():
        axis.broken_barh(intervals,(row-.35,.7),facecolors=colors[chrom],edgecolors='none',linewidth=0)
    return len(groups)


def chromosome_colors(labels,focal=()):
    """Contrast the focal pair consistently across its two haplotype orientations."""
    import colorsys
    import hashlib
    result={}
    for label in labels:
        number=label.removeprefix('chr').split('_')[0]
        seed=int(number) if number.isdigit() else int(hashlib.sha256(label.encode()).hexdigest()[:8],16)
        result[label]=colorsys.hsv_to_rgb((seed*.61803398875)%1,.25 if focal else .60,.65 if focal else .78)
    def order(label):
        number=label.removeprefix('chr').split('_')[0]
        return (0,int(number)) if number.isdigit() else (1,label)
    for label,color in zip(sorted(set(focal),key=order),('#0072B2','#D55E00')):
        result[label]=color
    return result


def shade_chromosomes(axes,regions,colors):
    """Only shade assigned regions; do not fill ambiguous/unmapped intervals."""
    for axis in axes:
        for region in regions:
            axis.axvspan(region['lo']/1e6,region['hi']/1e6,facecolor=colors[region['chrom']],alpha=.10,zorder=0)


def detail_window(event,selected_cuts=()):
    candidates=[t for t in event.get('alternatives',[]) if t.get('site_verdict')!='EXCLUDE']
    selected=[int(r['cut_bp']) for r in selected_cuts]
    focus=event.get('recommended_cut') or event.get('preferred_region') or (candidates[0] if candidates else event)
    lo=hi=selected[0] if selected else focus.get('cut_bp') or focus.get('probe_bp') or (focus['lo']+focus['hi'])//2
    if not selected and not focus.get('cut_bp') and focus['hi']-focus['lo']<=200000:lo,hi=focus['lo'],focus['hi']
    return max(0,lo-100000),hi+100000


def plot_event(event,out,detail=False,selected_cuts=()):
    started=time.monotonic();progress('Starting candidate plot '+event['id'])
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    grid=event.get('hifi',{});tracks=event['peer_tracks']
    regions=event.get('chromosome_regions',[])
    labels=sorted({b['chrom'] for t in tracks for b in t['blocks'] if b.get('chrom')}|{r['chrom'] for r in regions})
    focal=sorted(set((event['left'],event['right'])),key=lambda c: (not c.removeprefix('chr').isdigit(),int(c.removeprefix('chr')) if c.removeprefix('chr').isdigit() else c))
    labels=sorted(set(labels)|set(focal))
    colors=chromosome_colors(labels,focal)
    fig,axes=plt.subplots(5,1,figsize=(12,max(12,8+.32*len(tracks))),sharex=True,
                         gridspec_kw={'height_ratios':[max(2,len(tracks)*.25),1.5,1,1.3,1]})
    collections=sum(draw_chromosome_track(axes[0],t['blocks'],i,colors) for i,t in enumerate(tracks))
    progress(f"{event['id']}: {sum(sum(bool(b.get('chrom')) for b in t['blocks']) for t in tracks)} intervals in {collections} drawing collections")
    axes[0].set_yticks(range(len(tracks)),[t['peer'] for t in tracks]);axes[0].invert_yaxis()
    fig.legend(handles=[Patch(color=colors[c],label=c) for c in focal],ncol=len(focal),
               title='Chromosomes in this transition',loc='upper center',bbox_to_anchor=(.5,.973),fontsize=10,frameon=False)
    other=[c for c in labels if c not in focal]
    if other:
        fig.legend(handles=[Patch(color=colors[c],label=c) for c in other],ncol=min(8,len(other)),
                   title='Other chromosome assignments',loc='upper center',bbox_to_anchor=(.5,.927),fontsize=8,title_fontsize=8,frameon=False)
    axes[0].set_title('Chromosome assignments; blank regions have no usable assignment',fontsize=10,pad=24)
    ambiguous_axes=[]
    for name,color,label,style in [('hifi','#2166ac','MAPQ ≥20','-'),('hifi_strict','#b35806','MAPQ ≥30','-')]:
        probes=event.get(name+'_profile',event.get(name,{})).get('probes',[])
        if probes:axes[1].plot([p['position']/1e6 for p in probes],[p['molecules'] for p in probes],color=color,label=label,ls=style,lw=1)
    axes[1].axhline(2,color='gray',ls=':');axes[1].set_ylabel('HiFi spanning molecules')
    ambiguous=event.get('hifi_all_primary_profile',event.get('hifi_all_primary',{})).get('probes',[])
    if ambiguous:
        extra=axes[1].twinx();extra.plot([p['position']/1e6 for p in ambiguous],[p['molecules'] for p in ambiguous],color='gray',alpha=.6,ls=':',label='All primary placements · right axis')
        extra.set_ylabel('All primary molecules',color='gray',fontsize=8);extra.tick_params(axis='y',colors='gray',labelsize=8);extra.set_ylim(bottom=0);ambiguous_axes.append(extra)
    if grid.get('probes'):
        handles,legend=axes[1].get_legend_handles_labels()
        for extra in ambiguous_axes:
            h,l=extra.get_legend_handles_labels();handles+=h;legend+=l
        axes[1].legend(handles,legend,fontsize=8)
    else:axes[1].text(.5,.5,'HiFi evidence unavailable',ha='center',transform=axes[1].transAxes)
    if event.get('depth'):
        axes[2].plot([(p['start']+p['end'])/2e6 for p in event['depth']],[p['depth'] for p in event['depth']],color='#238b45',label='MAPQ ≥20',lw=1)
    if event.get('depth_all_primary'):
        extra=axes[2].twinx();extra.plot([(p['start']+p['end'])/2e6 for p in event['depth_all_primary']],[p['depth'] for p in event['depth_all_primary']],color='gray',alpha=.6,ls=':',label='All primary depth · right axis')
        extra.set_ylabel('All primary depth',color='gray',fontsize=8);extra.tick_params(axis='y',colors='gray',labelsize=8);extra.set_ylim(bottom=0)
        handles,legend=axes[2].get_legend_handles_labels();h,l=extra.get_legend_handles_labels();axes[2].legend(handles+h,legend+l,fontsize=8)
    elif event.get('depth'):axes[2].legend(fontsize=8)
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
    displayed=[t for t in event.get('alternatives',[]) if t.get('site_verdict')!='EXCLUDE'][:6]
    for ax in axes:
        ax.set_facecolor('white')
        ax.axvspan(event['lo']/1e6,event['hi']/1e6,color='orange',alpha=.025 if event['hi']-event['lo']>200000 else .12,zorder=0)
        for cut in event['cut_options']:
            if cut.get('site_verdict')=='EXCLUDE':continue
            ax.axvline(cut['cut_bp']/1e6,color='#b2182b',ls='--')
        for location in displayed:
            ax.axvspan(location['lo']/1e6,location['hi']/1e6,facecolor='#762a83',alpha=.08,zorder=0)
            ax.axvline(location['probe_bp']/1e6,color='#762a83',ls='--',lw=1)
        for selected in selected_cuts:ax.axvline(int(selected['cut_bp'])/1e6,color='#007A3D',ls='-',lw=1.5)
        ax.grid(alpha=.15)
    groups=[];separation=.03*(event['hi']-event['lo']+1000000)
    for location in sorted(displayed,key=lambda t:t['probe_bp']):
        if groups and location['probe_bp']-groups[-1][-1]['probe_bp']<separation:groups[-1].append(location)
        else:groups.append([location])
    for group in groups:
        label='/'.join(t['id'].removeprefix(event['id']) for t in group)
        axes[0].annotate(label,(sum(t['probe_bp'] for t in group)/len(group)/1e6,1),xycoords=('data','axes fraction'),fontsize=7,ha='center',va='bottom',xytext=(0,3),textcoords='offset points',color='#762a83',annotation_clip=True)
    for ax in axes[1:]:ax.set_ylim(bottom=0)
    window=detail_window(event,selected_cuts) if detail else (max(0,event['lo']-500000),event['hi']+500000)
    axes[-1].set_xlim(window[0]/1e6,window[1]/1e6)
    axes[4].legend(handles=[Line2D([],[],color='#b2182b',ls='--',label='Exact gap option'),Line2D([],[],color='#762a83',ls='--',label='R assay position (not a cut)')]+([Line2D([],[],color='#007A3D',label='User-selected cut')] if selected_cuts else []),fontsize=7,loc='upper right')
    fig.suptitle(event.get('assembly','Assessment')+' · '+event['id']+' · '+event['scaffold']+' · '+event['left']+' → '+event['right'],y=.995)
    fig.subplots_adjust(left=.20,right=.92,top=.83 if other else .89,bottom=.06,hspace=.44)
    # Explicit final labels and tight export protect panel labels across renderers.
    axes[1].set_ylabel('HiFi spanning molecules');axes[2].set_ylabel('Mean HiFi depth')
    axes[3].set_ylabel('Hi-C contact support\nCross / sqrt(within L x R)')
    axes[3].set_title('Binned Hi-C support: 100 kb bins, 300 kb flanks',fontsize=10)
    axes[4].set_ylabel('Telomere motif count\nper 1 kb');axes[4].set_xlabel('Original scaffold coordinate (Mb)')
    name=event['id']+('.detail' if detail else '')+'.png'
    fig.savefig(out/'plots'/name,dpi=150,facecolor='white',bbox_inches='tight',pad_inches=.15);plt.close(fig)
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
    assembly=next((e.get('assembly','Assessment') for e in events if e['scaffold']==chrom),'Assessment')
    axis.set_title(assembly+' · '+chrom+' · contact map',pad=28)
    axis.set_ylabel('Scaffold position (Mb)');axis.set_xlabel('Scaffold position (Mb)')
    fig.tight_layout();fig.savefig(out/'plots'/(chrom+'.contacts.png'),dpi=150,facecolor='white');plt.close(fig)
    progress(f'Finished contact map {chrom} in {time.monotonic()-started:.1f}s')

