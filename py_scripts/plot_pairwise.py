"""Dotplot and ribbon overview in final-assembly coordinates; invalid evidence fails."""
import argparse
import csv
import gzip
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon, Rectangle
from pairwise_paf import parse_paf


def read_fai(path):
    result={}
    with open(path) as handle:
        for line in handle:
            f=line.rstrip().split('\t')
            if len(f)<2 or not f[1].isdigit() or int(f[1])<=0 or f[0] in result:
                raise ValueError(f'Invalid or duplicate FAI entry: {path}')
            result[f[0]]=int(f[1])
    if not result:raise ValueError(f'Empty FAI: {path}')
    return result


def read_alignments(path,ref,query):
    opener=gzip.open if str(path).endswith('.gz') else open
    rows=[]
    with opener(path,'rt') as handle:
        for line in handle:
            if not line.strip():continue
            r=parse_paf(line)
            if ref.get(r['tname'])!=r['tlen'] or query.get(r['qname'])!=r['qlen']:
                raise ValueError('PAF sequence names/lengths disagree with final assembly indexes')
            rows.append(r)
    return rows


def layout(entries,minimum):
    def key(item):
        match=re.match(r'^chr(\d+)',item[0])
        return (int(match[1]) if match else 10**9,-item[1],item[0])
    offsets={};lengths={};total=0
    for name,length in sorted(entries.items(),key=key):
        if length<minimum:continue
        offsets[name]=total/1e6;lengths[name]=length/1e6;total+=length
    return offsets,lengths,total/1e6


def overview(rows,ref,query,out,ref_id,query_id,kind,minimum,length_min,width,height,alpha=.45):
    ro,rl,rt=layout(ref,minimum);qo,ql,qt=layout(query,minimum)
    chosen=[r for r in rows if r['length']>=length_min and r['tname'] in ro and r['qname'] in qo]
    fig,ax=plt.subplots(figsize=(width,height),facecolor='white');ax.set_facecolor('white')
    status='plotted' if chosen else ('no_sequences_above_threshold' if not ro or not qo else 'no_alignments_passing_plot_filters')
    if not chosen:ax.text(.5,.5,status.replace('_',' '),ha='center',va='center',transform=ax.transAxes)
    if kind=='dotplot':
        for strand,color in [('+','#0072B2'),('-','#D55E00')]:
            segments=[]
            for r in chosen:
                if r['strand']!=strand:continue
                xs=[r['qs'],r['qe']] if strand=='+' else [r['qe'],r['qs']]
                segments.append([(qo[r['qname']]+xs[0]/1e6,ro[r['tname']]+r['ts']/1e6),
                                 (qo[r['qname']]+xs[1]/1e6,ro[r['tname']]+r['te']/1e6)])
            if segments:ax.add_collection(LineCollection(segments,colors=color,linewidths=.7,label=f'{strand} strand'))
        for x in qo.values():ax.axvline(x,color='#dddddd',lw=.4,zorder=0)
        for y in ro.values():ax.axhline(y,color='#dddddd',lw=.4,zorder=0)
        ax.set_xlim(0,max(qt,.001));ax.set_ylim(0,max(rt,.001))
        ax.set_aspect('equal',adjustable='box')
        ax.set_xlabel(f'{query_id} — concatenated final sequence (Mb)')
        ax.set_ylabel(f'{ref_id} — concatenated final sequence (Mb)')
        if chosen:ax.legend(loc='best')
    else:
        colors={name:plt.get_cmap('tab20')((2*i)%20+(i//10)%2) for i,name in enumerate(ro)}
        for r in chosen:
            top=[ro[r['tname']]+r['ts']/1e6,ro[r['tname']]+r['te']/1e6]
            bottom=[qo[r['qname']]+r['qs']/1e6,qo[r['qname']]+r['qe']/1e6]
            if r['strand']=='-':bottom.reverse()
            ax.add_patch(Polygon([(top[0],.95),(top[1],.95),(bottom[1],.05),(bottom[0],.05)],
                                 facecolor=colors[r['tname']],edgecolor='none',alpha=alpha))
        for y,offsets,lengths in [(1,ro,rl),(0,qo,ql)]:
            for name,x in offsets.items():
                ax.add_patch(Rectangle((x,y-.025),lengths[name],.05,facecolor=colors.get(name,'#777777'),
                                       edgecolor='white',lw=.4,zorder=3))
                if len(offsets)<=50:ax.text(x+lengths[name]/2,y+.065,name,ha='center',rotation=35,fontsize=6)
        ax.set_xlim(0,max(rt,qt,.001));ax.set_ylim(-.2,1.4)
        ax.set_yticks([0,1],[query_id,ref_id]);ax.set_xlabel('Concatenated final sequence (Mb)')
    ax.set_title(f'{query_id} vs {ref_id}\n'+f'{len(chosen):,} alignments; sequences ≥ {minimum:,} bp; alignments ≥ {length_min:,} bp',fontsize=10)
    fig.tight_layout();fig.savefig(out,dpi=180,facecolor='white');plt.close(fig)
    return dict(plot=kind,status=status,alignments=len(chosen),ref_sequences=len(ro),query_sequences=len(qo),
                ref_axis_bp=round(rt*1e6),query_axis_bp=round(qt*1e6))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['paf','ref-fai','query-fai','ref','query','prefix']:p.add_argument('--'+name,required=True)
    p.add_argument('--dot-min-seq',type=int,default=1000000)
    p.add_argument('--ribbon-min-seq',type=int,default=1000000)
    p.add_argument('--ribbon-min-aln',type=int,default=150000)
    p.add_argument('--ribbon-alpha',type=float,default=.45)
    p.add_argument('--dot-width',type=float,default=10);p.add_argument('--dot-height',type=float,default=10)
    p.add_argument('--ribbon-width',type=float,default=14);p.add_argument('--ribbon-height',type=float,default=6)
    a=p.parse_args()
    if min(a.dot_min_seq,a.ribbon_min_seq,a.ribbon_min_aln)<0:p.error('Negative plotting thresholds')
    if min(a.dot_width,a.dot_height,a.ribbon_width,a.ribbon_height)<=0:p.error('Plot dimensions must be positive')
    if not 0<=a.ribbon_alpha<=1:p.error('Ribbon alpha must be between zero and one')
    try:
        ref=read_fai(a.ref_fai);query=read_fai(a.query_fai);rows=read_alignments(a.paf,ref,query)
        stats=[overview(rows,ref,query,a.prefix+'_dotplot.png',a.ref,a.query,'dotplot',
                        a.dot_min_seq,0,a.dot_width,a.dot_height),
               overview(rows,ref,query,a.prefix+'_riparian.png',a.ref,a.query,'ribbon',
                        a.ribbon_min_seq,a.ribbon_min_aln,a.ribbon_width,a.ribbon_height,a.ribbon_alpha)]
        with open(a.prefix+'.plot_status.tsv','w',newline='') as out:
            w=csv.DictWriter(out,stats[0].keys(),delimiter='\t');w.writeheader();w.writerows(stats)
    except (ValueError,OSError) as exc:p.exit(65,str(exc)+'\n')


if __name__=='__main__':main()
