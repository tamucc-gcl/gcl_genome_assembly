"""Compact, static evidence reports; detailed measurements stay in linked data files."""
import csv
import html
import json
from collections import defaultdict


def table(headers,rows):
    def cell(value):return str(value).replace('|','&#124;').replace('\n',' ')
    return ['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(cell,row))+' |' for row in rows]


def families(rows):
    groups=defaultdict(list)
    for row in rows:
        if row['assessment_status']!='retention_supported':groups[row['transition_id']].append(row)
    return list(groups.values())


def track_svg(row,measurement,interval,out):
    tracks=measurement.get('chromosome_tracks',[])
    if not tracks:return
    start=interval.get('start',int(row['review_start'] or 0))
    length=max((b['hi'] for t in tracks for b in t.get('bins',[])),default=1)
    chromosomes=sorted({b['chrom'] for t in tracks for b in t.get('bins',[]) if b.get('chrom')})
    colors=['#4477aa','#ee6677','#228833','#ccbb44','#66ccee','#aa3377','#332288','#44aa99','#882255']
    palette={c:colors[i%len(colors)] for i,c in enumerate(chromosomes)}
    height=100+28*len(tracks)+25*((len(chromosomes)+7)//6)
    svg=['<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="%d" viewBox="0 0 1200 %d">'%(height,height),
         '<rect width="1200" height="%d" fill="white"/>'%height,
         '<style>text{font:14px sans-serif;fill:#111}</style>']
    x=lambda pos:220+940*pos/length
    for i in range(5):
        pos=length*i/4;xx=x(pos)
        svg+=['<text x="%.1f" y="20" text-anchor="middle">%.3f Mb</text>'%(xx,(start+pos)/1e6)]
    for i,t in enumerate(tracks):
        y=35+28*i
        svg+=['<text x="0" y="%d">%s</text>'%(y+14,html.escape(t['peer'])),
              '<rect x="220" y="%d" width="940" height="18" fill="#f5f5f5"/>'%y]
        for b in t.get('bins',[]):
            if not b.get('aligned_bp') and not b.get('ambiguous_bp'):continue
            svg+=['<rect x="%.1f" y="%d" width="%.1f" height="18" fill="%s"><title>%s</title></rect>'%
                  (x(b['lo']),y,940*(b['hi']-b['lo'])/length,palette.get(b.get('chrom'),'#aaa'),html.escape(str(b)))]
    for boundary in (row['review_start'],row['review_end']):
        if boundary!='':svg+=['<line x1="%.1f" x2="%.1f" y1="30" y2="%d" stroke="black" stroke-dasharray="4 3"/>'%(x(int(boundary)-start),x(int(boundary)-start),35+len(tracks)*28)]
    y=65+28*len(tracks)
    for i,(chrom,color) in enumerate(list(palette.items())+[('Ambiguous / below threshold','#aaa'),('No assigned alignment','#f5f5f5')]):
        xx=(i%6)*200;yy=y+(i//6)*25
        svg+=['<rect x="%d" y="%d" width="12" height="12" fill="%s"/><text x="%d" y="%d">%s</text>'%(xx,yy,color,xx+17,yy+12,chrom)]
    svg+=['</svg>'];(out/(row['id']+'.tracks.svg')).write_text(''.join(svg),encoding='utf-8')


def local_plot(row,m,interval,context,out):
    depth=context/(m.get('packet_interval_id','')+'.depth.tsv')
    grid=m.get('continuity_grid',{})
    if not depth.exists() and not grid:return
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,1,figsize=(10,4.5),sharex=True)
    if depth.exists():
        points=[]
        with depth.open(encoding='utf-8') as handle:
            for i,line in enumerate(handle):
                if i%100==0:
                    f=line.split('\t');points.append((int(f[1])-1,int(f[2])))
        axes[0].plot([p/1e6 for p,d in points],[d for p,d in points],lw=.6)
    else:axes[0].text(.5,.5,'HiFi depth not measured',transform=axes[0].transAxes,ha='center')
    axes[0].set_ylabel('HiFi depth (MAPQ20)')
    sites=grid.get('probes',[])
    if sites and isinstance(sites[0],dict):
        axes[1].plot([s.get('position',s.get('cut_bp',0))/1e6 for s in sites],[s.get('molecules',s.get('spanning',0)) for s in sites],lw=.8)
    else:axes[1].text(.5,.5,'Local continuity: minimum molecules '+str(grid.get('minimum_molecules','not measured')),transform=axes[1].transAxes,ha='center')
    axes[1].set_ylabel('Spanning molecules')
    for ax in axes:
        if row['review_start']!='':ax.axvspan(int(row['review_start'])/1e6,int(row['review_end'])/1e6,color='orange',alpha=.25,label='Assessed interval')
        ax.grid(alpha=.15)
    axes[-1].set_xlabel('Original scaffold position (Mb)');axes[0].legend(loc='upper right')
    fig.suptitle(row['assembly']+' / '+row['scaffold']+' / '+row['id']);fig.tight_layout()
    fig.savefig(out/(row['id']+'.local_support.png'),dpi=130);plt.close(fig)


def control_plot(row,m,context,out):
    path=context/'contact_controls.tsv'
    if not path.exists() or not m.get('control_ids'):return
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    with path.open(encoding='utf-8') as handle:records=list(csv.DictReader(handle,delimiter='\t'))
    ids=set(m['control_ids']);key=m.get('packet_interval_id');libraries=sorted({r['library'] for r in records})
    role_path=context/'control_registry.json'
    roles=json.loads(role_path.read_text(encoding='utf-8')).get('intervals',{}) if role_path.exists() else {}
    fig,axes=plt.subplots(len(libraries),1,figsize=(10,3.3*len(libraries)),squeeze=False,sharex=True)
    maximum=.0001
    for ax,lib in zip(axes[:,0],libraries):
        controls=[r for r in records if r['library']==lib and r['candidate'] in ids]
        focal=next((r for r in records if r['library']==lib and r['candidate']==key),None)
        def ratio(r):
            denominator=(int(r['left_within'])*int(r['right_within']))**.5
            return int(r['cross'])/denominator if denominator else None
        controls=[r for r in controls if ratio(r) is not None]
        for role,y,color in [('gap_control',0,'tab:blue'),('continuous_control',.18,'tab:green')]:
            values=[r for r in controls if roles.get(r['candidate'],{}).get('role')==role]
            if values:ax.scatter([ratio(r) for r in values],[y]*len(values),alpha=.5,color=color,label=role.replace('_',' ')+' (n='+str(len(values))+')')
        if focal and ratio(focal) is not None:
            ax.scatter([ratio(focal)],[.5],s=90,marker='D',color='red',label='Tested join: '+focal['cross']+' cross-join pairs')
        elif focal:ax.text(.02,.8,'Candidate ratio undefined: zero within-flank count',transform=ax.transAxes)
        maximum=max([maximum]+[ratio(r) for r in controls]+([ratio(focal)] if focal and ratio(focal) is not None else []))
        ax.set_xscale('symlog',linthresh=.0001);ax.set_yticks([0,.18,.5],['Other gap joins','Continuous sequence','Tested join']);ax.set_ylim(-.15,.8);ax.set_title(lib);ax.legend(fontsize=8,loc='upper right')
        ax.set_xlabel('Normalized cross-join contacts (lower = fewer contacts between sides)')
        ax.set_xticks([0,.0001,.001,.01,.1],['0','0.0001','0.001','0.01','0.1'])
    axes[0,0].set_xlim(-.00001,maximum*1.4)
    fig.suptitle(row['assembly']+' / '+row['id']+' — exploratory Hi-C comparison\nAvailable controls are not necessarily matched or qualified')
    fig.tight_layout();fig.savefig(out/(row['id']+'.control_distribution.png'),dpi=130,facecolor='white');plt.close(fig)



def render(assembly,rows,measurements,provenance,sections,context,out):
    from chimera_decisions import assessment_label, advise, bp, decision_card, focus, location
    checksum=provenance.get('sha256') or next((r['assessment_sha256'] for r in rows),'')
    state=provenance.get('coordinate_audit',{}).get('alignment_transition_assessment','not recorded')
    if provenance.get('chromosome_scope',{}).get('unresolved'):state='chromosome scope unresolved'
    groups=families(rows);rejected=[r for r in rows if r['assessment_status']=='retention_supported']
    registry=dict(assembly=assembly,assessment_sha256=checksum,coordinate_stage='pre_finishing',candidate_count=len(rows),assessment_status=state,assessment_reason=state,
                  potential_chimeras=len(groups),rejected_signals=len(rejected))
    registry['assessment_reason']=assessment_label(registry,groups)
    (out/'registry.json').write_text(json.dumps(registry,indent=2),encoding='utf-8')
    sample=provenance.get('sample',assembly.rsplit('_hap',1)[0])
    md=['# '+assembly+' — chimera review','','[All assemblies](../README.md) · [Cut instructions](../cut-instructions.md)','',
        '**'+registry['assessment_reason']+'.**', '',
        'All selections start at NO. No cuts have been applied. Suggested actions below are advice for human review.','']
    if state=='unavailable':md+=['This assembly could not be screened because its reference alignment is unavailable. Zero detections here does not mean it is free of misjoins.','']
    elif state=='chromosome scope unresolved':md+=['The pipeline could not identify chromosome-scale scaffolds for this assembly. Chromosome-fusion screening was not completed.','']
    elif not groups:md+=['No chromosome-fusion decision is needed from this screen. This does not validate every assembly join.','']
    if groups:
        md+=['## Suggested actions','']+table(['Event','Scaffold / chromosomes','Suggested action','Cut or transition range'],[
            [' / '.join('['+r['id']+'](#candidate-'+r['id'].lower()+')' for r in g),g[0]['scaffold']+' · '+g[0]['detected_transition'],
             advise(focus(g,measurements,sample),measurements.get(focus(g,measurements,sample)['source_candidate_id'],{}),sample)[0]['recommendation'],
             location(focus(g,measurements,sample))] for g in groups])
        md+=['','Coordinates below use the original scaffold before gap filling or finishing. A range is where the chromosome assignment changes; it is not permission to cut at its midpoint.','']
    intervals=provenance.get('intervals',{})
    for group in groups:
        for r in group:md+=['<a name="candidate-'+r['id'].lower()+'"></a>']
        md+=['']+decision_card(group,measurements,provenance)
        primary=focus(group,measurements,sample)
        for r in group:
            m=measurements.get(r['source_candidate_id'],{});interval=intervals.get(m.get('packet_interval_id'),{})
            track_svg(r,m,interval,out);local_plot(r,m,interval,context,out);control_plot(r,m,context,out)
            if r is primary:
                for suffix,title in [('.tracks.svg','Which chromosomes align on each side'),('.local_support.png','HiFi read depth and local continuity'),('.control_distribution.png','Hi-C contacts compared with controls')]:
                    if (out/(r['id']+suffix)).exists():md+=['!['+title+']('+r['id']+suffix+')','']
        md+=['<details>','<summary>Supporting measurements, other assemblies and assay details</summary>','']
        for r in group:
            m=measurements.get(r['source_candidate_id'],{})
            md+=['#### '+r['id']+' — '+location(r),'']
            tracks=m.get('chromosome_tracks',[])
            if tracks:
                md+=table(['Other assembly','Chromosomes left → right','Local alignment test','Role'],[
                    [t['peer'],str(t['left'].get('chrom') or '?')+' → '+str(t['right'].get('chrom') or '?'),
                     'Enough alignment on both sides' if t['relationship']!='uninformative' else 'Too little or ambiguous alignment',
                     'Sister haplotype' if t['sample']==sample else 'Independent individual' if t.get('auto_evidence') else 'Context only'] for t in tracks])
                md+=['','The local test needs 100 kb of assigned alignment, 90% agreement and three informative 50 kb bins per side. Haplotypes of one individual are counted together. Sister-haplotype assignments that fail this test do not establish the local chromosome relationship.','']
            native=m.get('native_continuity')
            if native is not None:md+=['**Assembly structure:** '+r['structural_context']+' ('+str(m.get('native_left','?'))+' / '+str(m.get('native_right','?'))+').','']
            trials=m.get('farther_contact_evidence',{}).get('trials',[])
            if trials:
                md+=table(['Library','Distance from boundary','Cross pairs','Within pairs left / right','Matched controls continuous / gap','Cut test'],[
                    [t['library'],str(t['offset_bp']//1000)+' kb',t['raw_counts'].get('cross',0),str(t['raw_counts'].get('left_within',0))+' / '+str(t['raw_counts'].get('right_within',0)),
                     str(t['control_populations']['continuous_control'])+' / '+str(t['control_populations']['gap_control']),
                     'Supports cutting' if t.get('informative') and t.get('support_loss') else 'No contact loss detected' if t.get('informative') else 'Cannot classify'] for t in trials])
                md+=['','A measured zero is shown as 0. “Cannot classify” means coverage or matched controls are inadequate, rather than evidence for retaining the join.','']
            if r is not primary:
                for suffix in ('.tracks.svg','.local_support.png','.control_distribution.png'):
                    if (out/(r['id']+suffix)).exists():md+=['!['+r['id']+' supporting plot]('+r['id']+suffix+')','']
            packet=m.get('packet_interval_id')
            if packet and (context/(packet+'.igv.xml')).exists():md+=['[Open '+r['id']+' in IGV](../../sequence_context/'+assembly+'.sequence_context/'+packet+'.igv.xml). Update paths if the files have moved.','']
        md+=['</details>','']
    if rejected:
        md+=['<details>','<summary>Screened out: '+str(len(rejected))+' reference-alignment signals — no cuts suggested</summary>','',
             'These are false positives for the chromosome-fusion screen. The initial reference alignment changed chromosome label, but comparisons with other individuals place both sides on the same chromosome. Keep these joins. This does not establish whether the small segment is a repeat, a genuine insertion or an alignment error.','']
        for r in rejected:
            m=measurements.get(r['source_candidate_id'],{})
            md+=['<a name="candidate-'+r['id'].lower()+'"></a>','',
                 '- **'+r['id']+': keep.** '+r['scaffold']+', '+location(r)+'. Initial signal: '+r['detected_transition']+'. Local result: '+r['chromosome_context']+'.']
            detector=m.get('detector_context',{})
            if detector.get('arm_pattern')=='out_and_back':
                minor=min(int(detector.get('left_arm_aligned_bp',0)),int(detector.get('right_arm_aligned_bp',0)))
                md+=['  A '+bp(minor)+' bp minor alignment segment sits within the surrounding chromosome; it is not evidence for two chromosome-sized arms.']
        md+=['','These records remain in review.tsv for audit and deliberate manual overrides. They are excluded from the potential-chimera count and suggested-action list.','</details>','']
    images=sorted((out/'scaffold_evidence').glob('*.png'))
    if images:
        md+=['## Whole-scaffold evidence','','<details>','<summary>Contact maps and telomere profiles</summary>','',
             'These plots use the original scaffold and pool Hi-C libraries. Markers show tested regions, not selected cuts. Telomere motif peaks are descriptive evidence.','']
        for image in images:md+=['!['+image.stem+'](scaffold_evidence/'+image.name+')','']
        md+=['</details>','']
    elif rows and provenance.get('supplementary_status')!='disabled':md+=['**Whole-scaffold plots:** unavailable in this packet.','']
    md+=['## Record your decision','',
         'Use [review.tsv](review.tsv). To cut, enter review_disposition=CUT, selected=YES, reviewer and reason on the chosen row. For a boundary without an exact cut, supply an exact coordinate and the required cutting action first. To keep or defer, record RETAIN or DEFER and leave selected=NO.',
         '', '[Instructions and an example](../cut-instructions.md) explain how to apply decisions and add a cut that was not detected.','',
         '<details>','<summary>Assessment scope and source files</summary>','',
         '**Assessment SHA-256:** `'+checksum+'`.', '',
         '[Full measurements](measurements.json) · [Assessment provenance](assessment_provenance.json).',
         '', 'Screening targets chromosome-scale composite scaffolds. It does not assess every small scaffold or every join. Chromosome labels inside composite peer scaffolds come from reference alignments; unavailable labels are not evidence of biological absence.',
         '', '</details>']
    (out/'report.md').write_text('\n'.join(md)+'\n',encoding='utf-8')
