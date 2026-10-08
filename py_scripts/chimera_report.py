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


def family_row(group,cohort=False):
    first=group[0]
    proposed=[r for r in group if r['preferred_candidate']==r['id']]
    available=[r for r in group if r['cut_bp']!='']
    cut=(str(proposed[0]['cut_bp'])+' ('+proposed[0]['id']+', preferred)') if len(proposed)==1 else ('; '.join(r['id']+': '+str(r['cut_bp']) for r in available)+'; association unresolved' if available else 'Not localized')
    links=' / '.join('['+r['id']+']('+ (r['report_path'] if cohort else '#candidate-'+r['id'].lower())+')' for r in group)
    value=[links,first['scaffold'],first['detected_transition'] or first['chromosome_context'],
           '; '.join(r['id']+': '+r['review_range'] for r in group),cut,
           '; '.join(r['id']+': '+r['local_confirmation']+'; '+r['partial_confirmation'] for r in group),
           '; '.join(r['id']+': '+r['structural_context'] for r in group)]
    return [first['assembly']]+value if cohort else value


def track_svg(row,measurement,interval,out):
    tracks=measurement.get('chromosome_tracks',[])
    if not tracks:return
    start=interval.get('start',int(row['review_start'] or 0))
    length=max((b['hi'] for t in tracks for b in t.get('bins',[])),default=1)
    chromosomes=sorted({b['chrom'] for t in tracks for b in t.get('bins',[]) if b.get('chrom')})
    colors=['#4477aa','#ee6677','#228833','#ccbb44','#66ccee','#aa3377','#332288','#44aa99','#882255']
    palette={c:colors[i%len(colors)] for i,c in enumerate(chromosomes)}
    height=100+28*len(tracks)+25*((len(chromosomes)+5)//6)
    svg=['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 %d">'%height,
         '<style>text{font:13px sans-serif}</style>']
    x=lambda pos:220+980*pos/length
    for i in range(5):
        pos=length*i/4;xx=x(pos)
        svg+=['<text x="%.1f" y="20" text-anchor="middle">%.3f Mb</text>'%(xx,(start+pos)/1e6)]
    for i,t in enumerate(tracks):
        y=35+28*i
        svg+=['<text x="0" y="%d">%s</text>'%(y+14,html.escape(t['peer'])),
              '<rect x="220" y="%d" width="980" height="18" fill="#f5f5f5"/>'%y]
        for b in t.get('bins',[]):
            if not b.get('aligned_bp') and not b.get('ambiguous_bp'):continue
            svg+=['<rect x="%.1f" y="%d" width="%.1f" height="18" fill="%s"><title>%s</title></rect>'%
                  (x(b['lo']),y,980*(b['hi']-b['lo'])/length,palette.get(b.get('chrom'),'#aaa'),html.escape(str(b)))]
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
            ax.scatter([ratio(focal)],[.5],s=90,marker='D',color='red',label='Candidate: within L/R '+focal['left_within']+'/'+focal['right_within'])
        elif focal:ax.text(.02,.8,'Candidate ratio undefined: zero within-flank count',transform=ax.transAxes)
        maximum=max([maximum]+[ratio(r) for r in controls]+([ratio(focal)] if focal and ratio(focal) is not None else []))
        ax.set_xscale('symlog',linthresh=.0001);ax.set_yticks([]);ax.set_ylim(-.15,.8);ax.set_title(lib);ax.legend(fontsize=8,loc='upper right')
        ax.set_xlabel('Cross pairs / sqrt(within-left × within-right); zero is measured zero')
        ax.set_xticks([0,.0001,.001,.01,.1],['0','0.0001','0.001','0.01','0.1'])
    axes[0,0].set_xlim(-.00001,maximum*1.4)
    fig.suptitle(row['assembly']+' / '+row['id']+' — available controls; assay may be uncalibrated')
    fig.tight_layout();fig.savefig(out/(row['id']+'.control_distribution.png'),dpi=130);plt.close(fig)


def render(assembly,rows,measurements,provenance,sections,context,out):
    checksum=provenance.get('sha256') or next((r['assessment_sha256'] for r in rows),'')
    audit=provenance.get('coordinate_audit',{})
    state=audit.get('alignment_transition_assessment','not recorded')
    if provenance.get('chromosome_scope',{}).get('unresolved'):state='chromosome scope unresolved'
    reason=('No chromosome-scale composite passed the detector screen; local boundary assays were not needed' if not rows and state=='assessed' else
            'Reference alignment unavailable; absence of candidates is not a negative assessment' if state=='unavailable' else
            'Chromosome-scale inference unresolved: '+', '.join(provenance['chromosome_scope']['reasons']) if state=='chromosome scope unresolved' else
            'Assessment scope not recorded in this packet' if not rows else 'Candidate-scoped local evidence generated')
    registry=dict(assembly=assembly,assessment_sha256=checksum,coordinate_stage='pre_finishing',candidate_count=len(rows),assessment_status=state,assessment_reason=reason)
    (out/'registry.json').write_text(json.dumps(registry,indent=2),encoding='utf-8')
    groups=families(rows);retained=[r for r in rows if r['assessment_status']=='retention_supported']
    md=['# '+assembly+' — chimera evidence','','[Cohort report](../README.md) · [Manual cut instructions](../cut-instructions.md)','','## Assessment summary','',
        str(len(rows))+' candidate boundaries assessed; '+str(len(groups))+' transition families require review; '+str(len(retained))+' boundaries have evidence against a chromosome fusion.',
        'All selections start at NO. Human review disposition starts at PENDING. Evidence generation applies no cuts.','',
        '**Assessment scope:** '+state+'. '+reason+'.',
        'Coordinates are zero-based on the original pre-finishing scaffold. Chromosome identity, local confirmation and cut localization are separate assessments.','']
    if groups:
        md+=['## Transition families requiring review','']+table(['Measurements','Scaffold','Detected chromosomes','Measured intervals, bp','Available gap cut, bp','Independent chromosome observations','Native context'],[family_row(g) for g in groups])
    if retained:
        md+=['','## Boundaries with evidence against a chromosome fusion','']+table(['Measurement','Scaffold','Initial detection','Local chromosome context','Assessed gap / interval','Interpretation'],[
            ['['+r['id']+'](#candidate-'+r['id'].lower()+')',r['scaffold'],r['detected_transition'],r['chromosome_context'],r['review_range'],'Opposes this chromosome-fusion hypothesis; exact adjacency not validated'] for r in retained])
    for group in groups:
        bridge=[measurements[r['source_candidate_id']].get('transition_bridge_assessment') for r in group if r['source_candidate_id'] in measurements]
        for b in filter(None,bridge):
            md+=['','### Related measurements '+b['source']+' / '+b['gap'],'',
                 'Bridge '+format(b['start'],',')+'–'+format(b['end'],',')+' bp: '+b['status'].replace('_',' ')+'. Grouping records shared detector provenance; it does not establish a preferred cut.','']
            md+=table(['Peer','Assigned bridge bp','Bridge coverage','Discordant bp'],[
                [peer,bp,'%.1f%%'%(100*bp/max(1,b['end']-b['start'])),b.get('discordant_bp',{}).get(peer,0)] for peer,bp in b.get('coverage_bp',{}).items()])
            md+=['','Provisional precedence screen: at least 50% assigned bridge coverage in an independent peer without discordant segments; ≥10 kb discordant assignment or an intervening transition blocks precedence. These grouping thresholds are not an automated cut classifier.']
    intervals=provenance.get('intervals',{})
    for r in rows:
        m=measurements.get(r['source_candidate_id'],{});interval=intervals.get(m.get('packet_interval_id'),{})
        track_svg(r,m,interval,out);local_plot(r,m,interval,context,out);control_plot(r,m,context,out)
        md+=['','## Candidate '+r['id'],'', '**Location:** '+r['scaffold']+':'+r['review_range']+' bp. **Family:** '+r['transition_id']+'.',
             '**Detected chromosomes:** '+(r['detected_transition'] or 'Nearby gap hypothesis')+'. **Measured local context:** '+r['chromosome_context']+'.',
             '**Native context:** '+r['structural_context']+'. Native contig continuity describes the assembly hypothesis, not independent molecule support.',
             '**Local confirmation:** '+r['local_confirmation']+'.', '**Additional observations:** '+r['partial_confirmation']+'.',
             '**Sister haplotype:** '+r['sister_context']+'. Same-individual evidence is not an independent vote.',
             '**Evidence for reviewing a fusion:** '+r['evidence_for_cut']+'.', '**Evidence against a fusion / for continuity:** '+r['evidence_against_cut']+'.',
             '**Review priority:** '+r['review_priority']+'.', '**Localization:** '+r['localization_explanation']+'.',
             '**Available gap:** '+(str(r['gap_start'])+'–'+str(r['gap_end']) if r['gap_start']!='' else 'None established')+'. **Proposed cut:** '+(str(r['cut_bp']) if r['cut_bp']!='' else 'None')+'. **Selected:** NO. **Human review:** PENDING.', '',
             '### Peer chromosome assignments','']
        md+=table(['Peer','Role / eligibility','Left → right','Aligned kb L/R','Dominance L/R','Informative bins L/R','Qualification','Homologous placement'],[
            [t['peer'],'Same individual' if t['sample']==provenance.get('sample') else t.get('comparison_scope','Independent eligible' if t.get('auto_evidence') else 'Context only'),
             str(t['left'].get('chrom') or '?')+' → '+str(t['right'].get('chrom') or '?'),
             '%.1f / %.1f'%(t['left'].get('aligned_bp',0)/1000,t['right'].get('aligned_bp',0)/1000),
             '%.2f / %.2f'%(t['left'].get('dominance',0),t['right'].get('dominance',0)),
             str(t['left'].get('informative_bins',0))+' / '+str(t['right'].get('informative_bins',0)),
             ('Both sides qualify' if t['relationship']!='uninformative' else 'L '+str(t['left'].get('qualified',False))+' / R '+str(t['right'].get('qualified',False))),
             str(t.get('placement_relationship','Not measured'))+'; '+str(t.get('left_target') or '?')+' / '+str(t.get('right_target') or '?')]
            for t in m.get('chromosome_tracks',[])])
        md+=['','Local chromosome qualification requires ≥100 kb assigned bases, ≥90% dominance and ≥3 informative 50 kb bins on each side. Below-threshold identities remain observations. Unlabelled target blocks and filtered alignments do not establish biological absence.','']
        raw=m.get('hifi_raw',{})
        evidence_rows=[
            ['HiFi immediate flank molecules L/R',str(raw.get('left_molecules','Not measured'))+' / '+str(raw.get('right_molecules','Not measured'))],
            ['HiFi molecules spanning assessed interval',m.get('hifi_spanning_molecules','Not measured')],
            ['Flank observability',m.get('hifi_informative','Not measured')],
            ['Local continuity grid minimum molecules',m.get('continuity_grid',{}).get('minimum_molecules','Not measured')],
            ['Native primary contig continuity',m.get('native_continuity','Not measured')],
            ['Native contigs L/R',str(m.get('native_left','?'))+' / '+str(m.get('native_right','?'))],
            ['Direct native graph link',m.get('direct_native_link','Not measured')]]
        if not raw and 'hifi_spanning_molecules' not in m and not m.get('continuity_grid'):
            evidence_rows=evidence_rows[4:]
        md+=table(['Evidence','Measured result'],evidence_rows)
        if raw or 'hifi_spanning_molecules' in m:
            md+=['','Flank observability is not a test that a read could span the whole interval. Interpret zero spanning reads with interval width, molecule lengths, repeats and local continuity.','']
        contact_trials=m.get('farther_contact_evidence',{}).get('trials',[])
        if contact_trials:
            md+=table(['Library','Offset kb','Cross pairs','Within pairs L/R','Qualified controls sequence/gap','Calibration'],[
            [t['library'],t['offset_bp']//1000,t['raw_counts'].get('cross',0),str(t['raw_counts'].get('left_within',0))+' / '+str(t['raw_counts'].get('right_within',0)),
             str(t['control_populations']['continuous_control'])+' / '+str(t['control_populations']['gap_control']),
             'Informative' if t['informative'] else 'Insufficient observability or matched controls'] for t in contact_trials])
            md+=['','Zero counts are measured zeros after a completed scan. Unmeasured assays have no trial rows. Scaffolding Hi-C is corroboration; uncalibrated support loss does not establish a cut.','']
        for suffix,title in [('.tracks.svg','Chromosome assignments in original scaffold coordinates'),('.local_support.png','Local HiFi depth and continuity'),('.control_distribution.png','Candidate versus available control distributions')]:
            if (out/(r['id']+suffix)).exists():md+=['!['+title+']('+r['id']+suffix+')','']
        packet=m.get('packet_interval_id');link='../../sequence_context/'+assembly+'.sequence_context/'
        if packet and (context/(packet+'.igv.xml')).exists():md+=['[IGV session in original coordinates]('+link+packet+'.igv.xml). Paths may need updating when moving the data.','']
        md+=['**Limits:** '+r['evidence_limits']+'.',
             ('**Review disposition:** evidence favors retention for the chromosome-fusion question; optional adjacency review remains separate.' if r['assessment_status']=='retention_supported' else '**Decision:** review this family before selecting a cut. A reported gap location does not authorize cutting.'),
             '[Detailed machine-readable evidence](measurements.json) (key `'+r['source_candidate_id']+'`).','']
    supplement=out/'scaffold_evidence'
    images=sorted(supplement.glob('*.png')) if supplement.exists() else []
    if images:
        md+=['## Scaffold-wide contact maps, cross-contact profiles and telomeres','',
             'These panels use the pre-finishing scaffold and pool scaffolding libraries. Per-library measurements above remain separate. Markers indicate diagnostic locations, not selected cuts. Telomere motif counts are descriptive and do not establish a fusion.','']
        for image in images:
            md+=['!['+image.stem+'](scaffold_evidence/'+image.name+')','']
        md+=['### Scaffold evidence records','']
        for path in sorted(supplement.glob('*.tsv')):md+=['- ['+path.name+'](scaffold_evidence/'+path.name+')']
    elif rows and provenance.get('supplementary_status')!='disabled':
        md+=['## Scaffold-wide evidence availability','','No supplementary scaffold figures were supplied. See the pipeline evidence outputs/logs; this is not negative telomere or Hi-C evidence.']
    md+=['','## Source identity and review interface','', '**Assessment SHA-256:** `'+checksum+'`.',
         '[Editable rows](review.tsv) · [Full measurements](measurements.json) · [Assessment provenance and scope](assessment_provenance.json).',
         'Review disposition is PENDING, RETAIN, CUT or DEFER. Cut selection is a separate YES/NO field. Marking a review disposition alone never modifies sequence.',
         'Native primary-path graph assays do not measure unitig read support. Composite peer chromosome blocks come from reference alignment anchors, not from assigning one chromosome label to an entire composite scaffold.']
    (out/'report.md').write_text('\n'.join(md)+'\n',encoding='utf-8')
