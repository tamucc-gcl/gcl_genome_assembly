"""Display-only decomposition of retained GX taxonomy; no new search or cleaning rules."""
import csv
import json
import re
from pathlib import Path

CATEGORIES=('primary','repeat_expected','low_expected','repeat_other','inconclusive','taxonomy_flag','unreported')

def table(path):
    with Path(path).open() as f:return list(csv.DictReader(f,delimiter='\t'))

def write(path,rows,fields):
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)

def location(identifier,lengths):
    name=identifier;parts=[]
    while name not in lengths:
        m=re.fullmatch(r'(.+?)(?:~~|~)(\d+)\.\.(\d+)',name)
        if not m:raise ValueError('GX segment absent from exact metrics dictionary: '+identifier)
        parts.append((int(m[2]),int(m[3])));name=m[1]
    lo,hi=0,lengths[name]
    for a,b in reversed(parts):
        if not 1<=a<=b<=hi-lo:raise ValueError('GX segment outside parent frame: '+identifier)
        lo,hi=lo+a-1,lo+b
    return name,lo,hi

def category(result,division,expected):
    if result.startswith('primary-div'):return 'primary'
    if result.startswith('contaminant') or result=='same-kingdom-chimeric':return 'taxonomy_flag'
    if expected and division==expected:
        if result=='repeat':return 'repeat_expected'
        if result=='low-coverage':return 'low_expected'
    if expected and division not in ('','n/a','NA',expected) and result in ('repeat','low-coverage'):
        return 'repeat_other'
    return 'inconclusive'

def enrich(sequences,fcs,detail,stage):
    """Validate a non-overlapping span partition, retaining raw results and diagnostics."""
    segments=[]
    for assembly in sorted({r['cohort_assembly'] for r in sequences}):
        current=[r for r in sequences if r['cohort_assembly']==assembly]
        lengths={r['sequence']:int(r['length_bp']) for r in current}
        for r in current:
            r['is_mixed']=int(r.get('plot_taxon')=='mixed FCS classifications')
            for c in CATEGORIES:r['cat_'+c+'_bp']=0
        lines=Path(fcs[assembly]['taxonomy']).read_text().splitlines()
        expected=''
        if lines and lines[0].startswith('##[['):
            expected=json.loads(lines[0][2:])[1].get('run-info',{}).get('asserted-div','')
        header=next(x for x in lines if x.startswith('#seq-id\t'))
        spans={name:[] for name in lengths}
        for row in csv.DictReader([header[1:]]+[x for x in lines if x and not x.startswith('#')],delimiter='\t'):
            name,lo,hi=location(row['seq-id'],lengths)
            if hi-lo!=int(row['seq-len']):raise ValueError('GX segment length differs from coordinates: '+row['seq-id'])
            rec=dict(cohort_assembly=assembly,coordinate_stage=stage,sequence=name,start0=lo,end0=hi,
                     span_bp=hi-lo,result=row['result'],division=row['div'],expected_division=expected,
                     category=category(row['result'],row['div'],expected),alignment_bp='',
                     transposon_mask_bp='',low_complexity_mask_bp='',conserved_mask_bp='')
            if row.get('cvg-by-all','') not in ('','n/a'):
                rec['alignment_bp']=int(row['cvg-by-all'])
                if not 0<=rec['alignment_bp']<=hi-lo:raise ValueError('Invalid GX alignment coverage')
            mask=next((v for k,v in row.items() if k.startswith('(xp,')),None)
            if mask:
                for key,v in zip(('transposon_mask_bp','low_complexity_mask_bp','conserved_mask_bp'),mask.split(',')):
                    rec[key]=int(v)
                    if not 0<=rec[key]<=hi-lo:raise ValueError('Invalid GX masking track')
            spans[name].append(rec)
        byname={r['sequence']:r for r in current}
        for name,length in lengths.items():
            end=0;partition=[]
            for rec in sorted(spans[name],key=lambda x:x['start0']):
                if rec['start0']<end:raise ValueError('Overlapping GX spans would double-count bases: '+name)
                if rec['start0']>end:partition.append(dict(cohort_assembly=assembly,coordinate_stage=stage,sequence=name,
                    start0=end,end0=rec['start0'],span_bp=rec['start0']-end,result='not reported',division='',
                    expected_division=expected,category='unreported'))
                partition.append(rec);end=rec['end0']
            if end<length:partition.append(dict(cohort_assembly=assembly,coordinate_stage=stage,sequence=name,
                start0=end,end0=length,span_bp=length-end,result='not reported',division='',
                expected_division=expected,category='unreported'))
            for rec in partition:
                byname[name]['cat_'+rec['category']+'_bp']+=rec['span_bp']
                rec['is_mixed']=byname[name]['is_mixed']
                segments.append(rec)
            byname[name]['expected_division']=expected
            if sum(byname[name]['cat_'+c+'_bp'] for c in CATEGORIES)!=length:raise ValueError('Incomplete category partition')
    fields=list(dict.fromkeys(k for r in segments for k in r))
    write(detail/'cohort_fcs_segments.tsv',segments,fields)
    for index,assembly in enumerate(sorted({r['cohort_assembly'] for r in sequences}),1):
        prefix=f'assembly_{index:04d}'
        current=[r for r in sequences if r['cohort_assembly']==assembly]
        selected=sorted((r for r in current if r['is_mixed']),key=lambda r:(-int(r['length_bp']),r['sequence']))[:30]
        write(detail/(prefix+'_selected.tsv'),selected,list(current[0]))
        (detail/(prefix+'.md')).write_text('\n'.join([
            '# '+assembly+' — '+stage+' FCS evidence','',
            '![GC/read-depth category pies]('+prefix+'_gc_coverage.png)','',
            'One pie per record; area represents record length and wedges represent category base-pair fractions. '
            'GC and depth are for the complete record. Repeat/low GX coverage remain inconclusive.', '',
            '![Mixed-record composition]('+prefix+'_composition.png)','',
            f'Showing {len(selected)} of {sum(r["is_mixed"] for r in current)} mixed records: the 30 longest, ties by identifier. '
            'This is a display limit, not a biological filter.', '',
            '![GX alignment and masking]('+prefix+'_diagnostics.png)','',
            'Diagnostics summarize reported segments within mixed records. Masking tracks can overlap; do not add their percentages. '
            'GX masking is not a separate repeat annotation. Taxonomy flags are not cleaning recommendations.', '',
            '[All category/coordinate evidence](cohort_fcs_segments.tsv) · '
            '[All record evidence](cohort_fcs_sequence_evidence.tsv) · '
            '[Displayed records]('+prefix+'_selected.tsv)','']),encoding='utf-8')
    (detail/'fcs_category_methods.json').write_text(json.dumps(dict(schema='fcs-display-v1',
        categories=CATEGORIES,expected_division='GX run-info.asserted-div; no inferred fallback',
        coordinates='0-based half-open on exact screened record; nested suffixes composed',
        composition_denominator='Full mixed-record span, including unreported intervals',
        diagnostics_denominator='Reported GX segment span within mixed records',
        display_limit='30 longest mixed records; ties by identifier',
        biological_rules_changed=False),indent=2)+'\n')
    return segments
