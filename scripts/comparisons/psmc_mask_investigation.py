"""Removable Bash/Slurm investigation of retained Dipcall inputs; no production tuning."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tarfile

REPO=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(REPO/'py_scripts'))
from psmc_input import bed

STATES=['no_PAF_coverage','secondary_only','incomplete_scoring_only','below_10kb_only',
        'below_MAPQ5_only','overlapping_segments_one_query','overlapping_multiple_queries',
        'unique_10_to_50kb','unique_long_alignment']
LABELS=['No alignment span','Secondary alignments only','Missing scoring tag',
        'Only alignments shorter than 10 kb','Only low-MAPQ alignments at least 10 kb',
        'Overlapping segments from one query','Overlapping alignments from multiple queries',
        'One qualifying 10-50 kb alignment','One qualifying alignment at least 50 kb']


def read_json(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def write_json(path,value):Path(path).write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8')
def required(path):
    path=Path(path)
    if not path.is_file():raise ValueError('Required input missing: '+str(path))
    return path.resolve()
def digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()
def table(path,rows,fields):
    with Path(path).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
def load_bed(path):
    rows=defaultdict(list)
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as f:
        for line in f:
            if not line.strip() or line.startswith(('#','track','browser')):continue
            n,lo,hi,*_=line.split();rows[n].append((int(lo),int(hi)))
    return {n:merge(v) for n,v in rows.items()}
def merge(rows):
    out=[]
    for lo,hi in sorted(rows):
        if hi<=lo:raise ValueError('Invalid interval')
        if out and lo<=out[-1][1]:out[-1]=(out[-1][0],max(hi,out[-1][1]))
        else:out.append((lo,hi))
    return out
def intersect(a,b):
    out=[];i=j=0
    while i<len(a) and j<len(b):
        lo=max(a[i][0],b[j][0]);hi=min(a[i][1],b[j][1])
        if hi>lo:out.append((lo,hi))
        if a[i][1]<b[j][1]:i+=1
        else:j+=1
    return out
def size(rows):return sum(hi-lo for lo,hi in rows)


def task_name(job):
    return '.'.join([job['sample']]+([job['preset']] if job.get('preset') else [])+[job['haplotype']])


def paired_counts(scope,h1,h2,baseline):
    q=Counter()
    for name,regions in scope.items():
        left=intersect(regions,h1.get(name,[]));right=intersect(regions,h2.get(name,[]));both=intersect(left,right)
        old=intersect(regions,baseline.get(name,[]));common=intersect(both,old)
        q.update(scope_bp=size(regions),both_confident_bp=size(both),hap1_only_bp=size(left)-size(both),
            hap2_only_bp=size(right)-size(both),neither_bp=size(regions)-size(left)-size(right)+size(both),
            original_confident_bp=size(old),reproduced_only_bp=size(both)-size(common),original_only_bp=size(old)-size(common))
    return q


def prepare_presets(a):
    source=Path(a.source).resolve();m=read_json(required(source/'manifest.json'))
    if read_json(required(source/'reports/status.json'))['status']!='PASS':raise ValueError('Source mask investigation did not pass')
    out=Path(a.out).resolve()
    if out.exists() and any(out.iterdir()):raise ValueError('Use a new preset-comparison directory')
    selected=sorted((j for j in m['tasks'] if j['sample']==a.sample),key=lambda j:j['hap'])
    if [j['hap'] for j in selected]!=[1,2]:raise ValueError('Expected two completed haplotypes for selected individual')
    if len(set(a.presets))!=len(a.presets) or not a.presets:raise ValueError('Specify distinct trial presets')
    controls=[];trials=[]
    # Validate sources before creating or submitting any trial tasks.
    for j in selected:
        root=source/'tasks'/task_name(j)
        q=read_json(required(root/'audit.json'))
        if q['official_only_bp'] or q['predicted_only_bp']:raise ValueError('Source mask reconstruction disagrees')
        if q.get('preset','asm5')!='asm5':raise ValueError('Expected asm5 control')
        for p in [j['reference']['fasta'],j['reference']['fai'],*j['assemblies'],root/'raw'/f"diploid.{j['haplotype']}.paf.gz"]:required(p)
        for name in ['scope.bed.gz','haplotype.bed','classification.bed.gz','contig_audit.tsv','alignment_records.tsv','commands.json']:required(root/name)
        controls.append(dict(j,preset='asm5',source_directory=str(root),
            source_paf=str(required(root/'raw'/f"diploid.{j['haplotype']}.paf.gz"))))
    out.mkdir(parents=True,exist_ok=True)
    for c in controls:
        src=Path(c['source_directory']);dest=out/'tasks'/task_name(c);dest.mkdir(parents=True)
        for name in ['scope.bed.gz','haplotype.bed','classification.bed.gz','contig_audit.tsv','alignment_records.tsv','commands.json']:
            shutil.copyfile(src/name,dest/name)
        q=read_json(src/'audit.json');q.update(preset='asm5',source_control=str(src),source_control_audit_sha256=digest(src/'audit.json'))
        write_json(dest/'audit.json',q)
    for preset in a.presets:
        if preset not in ['asm10','asm20']:raise ValueError('Trials must use asm10 or asm20')
        for j in selected:trials.append(dict(j,index=len(trials),preset=preset,cached_paf=None))
    result=dict(schema_version=1,mode='alignment_presets',repository=str(REPO),out=str(out),
        baseline=m['baseline'],plot_env=m['plot_env'],sample=a.sample,source=str(source),
        controls=controls,tasks=trials)
    write_json(out/'manifest.json',result)
    table(out/'input_inventory.tsv',[dict(sample=j['sample'],preset=j['preset'],haplotype=j['haplotype'],
        role='reused_control' if j['preset']=='asm5' else 'new_alignment',env=j['env']) for j in controls+trials],
        ['sample','preset','haplotype','role','env'])
    print(f"Prepared {len(trials)} new alignment tasks for {a.sample}; two asm5 controls reused without remapping.")


def make_override(makefile,preset):
    if preset not in ['asm5','asm10','asm20']:raise ValueError('Unsupported assembly preset')
    if preset=='asm5':return []
    options=re.findall(r'^MM2_OPT=(.*)$',Path(makefile).read_text(),re.M)
    if len(options)!=1:raise ValueError('Expected one generated Dipcall MM2_OPT assignment')
    changed,n=re.subn(r'(?<!\S)-xasm5(?!\S)','-x'+preset,options[0])
    if n!=1:raise ValueError('Cannot isolate the standard asm5 preset in generated Dipcall options')
    return ['MM2_OPT='+changed]


def alignment_quality(paf,out):
    groups=defaultdict(lambda:Counter())
    with gzip.open(paf,'rt') as f:
        for line in f:
            if not line.strip():continue
            fields=line.rstrip().split('\t')
            if len(fields)<12:raise ValueError('Malformed PAF')
            if fields[5]=='*' or paf_level(fields)<2 or int(fields[11])<5:continue
            block=int(fields[10]);matches=int(fields[9])
            if not 0<=matches<=block or block==0:raise ValueError('Invalid PAF identity counts')
            label='<10kb' if block<10000 else '10-50kb' if block<50000 else '>=50kb'
            g=groups[label];g.update(records=1,block_bp=block,matching_bp=matches)
            values=[float(t[5:]) for t in fields[12:] if t.startswith('de:f:')]
            if values:
                if len(values)!=1 or not math.isfinite(values[0]) or not 0<=values[0]<=1:raise ValueError('Invalid PAF gap-compressed divergence')
                g['divergence_weighted_sum']+=values[0]*block;g['divergence_weight_bp']+=block
            else:g['divergence_missing_records']+=1
    rows=[]
    for label in ['<10kb','10-50kb','>=50kb']:
        g=groups[label]
        rows.append(dict(block_length=label,records=g['records'],block_bp=g['block_bp'],matching_bp=g['matching_bp'],
            weighted_identity=g['matching_bp']/g['block_bp'] if g['block_bp'] else '',
            weighted_gap_compressed_divergence=g['divergence_weighted_sum']/g['divergence_weight_bp'] if g['divergence_weight_bp'] else '',
            divergence_missing_records=g['divergence_missing_records']))
    table(out,rows,list(rows[0]));return rows


def prepare(a):
    run=Path(a.run).resolve();base=run/'results';out=Path(a.out).resolve()
    if out.exists() and any(out.iterdir()):raise ValueError('Use a new investigation output directory: '+str(out))
    plan=read_json(required(base/'post_assembly/psmc/psmc_plan.json'))
    with required(run/'trace.tsv').open() as f:trace=list(csv.DictReader(f,delimiter='\t'))
    contexts=read_json(required(run/'execution_environments.json'))['tasks']
    def environment(task,override):
        if override:prefix=Path(override).resolve()
        else:
            found={p for row in contexts if row['workdir']==task['workdir'] for p in row['conda_prefixes']}
            if len(found)!=1:raise ValueError('Expected one recorded environment for '+task['name'])
            prefix=Path(found.pop())
        required(prefix/'bin/python');return str(prefix)
    plots=[r for r in trace if r['name'].endswith('PSMC_INDEX')]
    if len(plots)!=1:raise ValueError('Missing/ambiguous current-run PSMC_INDEX task')
    plot_env=environment(plots[0],a.plot_env);required(Path(plot_env)/'bin/Rscript')
    tasks=[];jobs=[j for j in plan['jobs'] if j['mode']=='assemblies']
    if not jobs:raise ValueError('No assembly-based Dipcall jobs in this run')
    for j in sorted(jobs,key=lambda j:j['sample']):
        sample=j['sample']
        if not re.fullmatch(r'[A-Za-z0-9_.-]+',sample):raise ValueError('Unsafe sample ID')
        candidates=[r for r in trace if r['name'].endswith(':PSMC_ASSEMBLY_CALL ('+sample+')')]
        if len(candidates)!=1:raise ValueError('Missing/ambiguous Dipcall trace row: '+sample)
        task=candidates[0];prefix=environment(task,a.analysis_env)
        if task.get('status') not in ['COMPLETED','CACHED'] or task.get('exit')!='0':
            raise ValueError('Dipcall trace task did not succeed: '+sample)
        if len(j['assemblies'])!=2:raise ValueError('Expected two same-individual assemblies')
        sources=[str(required(p)) for p in j['assemblies']]
        ref=j['reference'];required(ref['fasta']);required(ref['fai'])
        if ref['sample']!=sample or Path(ref['fasta']).resolve()!=Path(sources[0]):raise ValueError('Expected own hap1 as reference-to-itself control')
        evidence=base/'post_assembly/psmc'/str(j['taxid'])/sample
        filters=read_json(required(evidence/'input/provenance.json'))['filters']
        qc=read_json(required(evidence/'input/input_qc.json'))
        if qc['scope_bp']<=0:raise ValueError('No selected non-gap reference scope: '+sample)
        if qc['reference_sha256']!=ref['sha256']:raise ValueError('Baseline reference and PSMC input fingerprints disagree')
        for path in [j.get('include_bed'),j.get('exclude_bed')]:
            if path:required(path)
        original_bed=str(required(evidence/'genotypes/diploid.dip.bed'))
        for hap in [1,2]:
            paf=Path(task['workdir'])/f'diploid.hap{hap}.paf.gz'
            tasks.append(dict(index=len(tasks),sample=sample,taxid=str(j['taxid']),haplotype=f'hap{hap}',
                hap=hap,reference=ref,assemblies=sources,env=prefix,
                cached_paf=str(paf.resolve()) if paf.is_file() else None,
                original_task=task['workdir'],original_diploid_bed=original_bed,
                filters=filters,expected_scope_bp=qc['scope_bp'],
                include_bed=j.get('include_bed'),exclude_bed=j.get('exclude_bed')))
    selected=[j for j in jobs if j['sample']==a.hifi_sample]
    if len(selected)!=1 or not selected[0].get('reads'):raise ValueError('Comparator requires a selected individual with complete HiFi reads')
    required(selected[0]['reads'])
    out.mkdir(parents=True,exist_ok=True)
    manifest=dict(schema_version=1,repository=str(REPO),baseline=str(base),out=str(out),plot_env=plot_env,
        hifi_sample=a.hifi_sample,tasks=tasks)
    write_json(out/'manifest.json',manifest)
    table(out/'hifi_manifest.tsv',[dict(sample=a.hifi_sample,mode='hifi')],['sample','mode'])
    table(out/'input_inventory.tsv',[dict(index=t['index'],sample=t['sample'],haplotype=t['haplotype'],
        alignment_source='cached' if t['cached_paf'] else 'regenerate_PAF_only',original_task=t['original_task'],
        env=t['env']) for t in tasks],['index','sample','haplotype','alignment_source','original_task','env'])
    print(f"Prepared {len(tasks)} independent haplotype tasks; {sum(bool(t['cached_paf']) for t in tasks)} cached PAFs available; the others regenerate PAF/masks only.")


def paf_level(fields):
    """Nested filter tiers from paftools.js 2.28 call, using PAF column 11 length."""
    tags={t[:2]:t[5:] for t in fields[12:] if t.startswith(('tp:A:','s1:i:','s2:i:'))}
    if tags.get('tp') in ['S','i']:return 0
    if 's1' in tags and 's2' not in tags:return 1
    if int(fields[10])<10000:return 2
    if int(fields[11])<5:return 3
    return 5 if int(fields[10])>=50000 else 4


def reason(counts,queries):
    if not counts[0]:return STATES[0]
    if not counts[1]:return STATES[1]
    if not counts[2]:return STATES[2]
    if not counts[3]:return STATES[3]
    if not counts[4]:return STATES[4]
    if counts[4]>1:return STATES[6] if len(queries)>1 else STATES[5]
    return STATES[8] if counts[5] else STATES[7]


def sweep(events):
    counts=[0]*6;queries={};scope=caller=0;previous=None
    for pos,group in itertools.groupby(sorted(events,key=lambda e:e[0]),lambda e:e[0]):
        if previous is not None and pos>previous and scope:
            yield previous,pos,reason(counts,queries),bool(caller)
        for _,delta,level,query in group:
            if level==-1:scope+=delta
            elif level==-2:caller+=delta
            else:
                for k in range(level+1):counts[k]+=delta
                if level>=4:
                    queries[query]=queries.get(query,0)+delta
                    if queries[query]==0:del queries[query]
        if scope not in [0,1] or caller<0 or min(counts)<0:raise ValueError('Invalid interval coverage state')
        previous=pos
    if scope or caller or any(counts) or queries:raise ValueError('Unbalanced interval endpoints')


def reference_scope(fasta,minimum,include,exclude):
    import numpy as np
    lengths=dict(zip(fasta.references,fasta.lengths))
    inclusion=bed(include,lengths);exclusion=bed(exclude,lengths);result={}
    for name,length in lengths.items():
        if length<minimum:continue
        pieces=[]
        for lo in range(0,length,2000000):
            seq=np.frombuffer(fasta.fetch(name,lo,min(length,lo+2000000)).upper().encode(),dtype='S1')
            known=np.isin(seq,[b'A',b'C',b'G',b'T'])
            changes=np.flatnonzero(np.diff(np.concatenate(([0],known.astype('int8'),[0]))))
            pieces.extend((lo+int(x),lo+int(y)) for x,y in zip(changes[::2],changes[1::2]))
        pieces=merge(pieces)
        if inclusion is not None:pieces=intersect(pieces,inclusion[name])
        if exclusion is not None:
            allowed=[];start=0
            for lo,hi in exclusion[name]:
                if start<lo:allowed.append((start,lo))
                start=hi
            if start<length:allowed.append((start,length))
            pieces=intersect(pieces,allowed)
        result[name]=pieces
    return result,lengths


def analyze(paf,caller_bed,fasta,job,out):
    scope,lengths=reference_scope(fasta,job['filters']['min_contig_bp'],job.get('include_bed'),job.get('exclude_bed'))
    if sum(size(v) for v in scope.values())!=job['expected_scope_bp']:raise ValueError('Diagnostic scope disagrees with original PSMC input')
    events=defaultdict(list);records=Counter();hist=Counter()
    with gzip.open(paf,'rt') as f:
        for line in f:
            if not line.strip():continue
            fields=line.rstrip().split('\t')
            if len(fields)<12:raise ValueError('Malformed PAF record')
            if fields[5]=='*':records['unmapped_records']+=1;continue
            name=fields[5];lo,hi=int(fields[7]),int(fields[8])
            if name not in lengths or int(fields[6])!=lengths[name] or not 0<=lo<hi<=lengths[name]:raise ValueError('PAF/reference dictionary disagreement')
            level=paf_level(fields);records['mapped_records']+=1
            hist[(level,'<10kb' if int(fields[10])<10000 else '10-50kb' if int(fields[10])<50000 else '>=50kb')]+=1
            if name in scope:
                events[name].extend([(lo,1,level,fields[0]),(hi,-1,level,fields[0])])
    caller=bed(caller_bed,lengths);totals=Counter();rows=[];mismatch=Counter()
    with gzip.open(out/'scope.bed.gz','wt') as sc,gzip.open(out/'classification.bed.gz','wt') as intervals:
        for name,pieces in scope.items():
            for lo,hi in pieces:
                sc.write(f'{name}\t{lo}\t{hi}\n');events[name].extend([(lo,1,-1,''),(hi,-1,-1,'')])
            for lo,hi in caller[name]:events[name].extend([(lo,1,-2,''),(hi,-1,-2,'')])
            counts=Counter();local=Counter()
            for lo,hi,state,official in sweep(events.pop(name,[])):
                bp=hi-lo;counts[state]+=bp
                if official:local['official_confident_bp']+=bp
                if state==STATES[-1]:local['predicted_confident_bp']+=bp
                if official and state!=STATES[-1]:local['official_only_bp']+=bp
                if not official and state==STATES[-1]:local['predicted_only_bp']+=bp
                intervals.write(f'{name}\t{lo}\t{hi}\t{state}\t{int(official)}\n')
            totals.update(counts);mismatch.update(local)
            rows.append(dict(contig=name,scope_bp=size(pieces),**{k:counts[k] for k in STATES},
                **{k:local[k] for k in ['official_confident_bp','predicted_confident_bp','official_only_bp','predicted_only_bp']}))
    table(out/'contig_audit.tsv',rows,list(rows[0]))
    table(out/'alignment_records.tsv',[dict(filter_tier=k[0],block_length=k[1],records=v) for k,v in sorted(hist.items())],['filter_tier','block_length','records'])
    result=dict(sample=job['sample'],haplotype=job['haplotype'],reference_sha256=job['reference']['sha256'],
        scope_bp=sum(totals.values()),classification_bp=dict(totals),alignment_records=dict(records),
        thresholds=dict(min_coverage_block_bp=10000,min_variant_block_bp=50000,min_mapq=5),
        alignment_record_scope='All mapped records, including reference sequences below the PSMC minimum length',
        **{k:mismatch[k] for k in ['official_confident_bp','predicted_confident_bp','official_only_bp','predicted_only_bp']})
    write_json(out/'analysis.json',result);return result


def run_task(a):
    import pysam
    manifest=read_json(a.manifest);job=manifest['tasks'][a.index]
    out=Path(manifest['out'])/'tasks'/task_name(job);out.mkdir(parents=True,exist_ok=True)
    spec=out/'task.json'
    if spec.exists() and read_json(spec)!=job:raise ValueError('Existing task directory belongs to different inputs')
    write_json(spec,job)
    if digest(job['reference']['fasta'])!=job['reference']['sha256']:raise ValueError('Reference SHA256 changed since baseline')
    work=out/'raw';work.mkdir(exist_ok=True)
    for name,source in [('reference.fa',job['reference']['fasta']),('reference.fa.fai',job['reference']['fai']),
                        ('hap1.fa',job['assemblies'][0]),('hap2.fa',job['assemblies'][1])]:
        dest=work/name
        if not dest.exists():dest.symlink_to(required(source))
    env=dict(os.environ,PATH=job['env']+'/bin:'+os.environ['PATH'],LC_ALL='C')
    commands=[]
    def run(cmd,stdout=None):
        commands.append(list(map(str,cmd)))
        with (out/'commands.log').open('a') as log:
            subprocess.run(cmd,cwd=work,env=env,stdout=stdout or log,stderr=log,check=True)
    with (work/'diploid.mak').open('w') as f:
        run(['run-dipcall','-a','-t',str(a.cpus),'diploid','reference.fa','hap1.fa','hap2.fa'],f)
    name=f"diploid.{job['haplotype']}";paf=work/(name+'.paf.gz')
    preset=job.get('preset','asm5')
    if preset!='asm5' and job['cached_paf']:raise ValueError('A new preset must not reuse an asm5 PAF')
    if job['cached_paf'] and not paf.exists():paf.symlink_to(required(job['cached_paf']))
    cmd=['make','-f','diploid.mak']
    cmd+=make_override(work/'diploid.mak',preset)
    if job['cached_paf']:cmd+=['-o',name+'.paf.gz']
    run(cmd+[name+'.bed'])
    shutil.copyfile(work/(name+'.bed'),out/'haplotype.bed')
    with pysam.FastaFile(str(work/'reference.fa')) as fa:result=analyze(paf,out/'haplotype.bed',fa,job,out)
    alignment_quality(paf,out/'alignment_quality.tsv')
    result.update(preset=preset,alignment_source='cached' if job['cached_paf'] else 'regenerated_PAF_only',
        diagnostic_sha256=digest(__file__),scope_helper_sha256=digest(REPO/'py_scripts/psmc_input.py'),
        query_sha256=digest(job['assemblies'][job['hap']-1]),paf_sha256=digest(paf),
        minimap2_version=subprocess.check_output(['minimap2','--version'],env=env,text=True).strip(),
        paftools_sha256=digest(Path(job['env'])/'bin/paftools.js'),
        run_dipcall_sha256=digest(Path(job['env'])/'bin/run-dipcall'))
    write_json(out/'audit.json',result);write_json(out/'commands.json',commands)
    print(f"Completed {job['sample']} {job['haplotype']}; mask disagreement {result['official_only_bp']+result['predicted_only_bp']} bp")


def report(a):
    m=read_json(a.manifest);out=Path(m['out']);reports=out/'reports';reports.mkdir(exist_ok=True)
    summaries=[];state_rows=[];pair_rows=[];missing=[]
    for job in m['tasks']:
        root=out/'tasks'/f"{job['sample']}.{job['haplotype']}"
        if not (root/'audit.json').is_file():missing.append(root.name);continue
        q=read_json(root/'audit.json');summaries.append(q)
        state_rows.extend(dict(sample=q['sample'],haplotype=q['haplotype'],state=s,label=l,bp=q['classification_bp'].get(s,0),
            percent=100*q['classification_bp'].get(s,0)/q['scope_bp']) for s,l in zip(STATES,LABELS))
    for sample in sorted({j['sample'] for j in m['tasks']}):
        roots=[out/'tasks'/f'{sample}.hap{i}' for i in [1,2]]
        if any(not (r/'audit.json').is_file() for r in roots):continue
        scope=load_bed(roots[0]/'scope.bed.gz');h1=load_bed(roots[0]/'haplotype.bed');h2=load_bed(roots[1]/'haplotype.bed')
        job=next(j for j in m['tasks'] if j['sample']==sample);original=load_bed(job['original_diploid_bed'])
        q=Counter()
        for name,regions in scope.items():
            left=intersect(regions,h1.get(name,[]));right=intersect(regions,h2.get(name,[]));both=intersect(left,right)
            old=intersect(regions,original.get(name,[]));common=intersect(both,old)
            q.update(scope_bp=size(regions),both_confident_bp=size(both),hap1_only_bp=size(left)-size(both),
                hap2_only_bp=size(right)-size(both),neither_bp=size(regions)-size(left)-size(right)+size(both),
                original_confident_bp=size(old),reproduced_only_bp=size(both)-size(common),original_only_bp=size(old)-size(common))
        pair_rows.append(dict(sample=sample,**q))
    table(reports/'classification_summary.tsv',state_rows,['sample','haplotype','state','label','bp','percent'])
    table(reports/'paired_masks.tsv',pair_rows,['sample','scope_bp','both_confident_bp','hap1_only_bp','hap2_only_bp','neither_bp',
        'original_confident_bp','reproduced_only_bp','original_only_bp'])
    disagreement=sum(q['official_only_bp']+q['predicted_only_bp'] for q in summaries)+sum(q['reproduced_only_bp']+q['original_only_bp'] for q in pair_rows)
    status='INCOMPLETE' if missing else 'MASK_REPRODUCTION_MISMATCH' if disagreement else 'PASS'
    write_json(reports/'status.json',dict(status=status,missing=missing,mask_disagreement_bp=disagreement,summary=summaries))
    md=['# Dipcall alignment and mask investigation','',f'Status: **{status}**.','',
        'A PASS verifies task completion and mask reproduction. Biological interpretation follows from the loss categories.','',
        'Diagnostic only. Production assemblies, alignments and confidence thresholds are unchanged. Missing raw PAFs are regenerated with the recorded environment and the standard Dipcall PAF target; SAM/BAM/VCF calling is not rerun.','',
        '![Alignment exclusions by haplotype](reports/mask_reasons.png)','',
        '| Individual | Hap1 mask, Mb | Hap2 mask, Mb | Both confident, Mb | Neither confident, Mb | Original-mask discrepancy, bp |',
        '| --- | ---: | ---: | ---: | ---: | ---: |']
    for r in pair_rows:
        h1=(r['both_confident_bp']+r['hap1_only_bp'])/1e6;h2=(r['both_confident_bp']+r['hap2_only_bp'])/1e6
        md.append(f"| {r['sample']} | {h1:.2f} | {h2:.2f} | {r['both_confident_bp']/1e6:.2f} | {r['neither_bp']/1e6:.2f} | {r['reproduced_only_bp']+r['original_only_bp']} |")
    md+=['','[Paired mask counts](reports/paired_masks.tsv) · [All exclusion counts](reports/classification_summary.tsv) · [Status and provenance](reports/status.json)','',
        'Hap1 is aligned to its own reference sequence. Its exclusions are an important control. Categories are a hierarchical explanation of observed PAF coverage, not proof of repeats, copy-number differences or biological absence. Only primary/inversion-primary, scoring-compatible records enter the 10 kb / MAPQ 5 coverage calculation; a single eligible alignment must also reach 50 kb to support the standard mask. PAF alignment-block length is used, not target span.','',
        'Overlap exclusions distinguish segments from the same query from alignments of different query sequences. The predicted single-long-alignment mask is checked against actual paftools output. The intersection of the two actual haplotype masks is checked against the original published diploid BED. Any discrepancy is reported rather than explained away.','',
        f"The complete-HiFi comparison is prepared for **{m['hifi_sample']}** in `hifi_manifest.tsv`; run it through the production workflow alongside this diagnostic. No quality thresholds have been relaxed."]
    if missing:md+=['','Missing or failed tasks: '+', '.join(missing)]
    for q in summaries:
        root=f"tasks/{q['sample']}.{q['haplotype']}"
        md+=['',f"## {q['sample']} / {q['haplotype']}",'',f"Alignment source: {q['alignment_source']}. Scope: {q['scope_bp']:,} bp. Predicted/actual mask disagreement: {q['official_only_bp']+q['predicted_only_bp']:,} bp.",'',
            '| Exclusion or support | Mb | Scope % |','| --- | ---: | ---: |']
        for s,l in zip(STATES,LABELS):
            bp=q['classification_bp'].get(s,0);md.append(f'| {l} | {bp/1e6:.2f} | {100*bp/q["scope_bp"]:.2f} |')
        md+=['',f'[Per-contig evidence]({root}/contig_audit.tsv) · [Classified intervals]({root}/classification.bed.gz) · [Commands]({root}/commands.json) · [Provenance]({root}/audit.json)']
    (out/'README.md').write_text('\n'.join(md)+'\n',encoding='utf-8')
    if summaries:
        subprocess.run([str(Path(m['plot_env'])/'bin/Rscript'),str(REPO/'r_scripts/plot_psmc_mask_investigation.R'),str(reports)],check=True)
    with tarfile.open(str(out)+'.review.tar.gz','w:gz') as t:
        for name in ['README.md','manifest.json','hifi_manifest.tsv','input_inventory.tsv','reports']:
            if (out/name).exists():t.add(out/name,arcname=name)
        for job in m['tasks']:
            root=out/'tasks'/f"{job['sample']}.{job['haplotype']}"
            if not root.exists():continue
            for p in root.iterdir():
                if p.is_file():t.add(p,arcname=str(p.relative_to(out)))
            for p in (root/'raw').glob('*.log'):
                if p.is_file():t.add(p,arcname=str(p.relative_to(out)))
            for p in (root/'raw').glob('*.vst'):
                if p.is_file():t.add(p,arcname=str(p.relative_to(out)))
    print('Report:',out/'README.md')
    if missing or disagreement:raise ValueError(f'Investigation {status}; inspect the report and task logs')


def report_presets(a):
    m=read_json(a.manifest);out=Path(m['out']);reports=out/'reports';reports.mkdir(exist_ok=True)
    jobs=m['controls']+m['tasks'];audits={};missing=[];issues=[];states=[];qualities=[];pairs=[]
    for j in jobs:
        root=out/'tasks'/task_name(j)
        if not (root/'audit.json').is_file():missing.append(task_name(j));continue
        q=read_json(root/'audit.json');audits[j['preset'],j['hap']]=q
        if q['official_only_bp'] or q['predicted_only_bp']:issues.append(task_name(j)+': mask reconstruction disagrees')
        control=read_json(out/'tasks'/task_name(next(c for c in m['controls'] if c['hap']==j['hap']))/'audit.json')
        for key in ['reference_sha256','query_sha256','scope_bp','thresholds','minimap2_version','paftools_sha256','run_dipcall_sha256','scope_helper_sha256']:
            if q.get(key)!=control.get(key):issues.append(task_name(j)+': changed '+key)
        paf=Path(j['source_paf']) if j['preset']=='asm5' else root/'raw'/f"diploid.{j['haplotype']}.paf.gz"
        if digest(required(paf))!=q['paf_sha256']:raise ValueError('PAF changed since audit: '+str(paf))
        quality=alignment_quality(paf,root/'alignment_quality.tsv')
        qualities.extend(dict(preset=j['preset'],haplotype=j['haplotype'],**r) for r in quality)
        states.extend(dict(preset=j['preset'],haplotype=j['haplotype'],state=s,label=l,bp=q['classification_bp'].get(s,0),
            percent=100*q['classification_bp'].get(s,0)/q['scope_bp']) for s,l in zip(STATES,LABELS))
    controls=[out/'tasks'/task_name(c) for c in m['controls']]
    scope=load_bed(controls[0]/'scope.bed.gz')
    old1=load_bed(controls[0]/'haplotype.bed');old2=load_bed(controls[1]/'haplotype.bed')
    oldboth={n:intersect(old1.get(n,[]),old2.get(n,[])) for n in scope}
    baseline=paired_counts(scope,old1,old2,load_bed(m['controls'][0]['original_diploid_bed']))
    if baseline['reproduced_only_bp'] or baseline['original_only_bp']:issues.append('asm5 control no longer reproduces published mask')
    for preset in dict.fromkeys(j['preset'] for j in jobs):
        if any((preset,h) not in audits for h in [1,2]):continue
        roots=[out/'tasks'/task_name(next(j for j in jobs if j['preset']==preset and j['hap']==h)) for h in [1,2]]
        if load_bed(roots[0]/'scope.bed.gz')!=scope or load_bed(roots[1]/'scope.bed.gz')!=scope:
            issues.append(preset+': changed scope intervals');continue
        q=paired_counts(scope,load_bed(roots[0]/'haplotype.bed'),load_bed(roots[1]/'haplotype.bed'),oldboth)
        pairs.append(dict(preset=preset,scope_bp=q['scope_bp'],hap1_confident_bp=q['both_confident_bp']+q['hap1_only_bp'],
            hap2_confident_bp=q['both_confident_bp']+q['hap2_only_bp'],paired_confident_bp=q['both_confident_bp'],
            paired_percent=100*q['both_confident_bp']/q['scope_bp'],gained_vs_asm5_bp=q['reproduced_only_bp'],lost_vs_asm5_bp=q['original_only_bp']))
    table(reports/'paired_masks.tsv',pairs,['preset','scope_bp','hap1_confident_bp','hap2_confident_bp','paired_confident_bp','paired_percent','gained_vs_asm5_bp','lost_vs_asm5_bp'])
    table(reports/'classification_summary.tsv',states,['preset','haplotype','state','label','bp','percent'])
    table(reports/'alignment_quality.tsv',qualities,['preset','haplotype','block_length','records','block_bp','matching_bp','weighted_identity','weighted_gap_compressed_divergence','divergence_missing_records'])
    status='INCOMPLETE' if missing else 'VALIDATION_FAILED' if issues else 'PASS'
    write_json(reports/'status.json',dict(status=status,missing=missing,issues=issues,summary=list(audits.values())))
    md=['# Assembly alignment preset comparison','',f"Individual: **{m['sample']}**. Test validation: **{status}**.",'',
        'Only the minimap2 assembly preset changes. The 50 kb variant/mask requirement, 10 kb coverage requirement, MAPQ 5 and overlap rules remain unchanged. The two asm5 controls reuse the completed mask investigation. No production defaults or ongoing HiFi jobs are changed.','',
        '| Preset | Hap1 mask, Mb | Hap2 mask, Mb | Paired mask, Mb | Paired scope % | Gained vs asm5, Mb | Lost vs asm5, Mb |',
        '| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for q in pairs:md.append(f"| {q['preset']} | {q['hap1_confident_bp']/1e6:.2f} | {q['hap2_confident_bp']/1e6:.2f} | {q['paired_confident_bp']/1e6:.2f} | {q['paired_percent']:.2f} | {q['gained_vs_asm5_bp']/1e6:.2f} | {q['lost_vs_asm5_bp']/1e6:.2f} |")
    md+=['','![Alignment support and exclusions](reports/preset_comparison.png)','',
        'PASS means the tasks completed, inputs/tools matched the control and predicted masks reproduced actual paftools masks. It does not select a biologically validated winning preset. Gains and losses against asm5 are expected measurements, not validation errors.','',
        '[Mask counts](reports/paired_masks.tsv) · [Exclusion counts](reports/classification_summary.tsv) · [Alignment lengths and identity](reports/alignment_quality.tsv) · [Validation and provenance](reports/status.json)','',
        'Alignment identity and gap-compressed divergence are weighted by alignment-block length across primary, scoring-compatible records with MAPQ ≥5. These record-level summaries can count overlapping placements; they are not unique genome coverage or a genome-wide heterozygosity estimate. More retained sequence alone does not establish correct orthology.','',
        'Inspect whether hap2 short blocks become longer, whether newly retained sequence is supported by plausible identity/divergence, and whether ambiguity increases. These mask-only tests do not produce replacement VCFs or PSMC fits. If a preset is adopted later, the full Dipcall alignment and variant-calling route must be rerun consistently.']
    for j in jobs:
        name=task_name(j)
        md+=['',f"- {j['preset']} / {j['haplotype']}: [per-contig counts](tasks/{name}/contig_audit.tsv), [commands](tasks/{name}/commands.json), [audit](tasks/{name}/audit.json)."]
    if missing:md+=['','Missing tasks: '+', '.join(missing)]
    if issues:md+=['','Validation issues: '+ '; '.join(issues)]
    (out/'README.md').write_text('\n'.join(md)+'\n',encoding='utf-8')
    if states:subprocess.run([str(Path(m['plot_env'])/'bin/Rscript'),str(REPO/'r_scripts/plot_psmc_alignment_presets.R'),str(reports)],check=True)
    if missing or issues:raise ValueError('Preset investigation '+status+'; inspect reports/status.json')
    print('Report:',out/'README.md')


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('prepare')
    for k in ['run','out','hifi-sample']:q.add_argument('--'+k,required=True)
    q.add_argument('--analysis-env');q.add_argument('--plot-env')
    q=sub.add_parser('prepare-presets')
    for k in ['source','out','sample']:q.add_argument('--'+k,required=True)
    q.add_argument('--presets',nargs='+',choices=['asm10','asm20'],default=['asm10','asm20'])
    q=sub.add_parser('task');q.add_argument('--manifest',required=True);q.add_argument('--index',type=int,required=True);q.add_argument('--cpus',type=int,default=8)
    q=sub.add_parser('report');q.add_argument('--manifest',required=True)
    q=sub.add_parser('report-presets');q.add_argument('--manifest',required=True)
    a=p.parse_args();{'prepare':prepare,'prepare-presets':prepare_presets,'task':run_task,'report':report,'report-presets':report_presets}[a.command](a)


if __name__=='__main__':main()
