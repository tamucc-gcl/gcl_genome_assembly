"""Pinned production Cactus adapters: guarded filtering, regional assignment and source accounting."""
import argparse
from collections import defaultdict
import gzip
import importlib.metadata
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import uuid
import xml.etree.ElementTree as ET

from pangenome_cactus_guard import node_owners, lane_paf, restore_names, sha256
from pangenome_source_audit import read_table, write_table, union_bp

POLICY='chromosome_guarded_regional_v1'
CACTUS_COMMIT='d7e5fc1b51327a5a4cdb494c85b0607e3d81a713'
TOOLS_COMMIT='f9b999e54c49a05ad6e9d44e1832c4a8421f28f2'
OLD='if (interval.value.first != ref || clip_candidates.empty()) {'
NEW='if (interval.value.second != ref || clip_candidates.empty()) {'


def dump(path,data):Path(path).write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')


def build(source,out):
    source=Path(source);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    manifest=json.loads((source/'manifest.json').read_text())
    if manifest['commit']!=TOOLS_COMMIT:raise ValueError('Unexpected regional source commit')
    src=out/'source';src.mkdir()
    for name,digest in manifest['files'].items():
        if Path(name).name!=name or sha256(source/name)!=digest:raise ValueError('Regional source checksum mismatch')
        shutil.copyfile(source/name,src/name)
    text=(src/'rgfa-split.cpp').read_text()
    if text.count(OLD)!=1:raise ValueError('Regional correction context mismatch')
    (src/'rgfa-split.cpp').write_text(text.replace(OLD,NEW),encoding='utf-8')
    cmd=shlex.split(os.environ.get('CXX','c++'))+shlex.split(os.environ.get('CPPFLAGS',''))+shlex.split(os.environ.get('CXXFLAGS',''))+[
        '-O2','-std=c++14','-include','limits','-include','cstdint','-pthread','-static-libstdc++','-static-libgcc',
        '-I',str(src),str(src/'rgfa-split_main.cpp'),str(src/'rgfa-split.cpp'),str(src/'pafcoverage.cpp'),
        '-o',str(out/'rgfa-split')]+shlex.split(os.environ.get('LDFLAGS',''))
    with (out/'compile.log').open('w') as log:subprocess.run(cmd,check=True,stdout=log,stderr=log)
    # Real CLI build check: one coherent secondary run composed of short intervals.
    check=out/'check';check.mkdir();paf=check/'input.paf'
    blocks=[(i*50000,i*50000+40000,'A') for i in range(60)]+[(3000000,8000000,'B')]
    paf.write_text(''.join(f'source\t8000000\t{lo}\t{hi}\t+\t{chrom}\t10000000\t0\t{hi-lo}\t{hi-lo}\t{hi-lo}\t60\n' for lo,hi,chrom in blocks))
    args=[str((out/'rgfa-split').resolve()),'-p',str(paf),'-b',str(check/'split_'),'-n','.25','-u','25000','-s','-a','_AMBIGUOUS_']
    with (check/'log').open('w') as log:subprocess.run(args,check=True,stdout=log,stderr=log)
    names=(check/'split_A.fa_contigs').read_text().splitlines()
    if names.count('source:1-2990000')!=1:raise ValueError('Production regional native regression failed')
    dump(out/'build.json',dict(policy=POLICY,source_manifest=manifest,binary_sha256=sha256(out/'rgfa-split'),
        command=cmd,patched_source_sha256=sha256(src/'rgfa-split.cpp'),native_regression='PASS',correction=dict(before=OLD,after=NEW)))


def fasta_lengths(path):
    result={};name=None;opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt',encoding='ascii') as f:
        for line in f:
            if line.startswith('>'):
                name=line[1:].split()[0]
                if name in result:raise ValueError('Duplicate production source FASTA ID')
                result[name]=0
            elif line.strip():
                if name is None:raise ValueError('Sequence before FASTA header')
                result[name]+=len(''.join(line.split()))
    if not result or any(n<=0 for n in result.values()):raise ValueError('Empty source sequence')
    return result


def prepare(tools,seqfile,runtime,out,chunk):
    import cactus
    from cactus.shared.version import cactus_commit
    if importlib.metadata.version('cactus')!='3.2.1' or cactus_commit!=CACTUS_COMMIT:
        raise ValueError('Production sequence policy requires audited Cactus 3.2.1 commit')
    if chunk<=0:raise ValueError('Regional chunk size must be positive')
    tools=Path(tools).resolve();runtime=Path(runtime).resolve();out=Path(out).resolve()
    runtime.mkdir();out.mkdir(parents=True,exist_ok=True);(runtime/'bin').mkdir()
    info=json.loads((tools/'build.json').read_text())
    if info['policy']!=POLICY or info['native_regression']!='PASS' or sha256(tools/'rgfa-split')!=info['binary_sha256']:
        raise ValueError('Production regional executable provenance mismatch')
    native=shutil.which('filter-paf-deletions')
    if not native:raise ValueError('Native deletion filter unavailable')
    config=Path(cactus.__file__).parent/'cactus_progressive_config.xml';xml=ET.parse(config).getroot()
    if xml.find('graphmap').get('collapse','none')!='none' or int(xml.find('graphmap').get('delFilter','-1'))<=0:
        raise ValueError('Production guarded policy requires positive deletion filtering without collapse')
    sources=[];seen=set()
    for line in Path(seqfile).read_text().splitlines():
        sample,fa=line.split()
        for sequence,length in fasta_lengths(fa).items():
            query=f'id={sample}|{sequence}'
            if query in seen:raise ValueError('Duplicate production source identity')
            seen.add(query);sources.append(dict(query=query,sample=sample,sequence=sequence,source_bp=length))
    write_table(out/'source_lengths.tsv',sources,['query','sample','sequence','source_bp'])
    policy=dict(policy=POLICY,cactus_commit=cactus_commit,chunk_bp=chunk,softclip=True,
        corrected_binary=str(tools/'rgfa-split'),binary_sha256=info['binary_sha256'],
        native_filter=str(Path(native).resolve()),native_filter_sha256=sha256(native),
        source_lengths=str(out/'source_lengths.tsv'),evidence=str(out),
        ambiguous_name=xml.find('graphmap_split').get('ambiguousName','_AMBIGUOUS_'),
        config_sha256=sha256(config),runtime_sha256=sha256(__file__))
    dump(runtime/'policy.json',policy);dump(out/'policy.json',policy);shutil.copyfile(tools/'build.json',out/'tool_build.json')
    for name,mode in [('filter-paf-deletions','filter'),('rgfa-split','split')]:
        script=(f'#!{sys.executable}\nimport os,sys\nos.execv({sys.executable!r},'
                f'[{sys.executable!r},{str(Path(__file__).resolve())!r},{mode!r},'
                f'"--policy",{str(runtime/"policy.json")!r},"--"]+sys.argv[1:])\n')
        file=runtime/'bin'/name;file.write_text(script);file.chmod(0o755)


def load_policy(path):
    data=json.loads(Path(path).read_text())
    if data['policy']!=POLICY or sha256(data['corrected_binary'])!=data['binary_sha256']:
        raise ValueError('Production policy/binary mismatch')
    return data


def guarded_filter(policy,args):
    data=load_policy(policy)
    if len(args)<3:raise ValueError('Unsupported native deletion filter invocation')
    vg,trans,paf=map(Path,args[:3]);gfa=Path(str(vg).removesuffix('.vg'))
    if str(gfa)==str(vg) or not gfa.is_file():raise ValueError('Pinned deletion filter GFA companion unavailable')
    if sha256(data['native_filter'])!=data['native_filter_sha256']:raise ValueError('Native filter checksum mismatch')
    with tempfile.TemporaryDirectory(prefix='gcl-filter-',dir=paf.parent) as directory:
        temp=Path(directory);node_owners(gfa,temp/'owners.tsv.gz')
        names=lane_paf(paf,temp/'owners.tsv.gz',temp/'lanes.paf',temp/'lanes.tsv')
        cmd=[data['native_filter'],str(vg),str(trans),str(temp/'lanes.paf'),*args[3:]]
        with (temp/'filtered.paf').open('wb') as output:subprocess.run(cmd,stdout=output,check=True)
        restore_names(temp/'filtered.paf',temp/'restored.paf',names)
        ident='filter_'+uuid.uuid4().hex
        with open(temp/'lanes.tsv','rb') as src,gzip.open(Path(data['evidence'])/(ident+'.lanes.tsv.gz'),'wb') as dst:
            shutil.copyfileobj(src,dst)
        dump(Path(data['evidence'])/(ident+'.json'),dict(kind='filter',policy=POLICY,lanes=len(names),
            input_sha256=sha256(paf),output_sha256=sha256(temp/'restored.paf'),native_command=cmd))
        with (temp/'restored.paf').open('rb') as src:shutil.copyfileobj(src,sys.stdout.buffer)


def regional_command(data,args):
    result=[];i=0
    while i<len(args):
        if args[i] in ('-u','--min-query-chunk'):
            if i+1==len(args):raise ValueError('Missing chunk option value')
            i+=2
        elif args[i] in ('-s','--allow-softclip'):i+=1
        else:result.append(args[i]);i+=1
    return [data['corrected_binary'],*result,'-u',str(data['chunk_bp']),'-s']


def source_fragment(name,lengths):
    if name in lengths:return name,0,lengths[name]
    match=re.fullmatch(r'(.*):(\d+)-(\d+)',name)
    if match and match[1] in lengths:
        base=match[1];lo=int(match[2])-1;hi=int(match[3])
        if not 0<=lo<hi<=lengths[base]:raise ValueError('Fragment outside source bounds')
        return base,lo,hi
    if name.startswith('id=') and not name.startswith('id=_MINIGRAPH_|'):
        raise ValueError('Unknown biological source in production partition: '+name)
    return None


def complement(spans,length):
    result=[];end=0
    for lo,hi in sorted(spans):
        if lo<end:raise ValueError('Overlapping or duplicate source assignments')
        if not 0<=lo<hi<=length:raise ValueError('Assigned region outside source bounds')
        if lo>end:result.append((end,lo))
        end=hi
    if end<length:result.append((end,length))
    return result


def audit_split(data,prefix,paf,ident):
    sources=read_table(data['source_lengths']);lengths={r['query']:int(r['source_bp']) for r in sources}
    by=defaultdict(list);allowed=set();region_rows=[]
    for path in Path(prefix).parent.glob(Path(prefix).name+'*.fa_contigs'):
        chrom=path.name[len(Path(prefix).name):-len('.fa_contigs')]
        for name in path.read_text().splitlines():
            parsed=source_fragment(name,lengths)
            if parsed is None:continue
            query,lo,hi=parsed;by[query].append((lo,hi,chrom));allowed.add((query,lo,hi,chrom))
            region_rows.append(dict(query=query,chromosome=chrom,start_0based=lo,end_exclusive=hi))
    summaries=[];excluded=[]
    for source in sources:
        query=source['query'];length=lengths[query];spans=by[query]
        gaps=complement([(lo,hi) for lo,hi,c in spans],length)
        named=union_bp([(lo,hi) for lo,hi,c in spans if c!=data['ambiguous_name']])
        ambiguous=union_bp([(lo,hi) for lo,hi,c in spans if c==data['ambiguous_name']])
        for lo,hi in gaps:excluded.append(dict(query=query,reason='outside_partition',start_0based=lo,end_exclusive=hi))
        for lo,hi,c in spans:
            if c==data['ambiguous_name']:excluded.append(dict(query=query,reason='ambiguous_not_in_chromfile',start_0based=lo,end_exclusive=hi))
        summaries.append(dict(**source,assigned_bp=named,ambiguous_bp=ambiguous,
            unreported_bp=sum(hi-lo for lo,hi in gaps),overlap_bp=0,status='PASS'))
    # PAF query coordinates become fragment-relative after rgfa-split.
    for path in Path(prefix).parent.glob(Path(prefix).name+'*.paf'):
        chrom=path.name[len(Path(prefix).name):-len('.paf')]
        with path.open() as f:
            for line in f:
                fields=line.rstrip().split('\t');parsed=source_fragment(fields[0],lengths)
                if parsed is None:continue
                query,lo,hi=parsed
                if (query,lo,hi,chrom) not in allowed or int(fields[1])!=hi-lo or not 0<=int(fields[2])<=int(fields[3])<=hi-lo:
                    raise ValueError('Output mapping outside its assigned source fragment')
    out=Path(data['evidence'])
    write_table(out/(ident+'.sources.tsv'),summaries,list(summaries[0]))
    for suffix,rows,columns in [('regions',region_rows,['query','chromosome','start_0based','end_exclusive']),
                                ('excluded',excluded,['query','reason','start_0based','end_exclusive'])]:
        plain=out/(ident+'.'+suffix+'.tsv');write_table(plain,rows,columns)
        with plain.open('rb') as src,gzip.open(str(plain)+'.gz','wb') as dst:shutil.copyfileobj(src,dst)
        plain.unlink()
    return dict(sources=len(sources),assigned_bp=sum(r['assigned_bp'] for r in summaries),
        ambiguous_bp=sum(r['ambiguous_bp'] for r in summaries),unreported_bp=sum(r['unreported_bp'] for r in summaries),
        coordinate_integrity='PASS',input_paf_sha256=sha256(paf))


def regional_split(policy,args):
    data=load_policy(policy);cmd=regional_command(data,args)
    if '-b' not in cmd or '-p' not in cmd:raise ValueError('Unsupported regional split invocation')
    # The production first-pass policy is regional; no change to source assembly coordinates.
    subprocess.run(cmd,check=True)
    ident='split_'+uuid.uuid4().hex
    stats=audit_split(data,cmd[cmd.index('-b')+1],cmd[cmd.index('-p')+1],ident)
    dump(Path(data['evidence'])/(ident+'.json'),dict(kind='split',policy=POLICY,chunk_bp=data['chunk_bp'],command=cmd,**stats))


def summarize(policy):
    data=load_policy(policy);out=Path(data['evidence'])
    events=[json.loads(p.read_text()) for p in out.glob('*.json') if p.name.startswith(('filter_','split_'))]
    filters=[r for r in events if r.get('kind')=='filter'];splits=[r for r in events if r.get('kind')=='split']
    if not filters or not splits:raise ValueError('Production wrappers were not invoked; graph policy not verified')
    if len(splits)!=1:raise ValueError('Unexpected multiple regional split passes; source accounting requires explicit reconciliation')
    stats=splits[0]
    text=['# Graph input sequence accounting','',f'Policy: `{POLICY}`. Regional chunk: {data["chunk_bp"]:,} bp.','',
        'Deletion filtering runs on separate query/chromosome lanes. Regional assignment uses the pinned grouping correction. '
        'Source coordinates and input assembly sequences are unchanged.','',
        '| Category | Bases |','|---|---:|',f'| Assigned to chromosome partitions | {stats["assigned_bp"]:,} |',
        f'| Ambiguous, omitted from chromosome subproblems | {stats["ambiguous_bp"]:,} |',
        f'| Outside reported partitions, including softclipped ends | {stats["unreported_bp"]:,} |','',
        'Assigned/ambiguous regions are in bounds and non-overlapping. Every emitted biological mapping belongs to its assigned fragment. '
        'Excluded intervals remain in the source assemblies and are recorded in `split_*.excluded.tsv.gz`; they are not silently assigned '
        'to a chromosome or included in an other bin. FULL and CLIP representation must be evaluated separately downstream.','']
    (out/'sequence_accounting.md').write_text('\n'.join(text),encoding='utf-8')
    dump(out/'sequence_accounting.json',dict(stats,policy=POLICY,filter_invocations=len(filters),split_invocations=len(splits)))


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='mode',required=True)
    b=s.add_parser('build');b.add_argument('--source',required=True);b.add_argument('--out',required=True)
    b=s.add_parser('prepare')
    for k in ('tools','seqfile','runtime','out'):b.add_argument('--'+k,required=True)
    b.add_argument('--chunk',type=int,default=25000)
    for mode in ('filter','split'):
        b=s.add_parser(mode);b.add_argument('--policy',required=True);b.add_argument('args',nargs=argparse.REMAINDER)
    b=s.add_parser('summarize');b.add_argument('--policy',required=True)
    a=vars(p.parse_args());mode=a.pop('mode')
    if mode in ('filter','split'):
        args=a.pop('args');a['args']=args[1:] if args[:1]==['--'] else args
    {'build':build,'prepare':prepare,'filter':guarded_filter,'split':regional_split,'summarize':summarize}[mode](**a)


if __name__=='__main__':main()
