#!/usr/bin/env python3
"""Bind the bounded FCS integration check to retained real assembly inputs."""
import argparse
import csv
import copy
import re
import hashlib
import json
from pathlib import Path
import shlex


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def task(project,trace,process,identifier):
    with trace.open() as handle:
        hashes={r['hash'] for r in csv.DictReader(handle,delimiter='\t')
                if r['name'].endswith(f'{process} ({identifier})')
                and r['status'] in ('CACHED','COMPLETED') and r['exit']=='0'}
    if len(hashes)!=1:raise ValueError(f'Expected one successful {process} ({identifier}) in {trace}, found {hashes}')
    prefix,short=hashes.pop().split('/')
    matches=list((project/'work'/prefix).glob(short+'*'))
    if len(matches)!=1:raise ValueError(f'Expected one retained task for {prefix}/{short}')
    return matches[0]


def staged_input(directory,name):
    """Resolve a scratch-only input using Nextflow's retained staging command.

    Scratch tasks need not retain input links in their work directory. Read the
    exact literal `ln -s SOURCE NAME` binding; never execute the task wrapper or
    search for a different FASTA with a matching basename.
    """
    requested=Path(name)
    path=directory/requested
    if path.is_file():return path.resolve(strict=True)
    candidates=set()
    wrapper=directory/'.command.run'
    if wrapper.is_file():
        for line in wrapper.read_text().splitlines():
            try:tokens=shlex.split(line)
            except ValueError:continue
            if (len(tokens)>=4 and tokens[0]=='ln' and '-s' in tokens[1:-2]
                    and Path(tokens[-1])==requested):
                source=Path(tokens[-2])
                if source.is_absolute() and source.is_file():
                    candidates.add(source.resolve(strict=True))
    if len(candidates)!=1:
        raise ValueError(f'Cannot resolve retained input {path}: the work-directory link is absent '
                         f'and .command.run has {len(candidates)} usable exact staging bindings. '
                         'Supply --source-fasta with the retained pre-FCS FASTA; '
                         'do not substitute a post-FCS cleaned assembly.')
    return candidates.pop()


def configure(args):
    if getattr(args,"cohort",False):return configure_cohort(args)
    project=args.project.resolve(strict=True)
    trace=(args.trace or project/'genome_assembly/pipeline/pipeline_trace.txt').resolve(strict=True)
    sample=args.assembly.rsplit('_hap',1)[0]
    if args.source_fasta:
        source=args.source_fasta
    else:
        screen=task(project,trace,'FCS_GX_SCREEN',args.assembly)
        # Read the exact original GX input from the successful task command.
        tokens=shlex.split((screen/'.command.sh').read_text().replace('\\\n',' '))
        if tokens.count('--fasta')!=1:raise ValueError('Cannot identify exact GX input; supply --source-fasta')
        source=staged_input(screen,tokens[tokens.index('--fasta')+1])
    final=args.final_fasta or task(project,trace,'FINALIZE_ASSEMBLY',args.assembly)/(args.assembly+'.fasta')
    reads=args.hifi_fastq or task(project,trace,'BAM_TO_FASTQ',sample)/(sample+'.fastq.gz')
    paths={k:Path(v).resolve(strict=True) for k,v in
           dict(source_fasta=source,final_fasta=final,hifi_fastq=reads).items()}
    for key,path in paths.items():
        if not path.is_file() or not path.stat().st_size:raise ValueError('Missing/empty '+key+': '+str(path))
    db=args.gxdb_dir.resolve(strict=True)
    if len(list(db.glob('*.gxi')))!=1:raise ValueError('Expected exactly one GX database prefix in '+str(db))
    args.out.mkdir(parents=True,exist_ok=True)
    params=dict(assembly=args.assembly,sample=sample,taxid=args.taxid,gxdb_dir=str(db),
                run_contamination_evidence=False,publish_dir_mode='copy',
                **{k:str(v) for k,v in paths.items()})
    (args.out/'params.json').write_text(json.dumps(params,indent=2)+'\n')
    manifest=dict(assembly=args.assembly,trace=str(trace),inputs={k:dict(path=str(v),sha256=digest(v))
                  for k,v in paths.items() if k!='hifi_fastq'},read_path=str(paths['hifi_fastq']),
                  scope='Complete single haplotype; no sequence sampling; no BlobTools searches')
    (args.out/'inputs.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(params,indent=2),flush=True)


def configure_cohort(args):
    project=args.project.resolve(strict=True)
    trace=(args.trace or project/'genome_assembly/pipeline/pipeline_trace.txt').resolve(strict=True)
    if args.source_fasta or args.final_fasta or args.hifi_fastq:
        raise ValueError('Single-assembly overrides cannot be used with --cohort')
    with trace.open() as f:
        trace_rows=list(csv.DictReader(f,delimiter='\t'))
    ids=sorted({re.search(r'FINALIZE_ASSEMBLY \(([^)]+)\)$',r['name'])[1]
                    for r in trace_rows
                    if re.search(r'FINALIZE_ASSEMBLY \(([^)]+)\)$',r['name'])
                    and r['status'] in ('CACHED','COMPLETED') and r['exit']=='0'})
    if not ids:raise ValueError('No finalized assemblies found in production trace')
    screened={re.search(r'FCS_GX_SCREEN \(([^)]+)\)$',r['name'])[1] for r in trace_rows
              if re.search(r'FCS_GX_SCREEN \(([^)]+)\)$',r['name'])
              and r['status'] in ('CACHED','COMPLETED') and r['exit']=='0'}
    if set(ids)!=screened:
        raise ValueError('Initial FCS and final assembly cohort membership differ; refusing a partial cohort')
    manifests=[];inputs=[]
    for identifier in ids:
        child=copy.copy(args);child.cohort=False;child.assembly=identifier
        child.out=args.out/'inputs'/identifier
        configure(child)
        manifests.append(json.loads((child.out/'inputs.json').read_text()))
        inputs.append(json.loads((child.out/'params.json').read_text()))
    # One matching HiFi source per individual; reject conflicting retained inputs.
    read_sources={}
    for item in inputs:
        if item['sample'] in read_sources and read_sources[item['sample']]!=item['hifi_fastq']:
            raise ValueError('Conflicting read inputs for '+item['sample'])
        read_sources[item['sample']]=item['hifi_fastq']
    params=dict(acceptance_assemblies=inputs,gxdb_dir=inputs[0]['gxdb_dir'],
                run_contamination_evidence=False,publish_dir_mode='copy')
    (args.out/'params.json').write_text(json.dumps(params,indent=2)+'\n')
    (args.out/'inputs.json').write_text(json.dumps(dict(assemblies=manifests,scope='All finalized production-trace assemblies; no sampling'),indent=2)+'\n')
    print(f'Configured {len(ids)} complete assemblies from {len(read_sources)} individuals',flush=True)


def verify(args):
    manifest=json.loads((args.out/'inputs.json').read_text())
    manifests=manifest.get('assemblies',[manifest])
    expected={m['assembly'] for m in manifests}
    for item in manifests:
        for value in item['inputs'].values():
            if digest(value['path'])!=value['sha256']:raise ValueError('Input FASTA changed: '+value['path'])
    summary={}
    for stage in ('contig','final'):
        base=args.out/'results'/stage/'decontam'
        report=json.loads((base/'report/report.json').read_text())
        if (not report['fcs_enabled'] or report['blobtools_enabled'] or
            report['plot_classification_source']!='FCS-GX' or
            {r['assembly'] for r in report['assemblies']}!=expected):
            raise ValueError('Incorrect report scope: '+stage)
        summary[stage]={}
        for item in manifests:
            assembly=item['assembly']
            metadata=json.loads((base/'fcs'/f'{assembly}.fcs_metrics.json').read_text())
            if metadata['assembly']!=assembly or metadata['coordinate_stage']!=stage or metadata['contigs']<=0:
                raise ValueError('Incorrect plot provenance: '+stage+'/'+assembly)
            source=item['inputs']['final_fasta']['path'] if stage=='final' else base/'fcs'/f'{assembly}.adaptor_cleaned.fasta'
            if metadata['source_sha256']!=digest(source):raise ValueError('Plot source differs: '+stage+'/'+assembly)
            if metadata['taxonomy_sha256']!=digest(base/'fcs'/f'{assembly}.taxonomy.rpt'):
                raise ValueError('Taxonomy provenance differs: '+stage+'/'+assembly)
            with (base/'fcs'/f'{assembly}.fcs_metrics.tsv').open() as f:count=sum(1 for _ in csv.DictReader(f,delimiter='\t'))
            if count!=metadata['contigs']:raise ValueError('Incomplete plot evidence')
            summary[stage][assembly]=dict(contigs=count,source_sha256=metadata['source_sha256'])
        for name,signature in [('cohort_gc_coverage.png',b'\x89PNG\r\n\x1a\n'),('cohort_gc_coverage.pdf',b'%PDF-'),('cohort_fcs_composition.png',b'\x89PNG\r\n\x1a\n')]:
            with (base/'report'/name).open('rb') as f:
                if f.read(len(signature))!=signature:raise ValueError('Missing rendered plot: '+stage+'/'+name)
        with (base/'report/cohort_fcs_sequence_evidence.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
        categories=('primary','repeat_expected','low_expected','repeat_other','inconclusive','taxonomy_flag','unreported')
        for row in rows:
            if sum(int(row['cat_'+c+'_bp']) for c in categories)!=int(row['length_bp']):raise ValueError('Incomplete category partition')
        if len(rows)!=sum(x['contigs'] for x in summary[stage].values()):raise ValueError('Incorrect cohort evidence census')
        for index in range(1,len(expected)+1):
            for suffix in ('gc_coverage','composition','diagnostics'):
                if not (base/'report'/f'assembly_{index:04d}_{suffix}.png').is_file():raise ValueError('Missing assembly detail figure')
    if not (args.out/'results/overview/decontam/cohort_fcs_composition.png').is_file():raise ValueError('Missing combined stage overview')
    result=dict(status='PASS',assemblies=sorted(expected),final_fastas_unchanged=True,stages=summary,
                scope='Integration and complete category partitions; review FCS findings and rendered figures separately')
    (args.out/'acceptance.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=['configure','verify'])
    p.add_argument('--project',type=Path,default=Path.cwd())
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--trace',type=Path)
    p.add_argument('--cohort',action='store_true')
    p.add_argument('--assembly',default='Sde-CTlk_104_hap1')
    p.add_argument('--taxid',type=int,default=373251)
    p.add_argument('--gxdb-dir',type=Path,default=Path('/work/birdlab/databases/fcs-gx'))
    for key in ('source-fasta','final-fasta','hifi-fastq'):p.add_argument('--'+key,type=Path)
    args=p.parse_args()
    if args.command=='configure':configure(args)
    else:verify(args)
