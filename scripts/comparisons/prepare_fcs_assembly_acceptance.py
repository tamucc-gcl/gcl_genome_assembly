#!/usr/bin/env python3
"""Bind the bounded FCS integration check to retained real assembly inputs."""
import argparse
import csv
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


def configure(args):
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
        source=screen/tokens[tokens.index('--fasta')+1]
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


def verify(args):
    manifest=json.loads((args.out/'inputs.json').read_text())
    assembly=manifest['assembly']
    for item in manifest['inputs'].values():
        if digest(item['path'])!=item['sha256']:raise ValueError('Input FASTA changed during acceptance: '+item['path'])
    summary={}
    for stage in ('contig','final'):
        base=args.out/'results'/stage/'decontam'
        report=json.loads((base/'report/report.json').read_text())
        if (not report['fcs_enabled'] or report['blobtools_enabled'] or
            report['plot_classification_source']!='FCS-GX' or
            len(report['assemblies'])!=1 or report['assemblies'][0]['assembly']!=assembly):
            raise ValueError('Incorrect report scope: '+stage)
        metadata=json.loads((base/'fcs'/f'{assembly}.fcs_metrics.json').read_text())
        if metadata['assembly']!=assembly or metadata['coordinate_stage']!=stage or metadata['contigs']<=0:
            raise ValueError('Incorrect plot provenance: '+stage)
        if stage=='final' and metadata['source_sha256']!=manifest['inputs']['final_fasta']['sha256']:
            raise ValueError('Final plot did not use the exact retained finalized FASTA')
        if stage=='contig' and metadata['source_sha256']!=digest(base/'fcs'/f'{assembly}.adaptor_cleaned.fasta'):
            raise ValueError('Contig plot did not use the exact adaptor-cleaned GX input')
        if metadata['taxonomy_sha256']!=digest(base/'fcs'/f'{assembly}.taxonomy.rpt'):
            raise ValueError('Plot taxonomy provenance differs from the published GX evidence')
        with (base/'fcs'/f'{assembly}.fcs_metrics.tsv').open() as handle:
            count=sum(1 for _ in csv.DictReader(handle,delimiter='\t'))
        if count!=metadata['contigs']:raise ValueError('Incomplete plot evidence: '+stage)
        for name,signature in [('cohort_gc_coverage.png',b'\x89PNG\r\n\x1a\n'),('cohort_gc_coverage.pdf',b'%PDF-')]:
            with (base/'report'/name).open('rb') as handle:
                if handle.read(len(signature))!=signature:raise ValueError('Missing rendered plot: '+stage+'/'+name)
        summary[stage]=dict(contigs=count,source_sha256=metadata['source_sha256'])
    result=dict(status='PASS',assembly=assembly,final_fasta_unchanged=True,stages=summary,
                scope='Integration and evidence checks; FCS findings still require review; figures require visual review')
    (args.out/'acceptance.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=['configure','verify'])
    p.add_argument('--project',type=Path,default=Path.cwd())
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--trace',type=Path)
    p.add_argument('--assembly',default='Sde-CTlk_104_hap1')
    p.add_argument('--taxid',type=int,default=373251)
    p.add_argument('--gxdb-dir',type=Path,default=Path('/work/birdlab/databases/fcs-gx'))
    for key in ('source-fasta','final-fasta','hifi-fastq'):p.add_argument('--'+key,type=Path)
    args=p.parse_args()
    if args.command=='configure':configure(args)
    else:verify(args)
